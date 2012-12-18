// Copyright (C) 2012 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Python.h>
#include <numpy/arrayobject.h>

#include <cstdlib>
#include <string>
#include <vector>

#include "ar.hpp"

/** @file
 * A Python extension module for exposing autoregressive model utilities.
 */

// FIXME Add window_T0-like option per arsel.cpp

// Compile-time defaults in the code also appearing in arsel docstring
#define DEFAULT_SUBMEAN   true
#define DEFAULT_ABSRHO    true
#define DEFAULT_CRITERION "CIC"
#define DEFAULT_MINORDER  0
#define DEFAULT_MAXORDER  512
#define STRINGIFY(x) STRINGIFY_HELPER(x)
#define STRINGIFY_HELPER(x) #x

// TODO The returned values are no longer keys but attributes
static const char ar_arsel_docstring[] =
"    Usage: M = arsel (data, submean, absrho, criterion, minorder, maxorder)\n"
"\n"
"    Automatically fit autoregressive models to input signals.\n"
"\n"
"    Use ar::burg_method and ar::best_model to fit an autoregressive process\n"
"    for signals contained in the rows of matrix data.  Sample means will\n"
"    be subtracted whenever submean is true.  Model orders minorder through\n"
"    min(columns(data), maxorder) will be considered.  A dictionary is\n"
"    returned where each key either contains a result indexable by the\n"
"    signal number (i.e. the row indices of input matrix data) or it contains\n"
"    a single scalar applicable to all signals.\n"
"\n"
"    The model order will be selected using the specified criterion.\n"
"    Criteria are specified using the following abbreviations:\n"
"        AIC  - Akaike information criterion\n"
"        AICC - asymptotically-corrected Akaike information criterion\n"
"        BIC  - consistent criterion BIC\n"
"        CIC  - combined information criterion\n"
"        FIC  - finite information criterion\n"
"        FSIC - finite sample information criterion\n"
"        GIC  - generalized information criterion\n"
"        MCC  - minimally consistent criterion\n"
"\n"
"    The number of samples in data (i.e. the number of rows) is returned\n"
"    in key 'N'.  The filter()-ready process parameters are returned\n"
"    in key 'AR', the sample mean in 'mu', and the innovation variance\n"
"    \\sigma^2_\\epsilon in 'sigma2eps'.  The process output variance\n"
"    \\sigma^2_\\x and process gain are returned in keys 'sigma2x' and\n"
"    'gain', respectively.  Autocorrelations for lags zero through the\n"
"    model order, inclusive, are returned in key 'autocor'.  The raw\n"
"    signals are made available for later use in field 'data'.\n"
"\n"
"    Given the observed autocorrelation structure, a decorrelation time\n"
"    'T0' is computed by ar::decorrelation_time and used to estimate\n"
"    the effective signal variance 'eff_var'.  The number of effectively\n"
"    independent samples is returned in 'eff_N'.  These effective values are\n"
"    combined to estimate the sampling error (i.e. the standard deviation\n"
"    of the sample mean) as field 'mu_sigma'.  The absolute value of the\n"
"    autocorrelation function will be used in computing the decorrelation\n"
"    times whenever absrho is true.\n"
"\n"
"    For example, given a sequence or *row-vector* of samples 'd', one can fit\n"
"    a process, obtain its poles, and simulate a realization of length M using\n"
"\n"
"        from ar import arsel\n"
"        from numpy import reciprocal, roots\n"
"        from numpy.random import randn\n"
"        from scipy.signal import lfilter\n"
"        a = arsel(d)\n"
"        p = [reciprocal(roots(model[::-1])) for model in a.AR]\n"
"        x = a.mu[0] + lfilter([1], a.AR[0], sqrt(a.sigma2eps[0])*randn(M))\n"
"\n"
"    Continuing the example, one can compute the model autocorrelation\n"
"    function at lags 0:M using\n"
"\n"
"        from numpy import array, ones\n"
"        from scipy.signal import lfiltic\n"
"        zi      = lfiltic([1], a.AR[0], a.autocor[0])\n"
"        rho     = ones(M+1)\n"
"        rho[1:] = lfilter([1], a.AR[0], zeros(M), zi=zi)[0]\n"
"\n"
"    When omitted, submean defaults to " STRINGIFY(DEFAULT_SUBMEAN) ".\n"
"    When omitted, absrho defaults to " STRINGIFY(DEFAULT_ABSRHO) ".\n"
"    When omitted, criterion defaults to " STRINGIFY(DEFAULT_CRITERION) ".\n"
"    When omitted, minorder defaults to " STRINGIFY(DEFAULT_MINORDER) ".\n"
"    When omitted, maxorder defaults to " STRINGIFY(DEFAULT_MAXORDER) ".\n"
;

// Return type for ar_arsel initialized within initar
static PyTypeObject *ar_ArselType = NULL;

extern "C" PyObject *ar_arsel(PyObject *self, PyObject *args)
{
    PyObject *ret_args = NULL, *ret = NULL;

    // Sanity check that initar could build ar_ArselType
    if (!PyType_Check(ar_ArselType)) {
        PyErr_SetString(PyExc_RuntimeError,
            "PyType_Check(ar_ArselType) failed");
        return NULL;
    }

    // Prepare argument storage with default values
    PyObject   *data_obj  = NULL;
    int         submean   = DEFAULT_SUBMEAN;
    int         absrho    = DEFAULT_ABSRHO;
    const char *criterion = DEFAULT_CRITERION;
    std::size_t minorder  = DEFAULT_MINORDER;
    std::size_t maxorder  = DEFAULT_MAXORDER;

    // Parse input tuple with second and subsequent arguments optional
    {
        unsigned long ul_minorder = minorder, ul_maxorder = maxorder;
        if (!PyArg_ParseTuple(args, "O|iiskk", &data_obj, &submean, &absrho,
                                               &criterion,
                                               &ul_minorder, &ul_maxorder)) {
            return NULL;
        }
        minorder = ul_minorder;
        maxorder = ul_maxorder;
    }

    // Lookup the desired model selection criterion
    typedef ar::best_model_function<
                ar::Burg,npy_intp,npy_intp,std::vector<double>
            > best_model_function;
    const best_model_function::type best_model
            = best_model_function::lookup(std::string(criterion), submean);
    if (!best_model) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unknown model selection criterion provided to arsel.");
        return NULL;
    }

    // Incoming data may be noncontiguous but should otherwise be well-behaved
    // On success, 'data' is returned so Py_DECREF is applied only on failure
    PyObject *data = PyArray_FROMANY(data_obj, NPY_DOUBLE, 1, 2,
            NPY_ALIGNED | NPY_ELEMENTSTRIDES | NPY_NOTSWAPPED);
    if (!data) {
        Py_DECREF(data);
        return NULL;
    }

    // Reshape any 1D ndarray into a 2D ndarray organized as a row vector.
    // Permits the remainder of the routine to uniformly worry about 2D only.
    if (1 == ((PyArrayObject *)(data))->nd) {
        npy_intp dim[2] = { 1, PyArray_DIM(data, 0) };
        PyArray_Dims newshape = { dim, sizeof(dim)/sizeof(dim[0]) };
        PyObject *olddata = data;
        data = PyArray_Newshape((PyArrayObject *)olddata,
                                &newshape, NPY_ANYORDER);
        Py_DECREF(olddata);
        if (!data) {
            Py_DECREF(data);
            PyErr_SetString(PyExc_RuntimeError,
                "Unable to reorganize data into matrix-like row vector.");
            return NULL;
        }
    }

    // How many data points are there?
    npy_intp M = PyArray_DIM(data, 0);
    npy_intp N = PyArray_DIM(data, 1);

    // Prepare per-signal storage locations to return to caller
    // TODO Ensure these invocations all worked as expected
    PyObject *_AR        = PyList_New(M);
    PyObject *_autocor   = PyList_New(M);
    PyObject *_eff_N     = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_eff_var   = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_gain      = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_mu        = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_mu_sigma  = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_sigma2eps = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_sigma2x   = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);
    PyObject *_T0        = PyArray_ZEROS(1, &M, NPY_DOUBLE, 0);

    // Prepare vectors to capture burg_method() output
    std::vector<double> params, sigma2e, gain, autocor;
    params .reserve(maxorder*(maxorder + 1)/2);
    sigma2e.reserve(maxorder + 1);
    gain   .reserve(maxorder + 1);
    autocor.reserve(maxorder + 1);

    // Prepare repeatedly-used working storage for burg_method()
    std::vector<double> f, b, Ak, ac;

    // Process each signal in turn...
    for (npy_intp i = 0; i < M; ++i)
    {
        // Use burg_method to estimate a hierarchy of AR models from input data
        params .clear();
        sigma2e.clear();
        gain   .clear();
        autocor.clear();
        ar::strided_adaptor<const double*> signal_begin(
                (const double*) PyArray_GETPTR2(data, i, 0),
                PyArray_STRIDES(data)[1] / sizeof(double));
        ar::strided_adaptor<const double*> signal_end  (
                (const double*) PyArray_GETPTR2(data, i, N),
                PyArray_STRIDES(data)[1] / sizeof(double));
        ar::burg_method(signal_begin, signal_end,
                        *(double*)PyArray_GETPTR1(_mu, i),
                        maxorder,
                        std::back_inserter(params),
                        std::back_inserter(sigma2e),
                        std::back_inserter(gain),
                        std::back_inserter(autocor),
                        submean, /* output hierarchy? */ true, f, b, Ak, ac);

        // Keep only best model per chosen criterion via function pointer
        try {
            best_model(N, minorder, params, sigma2e, gain, autocor);
        }
        catch (std::exception &e)
        {
            PyErr_SetString(PyExc_RuntimeError, e.what());
            goto fail;
        }

        // Compute decorrelation time from the estimated autocorrelation model
        ar::predictor<double> p = ar::autocorrelation(
                params.begin(), params.end(), gain[0], autocor.begin());
        *(double*)PyArray_GETPTR1(_T0, i)
            = ar::decorrelation_time(N, p, absrho);

        // Filter()-ready process parameters in field 'AR' with leading one
        PyObject *_ARi = PyList_New(params.size() + 1);
        PyList_SET_ITEM(_AR, i, _ARi);
        PyList_SET_ITEM(_ARi, 0, PyFloat_FromDouble(1));
        for (std::size_t k = 0; k < params.size(); ++k) {
            PyList_SET_ITEM(_ARi, k+1, PyFloat_FromDouble(params[k]));
        }

        // Field 'sigma2eps'
        *(double*)PyArray_GETPTR1(_sigma2eps, i) = sigma2e[0];

        // Field 'gain'
        *(double*)PyArray_GETPTR1(_gain, i) = gain[0];

        // Field 'sigma2x'
        *(double*)PyArray_GETPTR1(_sigma2x, i) = gain[0]*sigma2e[0];

        // Field 'autocor'
        PyObject *_autocori = PyList_New(autocor.size());
        PyList_SET_ITEM(_autocor, i, _autocori);
        for (std::size_t k = 0; k < autocor.size(); ++k) {
            PyList_SET_ITEM(_autocori, k, PyFloat_FromDouble(autocor[k]));
        }

        // Field 'eff_var'
        // Unbiased effective variance expression from [Trenberth1984]
        *(double*)PyArray_GETPTR1(_eff_var, i)
            = (N*gain[0]*sigma2e[0]) / (N - *(double*)PyArray_GETPTR1(_T0, i));

        // Field 'eff_N'
        *(double*)PyArray_GETPTR1(_eff_N, i)
            = N / *(double*)PyArray_GETPTR1(_T0, i);

        // Field 'mu_sigma'
        // Variance of the sample mean using effective quantities
        *(double*)PyArray_GETPTR1(_mu_sigma, i)
            = std::sqrt(   *(double*)PyArray_GETPTR1(_eff_var, i)
                         / *(double*)PyArray_GETPTR1(_eff_N,   i));
    }

    // Prepare build and return an ar_ArselType via tuple constructor
    // See initar(...) method for the collections.namedtuple-based definition
    ret_args = PyTuple_Pack(17, PyBool_FromLong(absrho),
                                _AR,
                                _autocor,
                                PyString_FromString(criterion),
                                data,
                                _eff_N,
                                _eff_var,
                                _gain,
                                PyInt_FromSize_t(maxorder),
                                PyInt_FromSize_t(minorder),
                                _mu,
                                _mu_sigma,
                                PyInt_FromLong(N),
                                _sigma2eps,
                                _sigma2x,
                                PyBool_FromLong(submean),
                                _T0);
    if (!ret_args) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unable to prepare arguments used to build return value.");
        goto fail;
    }
    ret = PyObject_CallObject((PyObject *)ar_ArselType, ret_args);
    if (!ret) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unable to construct return value.");
        goto fail;
    }

    return ret;

fail:
    Py_XDECREF(ret_args);
    Py_XDECREF(ret);
    Py_XDECREF(data);
    Py_XDECREF(_AR);
    Py_XDECREF(_autocor);
    Py_XDECREF(_eff_N);
    Py_XDECREF(_eff_var);
    Py_XDECREF(_gain);
    Py_XDECREF(_mu);
    Py_XDECREF(_mu_sigma);
    Py_XDECREF(_sigma2eps);
    Py_XDECREF(_sigma2x);
    Py_XDECREF(_T0);
    return NULL;
}

// Specification of methods available in the module
static PyMethodDef ar_methods[] = {
    {"arsel", ar_arsel, METH_VARARGS, ar_arsel_docstring},
    {NULL, NULL, 0, NULL}
};

// Module docstring
static const char ar_docstring[] = "Autoregressive process modeling tools";

extern "C" PyMODINIT_FUNC initar(void)
{
    // Initialize the module, including making NumPy available
    if (!Py_InitModule3("ar", ar_methods, ar_docstring)) return;
    import_array();

    // Use collections.namedtuple() to make type 'Arsel' for ar_arsel use
    // The tuple is ordered alphabetically in case-insensitive manner
    PyObject *modName = PyString_FromString((char *)"collections");
    PyObject *mod     = PyImport_Import(modName);
    PyObject *func    = PyObject_GetAttrString(mod,(char*)"namedtuple");
    PyObject *args    = PyTuple_Pack(2,
                                     PyString_FromString((char *) "Arsel"),
                                     PyString_FromString((char *) " absrho"
                                                                  " AR"
                                                                  " autocor"
                                                                  " criterion"
                                                                  " data"
                                                                  " eff_N"
                                                                  " eff_var"
                                                                  " gain"
                                                                  " maxorder"
                                                                  " minorder"
                                                                  " mu"
                                                                  " mu_sigma"
                                                                  " N"
                                                                  " sigma2eps"
                                                                  " sigma2x"
                                                                  " submean"
                                                                  " T0"));
    ar_ArselType = (PyTypeObject *) PyObject_CallObject(func, args);
    Py_DECREF(args);
    Py_DECREF(func);
    Py_DECREF(mod);
    Py_DECREF(modName);
}
