// Copyright (C) 2012 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Python.h>
#include <numpy/arrayobject.h>
#include "ar.hpp"

// Compile-time defaults in the code also appearing in arsel docstring
#define DEFAULT_SUBMEAN   true
#define DEFAULT_ABSRHO    false
#define DEFAULT_CRITERION "CIC"
#define DEFAULT_MAXORDER  512
#define STRINGIFY(x) STRINGIFY_HELPER(x)
#define STRINGIFY_HELPER(x) #x

static const char ar_arsel_docstring[] =
"    Usage: M = arsel (data, submean, absrho, criterion, maxorder)\n"
"\n"
"    Automatically fit autoregressive models to input signals.\n"
"\n"
"    Use ar::burg_method and ar::best_model to fit an autoregressive process\n"
"    for signals contained in the rows of matrix data.  Sample means will\n"
"    be subtracted whenever submean is true.  Model orders zero through\n"
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
"    When omitted, submean defaults to " STRINGIFY(DEFAULT_SUBMEAN) ".\n"
"    When omitted, absrho defaults to " STRINGIFY(DEFAULT_ABSRHO) ".\n"
"    When omitted, criterion defaults to " STRINGIFY(DEFAULT_CRITERION) ".\n"
"    When omitted, maxorder defaults to " STRINGIFY(DEFAULT_MAXORDER) ".\n"
;

static PyObject *ar_arsel(PyObject *self, PyObject *args)
{
    double m, b;
    PyObject *x_obj, *y_obj, *yerr_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "ddOOO", &m, &b, &x_obj, &y_obj,
                                         &yerr_obj))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *yerr_array = PyArray_FROM_OTF(yerr_obj, NPY_DOUBLE,
                                            NPY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (x_array == NULL || y_array == NULL || yerr_array == NULL) {
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(yerr_array);
        return NULL;
    }

    /* How many data points are there? */
    int N = (int)PyArray_DIM(x_array, 0);

    /* Get pointers to the data as C-types. */
    double *x    = (double*)PyArray_DATA(x_array);
    double *y    = (double*)PyArray_DATA(y_array);
    double *yerr = (double*)PyArray_DATA(yerr_array);

    /* Call the external C function to compute the chi-squared. */
    double value = chi2(m, b, x, y, yerr, N);

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(yerr_array);

    if (value < 0.0) {
        PyErr_SetString(PyExc_RuntimeError,
                    "Chi-squared returned an impossible value.");
        return NULL;
    }

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}

// Specification of methods available in the module
static PyMethodDef ar_methods[] = {
    {"arsel", ar_arsel, METH_VARARGS, ar_arsel_docstring},
    {NULL, NULL, 0, NULL}
};

// Module docstring
static const char ar_docstring[] = "Autoregressive process modeling tools";

// Initialize the module, including making NumPy available
PyMODINIT_FUNC init_ar(void)
{
    if (!Py_InitModule3("ar", ar_methods, ar_docstring)) return;
    import_array();
}
