// Copyright (C) 2012 Rhys Ulerich
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "ar.hpp"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include <octave/oct.h>
#include <octave/oct-map.h>
#include <octave/ov-struct.h>
#include <octave/Cell.h>

/** @file
 * A GNU Octave function estimating the best AR(p) model given signal input.
 * Compare \ref arsel.cpp.
 */

// FIXME Add window_T0-like option per arsel.cpp

// Compile-time defaults in the code also appearing in the help message
#define DEFAULT_SUBMEAN   true
#define DEFAULT_ABSRHO    true
#define DEFAULT_CRITERION "CIC"
#define DEFAULT_MINORDER  0
#define DEFAULT_MAXORDER  512
#define STRINGIFY(x) STRINGIFY_HELPER(x)
#define STRINGIFY_HELPER(x) #x

DEFUN_DLD(arsel, args, nargout,
"    M = arsel (data, submean, absrho, criterion, minorder, maxorder)\n"
"    Automatically fit autoregressive models to input signals.\n"
"\n"
"    Use ar::burg_method and ar::best_model to fit an autoregressive process\n"
"    for signals contained in the rows of matrix data.  Sample means will\n"
"    be subtracted whenever submean is true.  Model orders minorder through\n"
"    min(columns(data), maxorder) will be considered.  A structure is\n"
"    returned where each field either contains a result indexable by the\n"
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
"    in field 'N'.  The filter()-ready process parameters are returned\n"
"    in field 'AR', the sample mean in 'mu', and the innovation variance\n"
"    \\sigma^2_\\epsilon in 'sigma2eps'.  The process output variance\n"
"    \\sigma^2_\\x and process gain are returned in fields 'sigma2x' and\n"
"    'gain', respectively.  Autocorrelations for lags zero through the\n"
"    model order, inclusive, are returned in field 'autocor'.  The raw\n"
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
"    For example, given a *row-vector* of samples 'd', one can fit a\n"
"    process and then simulate a sample realization of length M using\n"
"\n"
"        a = arsel(d);\n"
"        x = a.mu(1) + filter([1], a.AR{1}, sqrt(a.sigma2eps(1)).*randn(1,M));\n"
"\n"
"    Continuing the example, if the Octave signal package has been loaded,\n"
"    one can compute the model autocorrelation function at lags 0:M using\n"
"\n"
"        si  = filtic([1], a.AR{1}, a.autocor{1});\n"
"        rho = [1, filter([1], a.AR{1}, zeros(1,M), si)];\n"
"\n"
"    When omitted, submean defaults to " STRINGIFY(DEFAULT_SUBMEAN) ".\n"
"    When omitted, absrho defaults to " STRINGIFY(DEFAULT_ABSRHO) ".\n"
"    When omitted, criterion defaults to " STRINGIFY(DEFAULT_CRITERION) ".\n"
"    When omitted, minorder defaults to " STRINGIFY(DEFAULT_MINORDER) ".\n"
"    When omitted, maxorder defaults to " STRINGIFY(DEFAULT_MAXORDER) ".\n"
)
{
    using std::size_t;
    using std::string;
    typedef Matrix::element_type element_type;
    typedef std::vector<element_type> vector;

    size_t maxorder  = DEFAULT_MAXORDER;
    size_t minorder  = DEFAULT_MINORDER;
    string criterion = DEFAULT_CRITERION;
    bool   absrho    = DEFAULT_ABSRHO;
    bool   submean   = DEFAULT_SUBMEAN;
    Matrix data;
    switch (args.length())
    {
        case 6: maxorder  = args(5).ulong_value();
        case 5: minorder  = args(4).ulong_value();
        case 4: criterion = args(3).string_value();
        case 3: absrho    = args(2).bool_value();
        case 2: submean   = args(1).bool_value();
        case 1: data      = args(0).matrix_value();
                if (!error_state) break;
        default:
            error("Invalid call to arsel.  See 'help arsel' for usage.");
            return octave_value();
        case 0:
            print_usage();
            return octave_value();
    }

    // Lookup the desired model selection criterion
    typedef ar::best_model_function<
                ar::Burg,octave_idx_type,octave_idx_type,vector
            > best_model_function;
    const best_model_function::type best_model
            = best_model_function::lookup(criterion, submean);
    if (!best_model)
    {
        error("Unknown model selection criterion provided to arsel.");
        return octave_value();
    }

    const octave_idx_type M = data.rows();  // Number of signals
    const octave_idx_type N = data.cols();  // Samples per signal

    // Prepare per-signal storage locations to return to caller
    Cell         _AR       (dim_vector(M,1));
    Cell         _autocor  (dim_vector(M,1));
    ColumnVector _eff_N    (M);
    ColumnVector _eff_var  (M);
    ColumnVector _gain     (M);
    ColumnVector _mu       (M);
    ColumnVector _mu_sigma (M);
    ColumnVector _sigma2eps(M);
    ColumnVector _sigma2x  (M);
    ColumnVector _T0       (M);

    // Prepare vectors to capture burg_method() output
    vector params, sigma2e, gain, autocor;
    params .reserve(maxorder*(maxorder + 1)/2);
    sigma2e.reserve(maxorder + 1);
    gain   .reserve(maxorder + 1);
    autocor.reserve(maxorder + 1);

    // Prepare repeatedly-used working storage for burg_method()
    vector f, b, Ak, ac;

    // Process each signal in turn...
    for (octave_idx_type i = 0; i < M; ++i)
    {
        // Use burg_method to estimate a hierarchy of AR models from input data
        params .clear();
        sigma2e.clear();
        gain   .clear();
        autocor.clear();
        ar::strided_adaptor<const element_type*> signal_begin(&data(i,0), M);
        ar::strided_adaptor<const element_type*> signal_end  (&data(i,N), M);
        ar::burg_method(signal_begin, signal_end, _mu(i), maxorder,
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
            error(e.what());
            return octave_value();
        }

        // Compute decorrelation time from the estimated autocorrelation model
        ar::predictor<element_type> p = ar::autocorrelation(
                params.begin(), params.end(), gain[0], autocor.begin());
        _T0(i) = ar::decorrelation_time(N, p, absrho);

        // Filter()-ready process parameters in field 'AR' with leading one
        {
            RowVector t(params.size() + 1);
            t(0) = 1;
            std::copy(params.begin(), params.end(), t.fortran_vec() + 1);
            _AR(i) = t;
        }

        // Field 'sigma2eps'
        _sigma2eps(i) = sigma2e[0];

        // Field 'gain'
        _gain(i) = gain[0];

        // Field 'sigma2x'
        _sigma2x(i) = gain[0]*sigma2e[0];

        // Field 'autocor'
        {
            RowVector t(autocor.size());
            std::copy(autocor.begin(), autocor.end(), t.fortran_vec());
            _autocor(i) = t;
        }

        // Field 'eff_var'
        // Unbiased effective variance expression from [Trenberth1984]
        _eff_var(i) = (N*gain[0]*sigma2e[0]) / (N - _T0(i));

        // Field 'eff_N'
        _eff_N(i) = N / _T0(i);

        // Field 'mu_sigma'
        // Variance of the sample mean using effective quantities
        _mu_sigma(i) = std::sqrt(_eff_var(i) / _eff_N(i));

        // Permit user to interrupt the computations at this time
        OCTAVE_QUIT;
    }

    // Provide no results whenever an error was detected
    if (error_state)
    {
        warning("arsel: error detected; no results returned");
        return octave_value_list();
    }

    // Build map containing return fields
    Octave_map retval;
    retval.assign("AR",        octave_value(_AR));
    retval.assign("absrho",    octave_value(absrho));
    retval.assign("autocor",   octave_value(_autocor));
    retval.assign("criterion", octave_value(criterion));
    retval.assign("data",      data);
    retval.assign("eff_N",     _eff_N);
    retval.assign("eff_var",   _eff_var);
    retval.assign("gain",      _gain);
    retval.assign("maxorder",  octave_value(maxorder));
    retval.assign("minorder",  octave_value(minorder));
    retval.assign("mu",        _mu);
    retval.assign("mu_sigma",  _mu_sigma);
    retval.assign("N",         octave_value(N));
    retval.assign("sigma2eps", _sigma2eps);
    retval.assign("sigma2x",   _sigma2x);
    retval.assign("submean",   octave_value(submean));
    retval.assign("T0",        _T0);

    return octave_value_list(retval);
}
