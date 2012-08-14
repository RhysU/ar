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
#include <iterator>
#include <limits>
#include <list>
#include <vector>

#include <octave/oct.h>
#include <octave/oct-map.h>
#include <octave/ov-struct.h>
#include <octave/Cell.h>

/** @file
 * A GNU Octave function estimating the best AR(p) model given signal input.
 * Compare \ref arsel.cpp.
 */

// Compile-time defaults in the code also appearing in the help message
#define DEFAULT_SUBMEAN  true
#define DEFAULT_ABSRHO   true
#define DEFAULT_MAXORDER 512
#define STRINGIFY(x) STRINGIFY_HELPER(x)
#define STRINGIFY_HELPER(x) #x

DEFUN_DLD(
    arsel, args, nargout,
    "\tM = arsel (data, submean, absrho, maxorder)\n"
    "\tAutomatically fit autoregressive models to input signals.\n"
    "\t\n"
    "\tUse ar::burg_method and ar::best_model<CIC<Burg<MeahHandling> >\n"
    "\tto fit an autoregressive process for signals contained in the rows\n"
    "\tof matrix data.  Sample means will be subtracted whenever submean\n"
    "\tis true.  Model orders zero through min(columns(data), maxorder)\n"
    "\twill be considered.  A structure is returned where each field\n"
    "\teither contains a result indexable by the signal number (i.e. the\n"
    "\trow indices of input matrix data) or it contains a single scalar\n"
    "\tapplicable to all signals.\n"
    "\t\n"
    "\tThe number of samples in data (i.e. the number of rows) is returned\n"
    "\tin field 'N'.  The filter()-ready process parameters are returned\n"
    "\tin field 'AR', the sample mean in 'mu', and the innovation variance\n"
    "\t\\sigma^2_\\epsilon in 'sigma2eps'.  The process gains are returned\n"
    "\tin 'gain' and the autocorrelation boundary conditions in 'autocor'\n"
    "\tfor lags zero through the model order, inclusive.  The raw signals\n"
    "\tare made available for later use in 'data'.\n"
    "\t\n"
    "\tGiven the observed autocorrelation structure, a decorrelation time\n"
    "\t'T0' is computed by ar::decorrelation_time and used to estimate\n"
    "\tthe effective signal variance 'eff_var'.  The number of effectively\n"
    "\tindependent samples is returned in 'eff_N'.  These effective values\n"
    "\tare combined to estimate the sampling error (i.e. the standard\n"
    "\tdeviation of the sample mean) as field 'mu_sigma'.  The absolute\n"
    "\tvalue of the autocorrelation function will be used in computing the\n"
    "\tdecorrelation times whenever absrho is true.\n"
    "\t\n"
    "\tOne may simulate N samples from a fitted process analogously to\n"
    "\t\n"
    "\t\tx = mu + filter([1], AR, sqrt(sigma2eps)*randn(N,1));\n"
    "\t\n"
    "\tWhen omitted, submean defaults to " STRINGIFY(DEFAULT_SUBMEAN) ".\n"
    "\tWhen omitted, absrho defaults to " STRINGIFY(DEFAULT_ABSRHO) ".\n"
    "\tWhen omitted, maxorder defaults to " STRINGIFY(DEFAULT_MAXORDER) ".\n"
)
{
    std::size_t maxorder = DEFAULT_MAXORDER;
    bool        absrho   = DEFAULT_ABSRHO;
    bool        submean  = DEFAULT_SUBMEAN;
    Matrix      data;
    switch (args.length())
    {
        case 4: maxorder = args(3).ulong_value();
        case 3: absrho   = args(2).bool_value();
        case 2: submean  = args(1).bool_value();
        case 1: data     = args(0).matrix_value();
                if (!error_state) break;
        default:
            error("Invalid call to arsel.  Correct usage is: ");
        case 0:
            print_usage();
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
    ColumnVector _T0       (M);

    // Prepare vectors to capture burg_method() output
    std::vector<double> params, sigma2e, gain, autocor;
    params .reserve(maxorder*(maxorder + 1)/2);
    sigma2e.reserve(maxorder + 1);
    gain   .reserve(maxorder + 1);
    autocor.reserve(maxorder + 1);

    // Prepare repeatedly-used working storage for burg_method()
    std::vector<double> f, b, Ak, ac;

    // Process each signal in turn...
    for (octave_idx_type i = 0; i < M; ++i)
    {
        // Use burg_method to estimate a hierarchy of AR models from input data
        params .clear();
        sigma2e.clear();
        gain   .clear();
        autocor.clear();
        ar::strided_adaptor<double*> signal_begin(&data(i,0), M);
        ar::strided_adaptor<double*> signal_end  (&data(i,N), M);
        ar::burg_method(signal_begin, signal_end, _mu(i), maxorder,
                        std::back_inserter(params),
                        std::back_inserter(sigma2e),
                        std::back_inserter(gain),
                        std::back_inserter(autocor),
                        submean, /* output hierarchy? */ true, f, b, Ak, ac);

        // Keep only best model according to CIC accounting for subtract_mean.
        // TODO Permit specifying the criterion as a function argument.
        if (submean)
        {
            ar::best_model<ar::CIC<ar::Burg<ar::mean_subtracted> > >(
                    N, params, sigma2e, gain, autocor);
        }
        else
        {
            ar::best_model<ar::CIC<ar::Burg<ar::mean_retained> > >(
                    N, params, sigma2e, gain, autocor);
        }

        // Compute decorrelation time from the estimated autocorrelation model
        ar::predictor<double> p = ar::autocorrelation(
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
    retval.assign("data",      data);
    retval.assign("eff_N",     _eff_N);
    retval.assign("eff_var",   _eff_var);
    retval.assign("gain",      _gain);
    retval.assign("maxorder",  octave_value(maxorder));
    retval.assign("mu",        _mu);
    retval.assign("mu_sigma",  _mu_sigma);
    retval.assign("N",         octave_value(N));
    retval.assign("sigma2eps", _sigma2eps);
    retval.assign("submean",   octave_value(submean));
    retval.assign("T0",        _T0);

    return octave_value_list(retval);
}
