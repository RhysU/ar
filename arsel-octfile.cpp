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
#include <iostream>
#include <limits>
#include <list>
#include <vector>

#include <octave/oct.h>

/** @file
 * A GNU Octave wrapper estimating the best AR(p) model given signal input.
 * Compare \ref arsel.cpp.
 */

// Compile-time defaults in the code also appearing in the help message
#define DEFAULT_SUBMEAN  true
#define DEFAULT_MAXORDER 512
#define STRINGIFY(x) STRINGIFY_HELPER(x)
#define STRINGIFY_HELPER(x) #x

DEFUN_DLD(
    arsel, args, nargout,
    "\t[A, mu, sigma2eps, eff_sigma2x, Neff, T0] = arsel (d, submean, maxorder)\n"
    "\tAutomatically fit an autoregressive model to an input signal.\n"
    "\t\n"
    "\tUse ar::burg_method followed by ar::best_model<CIC<Burg<MeanHandling > >\n"
    "\tto obtain the most likely autoregressive process for input signal d.\n"
    "\tThe sample mean of d will be subtracted whenever submean is true.\n"
    "\tModel orders zero through min(size(d), maxorder) will be considered.\n"
    "\t\n"
    "\tThe filter()-ready process parameters are returned in A, the sample mean\n"
    "\tin mu, and the innovation variance \\sigma^2_\\epsilon in sigma2eps.\n"
    "\tGiven the observed autocorrelation structure, a decorrelation time T0 is\n"
    "\tcomputed and used to estimate the effective signal variance eff_sigma2x.\n"
    "\tThe number of effectively independent samples is returned in Neff.\n"
    "\t\n"
    "\tOne may simulate N samples from the fitted process by calling:\n"
    "\t\n"
    "\t\tx = mu + filter([1], A, sqrt(sigma2eps)*randn(N,1));\n"
    "\t\n"
    "\tWhen not supplied, submean defaults to " STRINGIFY(DEFAULT_SUBMEAN) ".\n"
    "\tWhen not supplied, maxorder defaults to " STRINGIFY(DEFAULT_MAXORDER) ".\n"
    "\tWhen no output arguments are requested, all details are displayed.\n"
)
{
    std::size_t maxorder = DEFAULT_MAXORDER;
    bool        submean  = DEFAULT_SUBMEAN;
    RowVector   data;
    switch (args.length())
    {
        case 3:
            maxorder = args(2).ulong_value();
        case 2:
            submean = args(1).bool_value();
        case 1:
            data = args(0).row_vector_value();
            break;
        default:
            error("Invalid call to arsel.  Correct usage is: ");
        case 0:
            print_usage();
            return octave_value();
    }

    // How much data will we be processing?
    const size_t N = data.dims().numel();
    if (!N) {
        error("arsel: input signal d must have nonzero length");
        return octave_value();
    }

    // Use burg_method to estimate a hierarchy of AR models from input data
    double mu;
    std::vector<double> params, sigma2e, gain, autocor;
    params .reserve(maxorder*(maxorder + 1)/2);
    sigma2e.reserve(maxorder + 1);
    gain   .reserve(maxorder + 1);
    autocor.reserve(maxorder + 1);
    ar::burg_method(data.fortran_vec(), data.fortran_vec() + N,
                    mu, maxorder,
                    std::back_inserter(params),
                    std::back_inserter(sigma2e),
                    std::back_inserter(gain),
                    std::back_inserter(autocor),
                    submean, /* output hierarchy? */ true);

    // Keep only best model according to CIC accounting for subtract_mean.
    if (submean) {
        ar::best_model<ar::CIC<ar::Burg<ar::mean_subtracted> > >(
                N, params, sigma2e, gain, autocor);
    } else {
        ar::best_model<ar::CIC<ar::Burg<ar::mean_retained> > >(
                N, params, sigma2e, gain, autocor);
    }

    // Compute decorrelation time from the estimated autocorrelation model
    double T0 = std::numeric_limits<double>::quiet_NaN();
    if (nargout == 0 || nargout > 2) {
        ar::predictor<double> p = ar::autocorrelation(
                params.begin(), params.end(), gain[0], autocor.begin());
        T0 = ar::decorrelation_time(N, p, /* abs(rho) */ true);
    }

    // Prepare autoregressive process parameters for Octave's filter function
    RowVector A(params.size() + 1, /* leading one */ 1);
    std::copy(params.begin(), params.end(), A.fortran_vec() + 1);

    // Prepare output: [A, mu, sigma2eps, eff_sigma2x, Neff, T0]
    octave_value_list retval;
    std::list<std::string> names;
    switch (nargout == 0 ? std::numeric_limits<int>::max() : nargout)
    {
        default:
            if (nargout) warning("arsel: too many output values requested");
        case 6:
            retval.prepend(T0);
            names.push_front("T0");
        case 5:
            retval.prepend(N / T0);
            names.push_front("Neff");
        case 4:
            retval.prepend((N*gain[0]*sigma2e[0]) / (N - T0));
            names.push_front("eff_sigma2x");
        case 3:
            retval.prepend(sigma2e[0]);
            names.push_front("sigma2eps");
        case 2:
            retval.prepend(mu);
            names.push_front("mu");
        case 1:
            retval.prepend(A);
            names.push_front("A");
            break;
        case 0:
            panic_impossible();
    }
    retval.stash_name_tags(names);

    // Provide no results whenever errors are detected
    if (error_state) retval.resize(0);

    // FIXME Use Octave's default stream (?) instead of std::cout
    // Display results to user when no output was requested
    if (!nargout) {
        int i = 0;
        while (names.size()) {
            retval(i++).print_with_name(std::cout, names.front());
            names.pop_front();
        }
        retval.resize(0);
    }

    return retval;
}
