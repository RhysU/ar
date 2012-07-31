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
#include <vector>

#include <octave/oct.h>

/** @file
 * A GNU Octave wrapper estimating the best AR(p) model given signal input.
 * Compare \ref arsel.cpp.
 */

DEFUN_DLD(
    arsel, args, nargout,
    "[A, mu, sigma2eps, eff_sigma2x, T0] = arsel (d, submean=true, order=512)\n"
    "\n"
    "Use ar::burg_method followed by ar::best_model<CIC<Burg<MeanHandling > >\n"
    "to obtain the most likely autoregressive process for input signal d.\n"
    "The sample mean of d will be subtracted whenever submean is true.\n"
    "Model orders zero through min(size(d), order) will be considered.\n"
    "\n"
    "The filter()-ready process parameters are returned in A, the sample mean\n"
    "in mu, and the innovation variance \\sigma^2_\\epsilon in sigma2eps.\n"
    "Given the observed autocorrelation structure, a decorrelation time T0 is\n"
    "computed and used to estimate the effective signal variance eff_sigma2x.\n"
)
{
    std::size_t order = 512;
    bool        submean = true;
    RowVector   data;
    switch (args.length())  // Cases deliberately fall through
    {
        case 3:
            order = args(2).ulong_value();
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
    params .reserve(order*(order + 1)/2);
    sigma2e.reserve(order + 1);
    gain   .reserve(order + 1);
    autocor.reserve(order + 1);
    ar::burg_method(data.fortran_vec(), data.fortran_vec() + N,
                    mu, order,
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
    const double T0  = ar::decorrelation_time(N, ar::autocorrelation(
                params.begin(), params.end(), gain[0], autocor.begin()),
                /* abs(rho) */ true);

    // Prepare autoregressive process parameters for Octave's filter function
    RowVector A(params.size() + 1, /* leading one */ 1);
    std::copy(params.begin(), params.end(), A.fortran_vec() + 1);

    // Prepare output: [A, mu, sigma2eps, eff_sigma2x, T0]
    octave_value_list retval;
    switch (nargout)  // Cases deliberately fall through
    {
        default:
            warning("arsel: too many output values requested");
        case 5:
            retval.prepend(T0);
        case 4:
            retval.prepend((N*gain[0]*sigma2e[0]) / (N - T0));
        case 3:
            retval.prepend(sigma2e[0]);
        case 2:
            retval.prepend(mu);
        case 1:
            retval.prepend(A);
            break;
        case 0:
            warning("arsel: no output values requested");
    }

    // Return results iff no errors were detected
    return error_state ? octave_value_list() : retval;
}
