// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Estimate the best AR(p) model given data on standard input.  Illustrates
 * \ref autocorrelation, \ref predictor, and \ref decorrelation_time.
 */

#include "ar.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

int main(int argc, char *argv[])
{
    using namespace ar;
    using namespace std;

    // Process a possible --subtract-mean flag, shifting arguments if needed
    bool subtract_mean = false;
    if (argc > 1 && 0 == strcmp("--subtract-mean", argv[1])) {
        subtract_mean = true;
        argv[1] = argv[0];
        ++argv;
        --argc;
    }

    // Use burg_method to estimate a hierarchy of AR models from input data
    double mean;
    size_t order = 512;
    vector<double> params, sigma2e, gain, autocor;
    params .reserve(order*(order + 1)/2);
    sigma2e.reserve(order + 1);
    gain   .reserve(order + 1);
    autocor.reserve(order + 1);
    const size_t N = burg_method(istream_iterator<double>(cin),
                                 istream_iterator<double>(),
                                 mean,
                                 order,
                                 back_inserter(params),
                                 back_inserter(sigma2e),
                                 back_inserter(gain),
                                 back_inserter(autocor),
                                 subtract_mean,
                                 /* output hierarchy? */ true);

    // Keep only best model according to CIC accounting for subtract_mean.
    if (subtract_mean) {
        best_model<CIC<Burg<mean_subtracted> > >(
                N, params, sigma2e, gain, autocor);
    } else {
        best_model<CIC<Burg<mean_retained> > >(
                N, params, sigma2e, gain, autocor);
    }

    // Compute decorrelation time from the estimated autocorrelation model
    const double T0  = decorrelation_time(N, autocorrelation(params.begin(),
                params.end(), gain[0], autocor.begin()), /* abs(rho) */ true);

    // Output details about the best model and derived information
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos
         <<   "# N                   "   << N
         << "\n# AR(p)               "   << params.size()
         << "\n# Mean                "   << mean
         << "\n# \\sigma^2_\\epsilon   " << sigma2e[0]
         << "\n# Gain                "   << gain[0]
         << "\n# \\sigma^2_x          "  << gain[0]*sigma2e[0]
         << "\n# t_decorrelation     "   << T0
         << "\n# N_effective         "   << N / T0
         << "\n# Variance_effective  "   << (N*gain[0]*sigma2e[0]) / (N - T0)
         << '\n';
    copy(params.begin(), params.end(), ostream_iterator<double>(cout,"\n"));
    cout.flush();

    return 0;
}
