// Copyright (C) 2012 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Estimate the best AR(p) model given data on standard input.  Illustrates
 * \ref ar::autocorrelation, \ref ar::predictor,
 * and \ref ar::decorrelation_time.
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
    // Along the way, save the variance as computed for model order zero.
    double var;
    if (subtract_mean) {
        var = sigma2e[0];              // Already centered
        best_model<CIC<Burg<mean_subtracted> > >(
                N, params, sigma2e, gain, autocor);
    } else {
        var = sigma2e[0] - mean*mean;  // Uncentered so remove mean^2
        best_model<CIC<Burg<mean_retained> > >(
                N, params, sigma2e, gain, autocor);
    }

    // Compute decorrelation time from the estimated autocorrelation model
    const double T0  = decorrelation_time(N, autocorrelation(params.begin(),
                params.end(), gain[0], autocor.begin()), /* abs(rho) */ true);

    // Output details about the best model and derived information
    // Unbiased effective variance expression from [Trenberth1984]
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos
         <<   "# N                   "   << N
         << "\n# AR(p)               "   << params.size()
         << "\n# Mean                "   << mean
         << "\n# \\sigma^2_\\epsilon   " << sigma2e[0]
         << "\n# Gain                "   << gain[0]
         << "\n# Variance            "   << var
         << "\n# t_decorrelation     "   << T0
         << "\n# N_effective         "   << N / T0
         << "\n# Variance_effective  "   << (N*var) / (N - T0)
         << '\n';
    copy(params.begin(), params.end(), ostream_iterator<double>(cout,"\n"));
    cout.flush();

    return 0;
}
