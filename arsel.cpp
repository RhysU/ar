// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "burg.hpp"

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

// Estimate the best AR(p) model given data on standard input
int main(int argc, char *argv[])
{
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
    const std::size_t N = burg_method(std::istream_iterator<double>(cin),
                                      std::istream_iterator<double>(),
                                      mean,
                                      order,
                                      std::back_inserter(params),
                                      std::back_inserter(sigma2e),
                                      std::back_inserter(gain),
                                      std::back_inserter(autocor),
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

    // Output some details about the best model
    cout.precision(std::numeric_limits<double>::digits10 + 2);
    cout << showpos
         <<   "# N                   "   << N
         << "\n# AR(p)               "   << params.size()
         << "\n# Mean                "   << mean
         << "\n# \\sigma^2_\\epsilon   " << sigma2e[0]
         << "\n# Gain                "   << gain[0]
         << "\n# \\sigma^2_x          "  << gain[0]*sigma2e[0]
         << '\n';
    copy(params.begin(), params.end(), ostream_iterator<double>(cout,"\n"));
    cout.flush();

    return 0;
}
