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

template <class Real, class In>
std::size_t process(In& in,
                    Real& mean,
                    std::size_t &order,
                    std::vector<Real>& params,
                    std::vector<Real>& sigma2e,
                    std::vector<Real>& gain,
                    std::vector<Real>& autocor,
                    bool subtract_mean)
{
    // Use burg_method to estimate a hierarchy of AR models from input data
    params .reserve(order*(order + 1)/2);
    sigma2e.reserve(order + 1);
    gain   .reserve(order + 1);
    autocor.reserve(order + 1);
    const std::size_t N = burg_method(std::istream_iterator<Real>(in),
                                      std::istream_iterator<Real>(),
                                      mean,
                                      order,
                                      std::back_inserter(params),
                                      std::back_inserter(sigma2e),
                                      std::back_inserter(gain),
                                      std::back_inserter(autocor),
                                      subtract_mean,
                                      /* output hierarchy? */ true);

    // Find the best model according to CIC accounting for subtract_mean.
    typename std::vector<Real>::difference_type best;
    if (subtract_mean) {
        best = select_model<CIC<Burg<mean_subtracted> > >(
                    N, 0u, sigma2e.begin(), sigma2e.end());
    } else {
        best = select_model<CIC<Burg<mean_retained> > >(
                    N, 0u, sigma2e.begin(), sigma2e.end());
    }

    // Trim away everything but the best model (ranges might overlap)
    std::copy_backward(params.begin() + best*(best+1)/2,
                       params.begin() + best*(best+1)/2 + best,
                       params.begin() + best);
    params.resize(best);
    sigma2e[0] = sigma2e[best]; sigma2e.resize(1);
    gain   [0] = gain   [best]; gain   .resize(1);
    autocor.resize(best + 1);

    return N;
}

// Test burg_method against synthetic data
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

    double mean;
    size_t order = 512;
    vector<double> params, sigma2e, gain, autocor;
    size_t N = process(cin, mean, order, params,
                       sigma2e, gain, autocor, subtract_mean);

    cout.precision(std::numeric_limits<double>::digits10 + 2);
    cout << showpos
         <<   "# N                   "   << N
         << "\n# Mean                "   << mean
         << "\n# \\sigma^2_\\epsilon   " << sigma2e[0]
         << "\n# Gain                "   << gain[0]
         << "\n# \\sigma^2_x          "  << gain[0]*sigma2e[0]
         << '\n';
    copy(params.begin(), params.end(), ostream_iterator<double>(cout,"\n"));
    cout.flush();

    return 0;
}
