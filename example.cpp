// Except for any way in which it interferes with Cedrick Collomb's 2009
// copyright assertion in the article "Burg’s Method, Algorithm and Recursion":
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
#include "burg.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>

// Example program using burg_method and model selection
int main(int argc, char *argv[])
{
    using namespace std;

    // Maximum model estimation order to perform (limited by N)
    size_t maxorder = argc > 1 && atol(argv[1]) >= 0
                      ? atol(argv[1])
                      : 7;

    // Default sample size shows model selection criteria differences
    const std::size_t N = argc > 2 && atol(argv[2]) > 0
                        ? atol(argv[2])
                        : 10;


    // Create data to approximate from [Collomb2009]'s example.
    //
    // Notice that this isn't the most nicely conditioned problem
    // and truncating coefficients (as we'll do for display below)
    // causes roots to appear outside the unit circle.
    vector<long double> data(N, 0.0);
    for (size_t i = 0; i < N; i++)
    {
        data[i] =     cos(i*0.01) + 0.75*cos(i*0.03)
                + 0.5*cos(i*0.05) + 0.25*cos(i*0.11);
    }

    // Estimate process parameters using Burg's method
    printf("Estimating at most an AR(%lu) model using %lu samples\n\n",
           maxorder, N);

    long double mean;
    vector<long double> params, sigma2e, gain, autocor;
    burg_method(data.begin(), data.end(), mean, maxorder,
                back_inserter(params), back_inserter(sigma2e),
                back_inserter(gain), back_inserter(autocor),
                /* subtract mean?    */ false,
                /* output hierarchy? */ true);

    // Display orders, mean squared discrepancy, and model coefficients
    printf("%2s  %9s %9s %s\n", "AR", "RMS/N", "Gain", "Filter Coefficients");
    printf("%2s  %9s %9s %s\n", "--", "-----", "----", "-------------------");
    for (size_t p = 0, c = 0; p <= maxorder; ++p)
    {
        printf("%2lu  %9.2Le %9.2Le [ 1 ", p, sigma2e[p], gain[p]);
        for (size_t i = 0; i < p; ++i)
            printf(" %8.4Lg", params[c++]);
        printf(" ]\n");
    }

    // Display model selection results
    printf("\n");
    vector<long double>::difference_type best;
    best = evaluate_models<AIC>(N, 0u, sigma2e.begin(), sigma2e.end());
    printf("AIC  selects model order %d as best\n", (int) best);
    best = evaluate_models<AICC>(N, 0u, sigma2e.begin(), sigma2e.end());
    printf("AICC selects model order %d as best\n", (int) best);
    best = evaluate_models<CIC<Burg<mean_retained> > >(
                N, 0u, sigma2e.begin(), sigma2e.end());
    printf("CIC  selects model order %d as best\n", (int) best);

    return 0;
}
