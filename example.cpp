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

// Example program using burg_method modified from Collomb's sample
int main()
{
    using namespace std;

    // Create data to approximate.  Sample size decreased to show selection.
    vector<long double> data(10, 0.0);
    for (size_t i = 0; i < data.size(); i++)
    {
        data[i] =     cos(i*0.01) + 0.75*cos(i*0.03)
                + 0.5*cos(i*0.05) + 0.25*cos(i*0.11);
    }

    // Get linear prediction coefficients for orders 1 through order
    size_t maxorder = 9;
    long double mean;
    vector<long double> params, sigma2e, gain, autocor;
    burg_method(data.begin(), data.end(), mean, maxorder,
                back_inserter(params), back_inserter(sigma2e),
                back_inserter(gain), back_inserter(autocor),
                false, true);

    // Display orders, mean squared discrepancy, and model coefficients
    printf("%2s  %8s %8s %8s\n", "AR", "RMS/N", "Gain", "Coefficients");
    printf("%2s  %8s %8s %8s\n", "--", "-----", "----", "------------");
    for (size_t p = 0, c = 0; p < maxorder; ++p) {
        printf("%2lu  %8.3Lg %8.3Lg", p+1, sigma2e[p], gain[p]);
        for (size_t i = 1; i < p+1; ++i) printf(" %8.3Lg", params[c++]);
        putchar('\n');
    }

    // Display model selection results
    printf("\n");
    vector<long double>::difference_type best;
    best = 1 + select_model<AIC>(data.size(), sigma2e.begin(), sigma2e.end());
    printf("AIC  selects model order %d as best\n", (int) best);
    best = 1 + select_model<AICC>(data.size(), sigma2e.begin(), sigma2e.end());
    printf("AICC selects model order %d as best\n", (int) best);
    best = 1 + select_model<CIC<Burg<mean_retained> > >(
            data.size(), sigma2e.begin(), sigma2e.end());
    printf("CIC  selects model order %d as best\n", (int) best);

    return 0;
}
