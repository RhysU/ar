// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "burg.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>

// Computes percent difference of \c b against theoretical result \c a.
template<typename FPT> FPT pdiff(FPT a, FPT b) { return (b - a) / a * 100; }

int main(int argc, char *argv[])
{
    using namespace std;

    // Lag 1 through lag 4 autocorrelations specified
    vector<double> rho;
    rho.push_back(-1./2);
    rho.push_back( 1./2);
    rho.push_back( 1./10);
    rho.push_back( 1./20);

    // AR(4) process parameters computed in Octave using [a,v,k]=arburg(rho,4)
    vector<double> params;
    params.push_back( 125./154);
    params.push_back(- 23./ 70);
    params.push_back(- 61./ 70);
    params.push_back(- 31./ 77);

    // Compute autocorrelation functions by way of reflection coefficients
    //
    // This is even less numerically stable than the usual Yule-Walker case
    // because we use reverse Levinson recursion followed by forward recursion
    // to avoid needing auxiliary storage.
    vector<double> k(params);
    reflection_coefficients(k.begin(), k.end());
    vector<double> rho_est(rho.size());
    autocorrelations(k.begin(), k.end(), rho_est.begin());  // k overwritten

    printf("%22s %22s %22s\n", "Autocorrelation", "algorithm", "PercentDiff");
    printf("%22s %22s %22s\n", "---------------", "---------", "-----------");
    for (vector<double>::const_iterator i = rho.begin(), j = rho_est.begin();
         i != rho.end(); ++i, ++j)
    {
        printf("%22.14g %22.14g %22.14g\n", *i, *j, pdiff(*i, *j));
    }

    return 0;
}
