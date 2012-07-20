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

// Test burg_algorithm against synthetic data and NR routine memcof()
int main(int argc, char *argv[])
{
    using namespace std;

    vector<double> exact;
    {
        ifstream f;
        f.exceptions(ifstream::badbit);
        f.open(argc > 1 ? argv[1] : "test1.coeff");
        copy(istream_iterator<double>(f), istream_iterator<double>(),
             back_inserter(exact));
    }
    vector<double> est(exact.size()), nr(exact.size());

    vector<double> data;
    {
        ifstream f;
        f.exceptions(ifstream::badbit);
        f.open(argc > 2 ? argv[2] : "test1.dat");
        copy(istream_iterator<double>(f), istream_iterator<double>(),
             back_inserter(data));
    }

    // Compute burg_algorithm's answer
    double est_msd;
    burg_algorithm(data.begin(), data.end(), est.begin(), est.end(), &est_msd);

    printf("%22s %22s %22s\n", "Coefficient", "burg_algorithm", "PercentDiff");
    printf("%22s %22s %22s\n", "-----------", "--------------", "-----------");
    for (vector<double>::const_iterator i = exact.begin(), j = est.begin();
         i != exact.end(); ++i, ++j)
    {
        printf("%22.14g %22.14g %22.14g\n", *i, *j, pdiff(*i, *j));
    }
    printf("\n%22s %22.14g\n", "Mean^2 Discrepancy", est_msd);

// FIXME
//  // Compute autocorrelation functions by way of reflection coefficients 
//  vector<double> rho(est);
//  reflection_coefficients(rho.begin(), rho.end());
//  autocorrelations(est.begin(), est.end(), rho.begin(), rho.begin());
//  for (vector<double>::const_iterator i = rho.begin(); i != rho.end(); ++i) {
//      printf("%22.14g\n", *i);
//  }

    return 0;
}
