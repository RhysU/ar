// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "burg.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>

// Computes percent difference of \c b against theoretical result \c a.
template<typename FPT> FPT pdiff(FPT a, FPT b) { return (b - a) / a * 100; }

template<typename T> struct sum_error : public std::binary_function<T,T,T> {
    T operator() (T a, T b) {using std::abs; return a + abs(b);}
};

// Set working precision
typedef double real;

// Test burg_algorithm against synthetic data
int main(int argc, char *argv[])
{
    using namespace std;

    vector<real> exact;
    {
        ifstream f;
        f.exceptions(ifstream::badbit);
        f.open(argc > 1 ? argv[1] : "test1.coeff");
        copy(istream_iterator<real>(f), istream_iterator<real>(),
             back_inserter(exact));
    }
    vector<real> est(exact.size()), cor(exact.size());

    vector<real> data;
    {
        ifstream f;
        f.exceptions(ifstream::badbit);
        f.open(argc > 2 ? argv[2] : "test1.dat");
        copy(istream_iterator<real>(f), istream_iterator<real>(),
             back_inserter(data));
    }

    // Use burg_algorithm to fit an AR model and characterize it completely
    real sigma2e, gain;
    burg_algorithm(data.begin(), data.end(), est.begin(), est.end(),
                   &sigma2e, &gain, cor.begin());

    // Solve Yule-Walker equations using Zohar's algorithm as consistency check
    // Given right hand side containing rho_1, ..., rho_p the solution should
    // be -a_1, ..., -a_p on success so adding to it a_1, ..., a_p gives errors.
    vector<real> rhs(cor);
    zohar_linear_solve(cor.begin(), --cor.end(), cor.begin(), rhs.begin());
    transform(rhs.begin(), rhs.end(), est.begin(), rhs.begin(), plus<real>());
    real res = accumulate(rhs.begin(), rhs.end(), real(0), sum_error<real>());

    printf("%22s %22s %22s\n", "Coefficient", "burg_algorithm", "PercentDiff");
    printf("%22s %22s %22s\n", "-----------", "--------------", "-----------");
    for (vector<real>::const_iterator i = exact.begin(), j = est.begin();
         i != exact.end(); ++i, ++j)
    {
        printf("%22.14g %22.14g %22.14g\n", *i, *j, pdiff(*i, *j));
    }
    printf("\n");
    printf("%22s %22.14g\n", "Yule-Walker residual:", res);
    printf("%22s %22.14g\n", "\\sigma^2_\\epsilon:",  sigma2e);
    printf("%22s %22.14g\n", "signal gain:",          gain);

    return 0;
}
