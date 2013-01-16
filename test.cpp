// Copyright (C) 2012, 2013 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Test \ref ar::burg_method against synthetic data.
 */

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

#include "ar.hpp"
#include "optionparser.h"
#include "real.hpp"

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

// Command line argument declarations for optionparser.h usage
enum OptionIndex {
    UNKNOWN, HELP, SUBMEAN
};
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "",      option::Arg::None,
     "Usage: test COEFFICIENTS DATA\n"
     "Compute ar::burg_method deviation from known coefficients on test data ("
         /* Working precision */ STRINGIFY(REAL) ").\n"
     "\n"
     "Options:" },
    {0,0,"","",option::Arg::None,0}, // table break
    {HELP,    0,  "h", "help",          option::Arg::None,
     "  -h \t--help   \tDisplay this help message and immediately exit" },
    {SUBMEAN, 0,  "s", "subtract-mean", option::Arg::None,
     "  -s \t--subtract-mean  \tSubtract the sample mean from the incoming data" },
    {0,0,0,0,0,0}
};

// Computes percent difference of \c b against theoretical result \c a.
template<typename FPT> FPT pdiff(FPT a, FPT b) { return (b - a) / a * 100; }

template<typename T> struct sum_error : public std::binary_function<T,T,T> {
    T operator() (T a, T b) {using std::abs; return a + abs(b);}
};

// Test burg_method against synthetic data
int main(int argc, char *argv[])
{
    using namespace ar;
    using namespace std;

    string filename_coeffs, filename_data;
    bool subtract_mean = false;
    {
        option::Stats stats(usage, argc-(argc>0), argv+(argc>0));

        vector<option::Option> options(stats.options_max + stats.buffer_max);

        option::Parser parse(usage, argc-(argc>0), argv+(argc>0),
                             &options[0], &options[stats.options_max]);

        if (parse.error() || options[UNKNOWN]) {
            for (option::Option* o = options[UNKNOWN]; o; o = o->next()) {
                cerr << "Unknown option: " << o->name << "\n";
            }
            return EXIT_FAILURE;
        }

        if (options[HELP]) {
            printUsage(cout, usage);
            return EXIT_SUCCESS;
        }

        if (options[SUBMEAN])
            subtract_mean = true;

        if (parse.nonOptionsCount() != 2) {
            cerr << "Missing one or more file operands.  Try --help.\n";
            return EXIT_FAILURE;
        }
        filename_coeffs = parse.nonOption(0);
        filename_data = parse.nonOption(1);
    }

    // Read "exact" coefficients
    vector<real> exact;
    {
        ifstream f;
        f.exceptions(ifstream::badbit);
        f.open(filename_coeffs.c_str());
        copy(istream_iterator<real>(f), istream_iterator<real>(),
             back_inserter(exact));
    }
    vector<real> est(exact.size()), cor(exact.size() + 1);

    // Read time series data
    ifstream f;
    f.exceptions(ifstream::badbit);
    f.open(filename_data.c_str());

    // Use burg_method to fit an AR model and characterize it completely
    size_t maxorder = exact.size();
    real mean, sigma2e, gain;
    burg_method(istream_iterator<real>(f), istream_iterator<real>(),
                mean, maxorder, est.begin(), &sigma2e, &gain, cor.begin(),
                subtract_mean, false);

    // Close input file
    f.close();

    // Solve Yule-Walker equations using Zohar's algorithm as consistency check
    // Given right hand side containing rho_1, ..., rho_p the solution should
    // be -a_1, ..., -a_p on success so adding to it a_1, ..., a_p gives errors.
    vector<real> rhs(++cor.begin(), cor.end());
    zohar_linear_solve(++cor.begin(), --cor.end(), rhs.begin());
    transform(rhs.begin(), rhs.end(), est.begin(), rhs.begin(), plus<real>());
    real res = accumulate(rhs.begin(), rhs.end(), real(0), sum_error<real>());

    printf("%22s %22s %22s\n", "Nominal Value", "burg_method", "Percent Diff");
    printf("%22s %22s %22s\n", "-------------", "-----------", "------------");
    for (vector<real>::const_iterator i = exact.begin(), j = est.begin();
         i != exact.end(); ++i, ++j)
    {
        printf("%22.14g %22.14g %22.14g\n",
               (double) *i, (double) *j, (double) pdiff(*i, *j));
    }
    printf("\n");
    printf("%22s %22.14g\n", "mean of data:",         (double) mean          );
    printf("%22s %22.14g\n", "\\sigma^2_\\epsilon:",  (double) sigma2e       );
    printf("%22s %22.14g\n", "signal gain:",          (double) gain          );
    printf("%22s %22.14g\n", "\\sigma^2_x:",          (double) (gain*sigma2e));
    printf("%22s %22.14g\n", "Yule-Walker residual:", (double) res           );

    return 0;
}
