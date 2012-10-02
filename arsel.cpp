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
#include "optionparser.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

// Command line argument declarations for optionparser.h usage
enum OptionIndex { UNKNOWN, HELP, SUBMEAN, SIGNRHO };
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "",      option::Arg::None,
     "Usage: arsel [OPTION]...\n"
     "Fit an optimal autoregressive model to data from standard input.\n"
     "\n"
     "Options:\n" },
    {SUBMEAN, 0,  "",  "subtract-mean", option::Arg::None,
     "\t--subtract-mean  \tSubtract the sample mean from the data" },
    {SIGNRHO, 0,  "",  "signed-rho",    option::Arg::None,
     "\t--signed-rho     \tUse autocorrelation, including sign, when computing T0" },
    {HELP,    0,  "h", "help",        option::Arg::None,
     "\t-h,--help        \tDisplay this help and exit" },
    {0,0,0,0,0,0}
};

int main(int argc, char *argv[])
{
    using namespace ar;
    using namespace std;

    // TODO Add processing of files specified on the command line

    // Process any incoming command line arguments using optionparser.h
    argc-=(argc>0); argv+=(argc>0); // Skip program name argv[0] if present
    option::Stats stats(usage, argc, argv);
    struct Guard { // scoped_array or unique_ptr would simply RAII here
        option::Option *options, *buffer;
        Guard(option::Option *options, option::Option *buffer)
            : options(options), buffer(buffer) {}
        ~Guard() { delete[] options; delete[] buffer; }
    } guard(new option::Option[stats.options_max],
            new option::Option[stats.buffer_max ]);
    option::Parser parse(usage, argc, argv, guard.options, guard.buffer);
    if (parse.error() || guard.options[UNKNOWN]) {
        for (option::Option* o = guard.options[UNKNOWN]; o; o = o->next()) {
            cerr << "Unknown option: " << o->name << "\n";
        }
        return EXIT_FAILURE;
    }
    if (guard.options[HELP]) {
        option::printUsage(std::cout, usage);
        return EXIT_SUCCESS;
    }
    const bool subtract_mean = guard.options[SUBMEAN] ? true : false;
    const bool absolute_rho  = guard.options[SIGNRHO] ? false : true;

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
                params.end(), gain[0], autocor.begin()), absolute_rho);

    // Output details about the best model and derived information
    // Unbiased effective variance expression from [Trenberth1984]
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos
         <<   "# N                   "   << N
         << "\n# AR(p)               "   << params.size()
         << "\n# Mean                "   << mean
         << "\n# \\sigma^2_\\epsilon   " << sigma2e[0]
         << "\n# Variance_effective  "   << (N*gain[0]*sigma2e[0]) / (N - T0)
         << "\n# N_effective         "   << N / T0
         << "\n# t_decorrelation     "   << T0
         << '\n';
    copy(params.begin(), params.end(), ostream_iterator<double>(cout,"\n"));
    cout.flush();

    return EXIT_SUCCESS;
}
