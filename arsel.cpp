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

// Forward declarations for argument checking logic
struct Arg : public option::Arg
{
    static option::ArgStatus NonEmpty(const option::Option& opt, bool msg);
    static option::ArgStatus Numeric (const option::Option& opt, bool msg);
};

// Command line argument declarations for optionparser.h usage
enum OptionIndex { UNKNOWN, CRITERION, HELP, SUBMEAN, SIGNRHO };
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "",      option::Arg::None,
     "Usage: arsel [OPTION]...\n"
     "Fit an optimal autoregressive model to data from standard input.\n"
     "\n"
     "Options:" },
    {0,0,"","",Arg::None,0}, // table break
    {CRITERION, 0,  "c",  "criterion",     Arg::NonEmpty,
     "  -c \t--criterion=ABBREV  \tUse the specified model order selection criterion" },
    {HELP,      0,  "h", "help",           Arg::None,
     "  -h \t--help        \tDisplay this help message and immediately exit" },
    {SUBMEAN,   0,  "m",  "subtract-mean", Arg::None,
     "  -m \t--subtract-mean  \tSubtract the sample mean from the incoming data" },
    {SIGNRHO,   0,  "r",  "signed-rho",    Arg::None,
     "  -r \t--signed-rho     \tUse signed autocorrelation values when computing T0" },
    {0,0,0,0,0,0}
};

int main(int argc, char *argv[])
{
    using namespace std;

    // TODO Add processing of files specified on the command line

    // Parse and process any command line arguments using optionparser.h
    string criterion   = "";
    bool subtract_mean = false;
    bool absolute_rho  = true;
    {
        using namespace option;

        Stats stats(usage, argc-(argc>0), argv+(argc>0));

        // scoped_array or unique_ptr would simply RAII but add dependencies
        struct Guard {
            Option *options, *buffer;
            Guard(Option *options, Option *buffer)
                : options(options), buffer(buffer) {}
            ~Guard() { delete[] options; delete[] buffer; }
        } g(new Option[stats.options_max], new Option[stats.buffer_max]);

        Parser parse(usage, argc-(argc>0), argv+(argc>0), g.options, g.buffer);

        if (parse.error() || g.options[UNKNOWN]) {
            for (Option* o = g.options[UNKNOWN]; o; o = o->next()) {
                cerr << "Unknown option: " << o->name << "\n";
            }
            return EXIT_FAILURE;
        }

        if (g.options[HELP]) {
            printUsage(cout, usage);
            return EXIT_SUCCESS;
        }

        if (g.options[CRITERION])
            criterion = g.options[CRITERION].last()->arg;

        if (g.options[SUBMEAN])
            subtract_mean = true;

        if (g.options[SIGNRHO])
            absolute_rho = false;
    }

    // Look up desired model selection criterion using ar::best_model_function
    // best_model_function template parameters fit ar::burg_method usage below
    typedef ar::best_model_function<
            ar::Burg,size_t,vector<double> > best_model_function;
    const best_model_function::type best_model
            = best_model_function::lookup(criterion, subtract_mean);
    if (!best_model) {
        cerr << "Unknown model selection criterion: " << criterion << "\n";
        return EXIT_FAILURE;
    }

    // Use burg_method to estimate a hierarchy of AR models from input data
    double mean;
    size_t order = 512;
    vector<double> params, sigma2e, gain, autocor;
    params .reserve(order*(order + 1)/2);
    sigma2e.reserve(order + 1);
    gain   .reserve(order + 1);
    autocor.reserve(order + 1);
    const size_t N = ar::burg_method(istream_iterator<double>(cin),
                                     istream_iterator<double>(),
                                     mean,
                                     order,
                                     back_inserter(params),
                                     back_inserter(sigma2e),
                                     back_inserter(gain),
                                     back_inserter(autocor),
                                     subtract_mean,
                                     /* output hierarchy? */ true);


    // Keep only best model according to selected criterion
    best_model(N, params, sigma2e, gain, autocor);

    // Compute decorrelation time from the estimated autocorrelation model
    const double T0  = ar::decorrelation_time(N, ar::autocorrelation(
                params.begin(), params.end(), gain[0], autocor.begin()),
                absolute_rho);

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

// Argument checking logic suggested by optionparser.h example routines
option::ArgStatus Arg::NonEmpty(const option::Option& opt, bool msg)
{
    if (opt.arg != 0 && opt.arg[0] != 0) {
        return option::ARG_OK;
    }

    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                    << " requires a non-empty argument\n";
    }
    return option::ARG_ILLEGAL;
}

// Argument checking logic suggested by optionparser.h example routines
option::ArgStatus Arg::Numeric(const option::Option& opt, bool msg)
{
    char *endptr;
    if (opt.arg != 0 && strtol(opt.arg, &endptr, 10) && *endptr == 0) {
        return option::ARG_OK;
    }

    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                    << " requires a numeric argument\n";
    }
    return option::ARG_ILLEGAL;
}
