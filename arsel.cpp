// Copyright (C) 2012, 2013 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Estimate the best AR(p) model given data on standard input.  Illustrates
 * \ref ar::autocorrelation, \ref ar::predictor,
 * and \ref ar::decorrelation_time.
 */

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

#include "ar.hpp"
#include "optionparser.h"
#include "real.hpp"

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

// Forward declarations for argument checking logic
struct Arg : public option::Arg
{
    static option::ArgStatus NonEmpty   (const option::Option& opt, bool msg);
    static option::ArgStatus NonNegative(const option::Option& opt, bool msg);
};

// Command line argument declarations for optionparser.h usage
enum OptionIndex {
    UNKNOWN, CRITERION, HELP, MAXORDER, MINORDER, NONABSRHO, SUBMEAN, WINT0
};
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "",      option::Arg::None,
     "Usage: arsel [OPTION]...\n"
     "Fit an optimal autoregressive model to data from standard input ("
         /* Working precision */ STRINGIFY(REAL) ").\n"
     "\n"
     "Options:" },
    {0,0,"","",Arg::None,0}, // table break
    {CRITERION, 0,  "c",  "criterion",         Arg::NonEmpty,
     "  -c \t--criterion=ABBREV  \tUse the specified model selection criterion" },
    {HELP,      0,  "h", "help",               Arg::None,
     "  -h \t--help   \tDisplay this help message and immediately exit" },
    {MINORDER,  0,  "M",  "minorder",          Arg::NonNegative,
     "  -M \t--minorder=MIN  \tConsider only models of at least order AR(p=MIN)" },
    {MAXORDER,  0,  "m",  "maxorder",          Arg::NonNegative,
     "  -m \t--maxorder=MAX  \tConsider only models of at most order AR(p=MAX)" },
    {NONABSRHO, 0,  "n",  "non-absolute-rho",  Arg::None,
     "  -n \t--non-absolute-rho   \tUse non-absolute autocorrelation when computing T0" },
    {SUBMEAN,   0,  "s",  "subtract-mean",     Arg::None,
     "  -s \t--subtract-mean  \tSubtract the sample mean from the incoming data" },
    {WINT0,   0,    "w",  "window-T0",         Arg::NonNegative,
     "  -w \t--window-T0=W  \tIntegrate T0 until W times the data length (default 1)" },
    {0,0,0,0,0,0}
};

int main(int argc, char *argv[])
{
    using namespace std;

    // Parse and process any command line arguments using optionparser.h
    string criterion     = "CIC";
    bool   subtract_mean = false;
    size_t minorder      = 0;
    size_t maxorder      = 512;
    bool   absolute_rho  = true;
    double window_T0     = 1.0;
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

        if (options[CRITERION])
            criterion = options[CRITERION].last()->arg;

        if (options[MAXORDER])
            maxorder = (size_t) strtol(options[MAXORDER].last()->arg, NULL, 10);

        if (options[MINORDER])
            minorder = (size_t) strtol(options[MINORDER].last()->arg, NULL, 10);

        if (options[NONABSRHO])
            absolute_rho = false;

        if (options[SUBMEAN])
            subtract_mean = true;

        if (options[WINT0])
            window_T0 = strtod(options[WINT0].last()->arg, NULL);
    }

    // Look up desired model selection criterion using ar::best_model_function
    // best_model_function template parameters fit ar::burg_method usage below
    typedef ar::best_model_function<
                ar::Burg,size_t,size_t,vector<real>
            > best_model_function;
    const best_model_function::type best_model
            = best_model_function::lookup(criterion, subtract_mean);
    if (!best_model) {
        cerr << "Unknown model selection criterion: " << criterion << "\n";
        return EXIT_FAILURE;
    }

    // Use burg_method to estimate a hierarchy of AR models from input data
    real mu;
    vector<real> params, sigma2e, gain, autocor;
    params .reserve(maxorder*(maxorder + 1)/2);
    sigma2e.reserve(maxorder + 1);
    gain   .reserve(maxorder + 1);
    autocor.reserve(maxorder + 1);
    const size_t N = ar::burg_method(istream_iterator<real>(cin),
                                     istream_iterator<real>(),
                                     mu,
                                     maxorder,
                                     back_inserter(params),
                                     back_inserter(sigma2e),
                                     back_inserter(gain),
                                     back_inserter(autocor),
                                     subtract_mean,
                                     /* output hierarchy? */ true);


    // Keep only best model according to selected criterion
    best_model(N, minorder, params, sigma2e, gain, autocor);

    // Compute decorrelation time from the estimated autocorrelation model
    const real T0 = ar::decorrelation_time(
                static_cast<size_t>(window_T0*N),
                ar::autocorrelation(params.begin(), params.end(),
                                    gain[0], autocor.begin()),
                absolute_rho);

    // Compute effective variance and effective number of independent samples
    const real eff_var  = (N*gain[0]*sigma2e[0]) / (N - T0); // Trenberth1984
    const real eff_N    = N / T0;
    const real mu_sigma = sqrt(eff_var / eff_N);

    // Output details about the best model and derived information
    // Naming conventions here match the output of arsel-octfile by design
    cout.precision(numeric_limits<real>::digits10 + 2);
    cout << boolalpha
         <<   "# absrho    " << absolute_rho
         << "\n# criterion " << criterion
         << "\n# eff_N     " << eff_N
         << "\n# eff_var   " << eff_var
         << "\n# gain      " << gain[0]
         << "\n# maxorder  " << maxorder
         << "\n# minorder  " << minorder
         << "\n# mu        " << mu
         << "\n# mu_sigma  " << mu_sigma
         << "\n# N         " << N
         << "\n# AR(p)     " << params.size()
         << "\n# sigma2eps " << sigma2e[0]
         << "\n# sigma2x   " << gain[0]*sigma2e[0]
         << "\n# submean   " << subtract_mean
         << "\n# T0        " << T0
         << "\n# window_T0 " << window_T0
         << noboolalpha
         << showpos                               // Line up signs
         << '\n'             << real(1)           // Leading one coefficient
         << '\n';
    copy(params.begin(), params.end(), ostream_iterator<real>(cout,"\n"));
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

// Argument checking logic based upon optionparser.h example routines
option::ArgStatus Arg::NonNegative(const option::Option& opt, bool msg)
{
    char *p = 0;
    if (opt.arg) {
        double v = strtol(opt.arg, &p, 10);
        if (p != opt.arg && !*p && v >= 0) {
            return option::ARG_OK;
        }
    }

    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                    << " requires a nonnegative numeric argument\n";
    }
    return option::ARG_ILLEGAL;
}
