// Copyright (C) 2012, 2013 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Generate (t, x) data for the AR(6) process defined in section 4.1
 * of https://arxiv.org/abs/1802.01056v1.  Intended to be piped into
 * another program, e.g. arsel.cpp, for analysis.
 */

#include <sys/time.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "ar.hpp"
#include "optionparser.h"

// Declarations for argument checking logic based upon optionparser.h examples
struct Arg : public option::Arg
{
    static option::ArgStatus IntNonNeg(const option::Option& opt, bool msg);
    static option::ArgStatus Double   (const option::Option& opt, bool msg);
};

// Command line argument declarations for optionparser.h usage
enum OptionIndex {
    UNKNOWN, BURN, SEED, TFINAL, HELP
};
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "", Arg::None,
     "Usage: ar6 [OPTION]...\n"
     "Output tab-separated (t, x) data for an AR(6) problem\n"
     "Details in section 4.1 of https://arxiv.org/abs/1802.01056v1.\n"
     "Individual columns may be extracted by piping to cut(1) utility.\n"
     "\n"
     "Options:" },
    {0,0,"","",Arg::None,0}, // table break
    {BURN,     0, "b", "burn",     Arg::IntNonNeg,
     "  -b \t--burn=BURN  \t \"Burn-in\" for 0 <= t < BURN defaulting to 500" },
    {SEED,     0, "s", "seed",     Arg::Double,
     "  -s \t--seed=SEED  \t Random seed defaulting to gettimeofday tv_usec"  },
    {TFINAL,   0, "t", "tfinal",   Arg::IntNonNeg,
     "  -t \t--tfinal=T   \t Advance time until t >= T defaulting to 3000"    },
    {HELP,     0, "h", "help",     Arg::None,
     "  -h \t--help       \t Display this help message and immediately exit"  },
    {0,0,0,0,0,0}
};

// Box-Muller logic adjusted from
// https://github.com/rflynn/c/blob/master/rand-normal-distribution.c
static double rand01() {
    const double x = (double) random() / RAND_MAX;
    const double y = (double) random() / RAND_MAX;
    const double z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
    return z;
}

// Scale to obtain variance 0.1.
static double rand0point1() {
    return sqrt(0.1) * rand01();
}

int main(int argc, char *argv[])
{
    using namespace std;

    long t      = 0;
    long burn   = 500;
    long tfinal = 3000;

    {
        option::Stats stats(usage, argc-(argc>0), argv+(argc>0));

        vector<option::Option> opts(stats.options_max + stats.buffer_max);

        option::Parser parse(usage, argc-(argc>0), argv+(argc>0),
                             &opts[0], &opts[stats.options_max]);

        if (parse.error() || opts[UNKNOWN]) {
            for (option::Option* o = opts[UNKNOWN]; o; o = o->next()) {
                cerr << "Unknown option: " << o->name << "\n";
            }
            return EXIT_FAILURE;
        }

        if (opts[HELP]) {
            printUsage(cout, usage);
            return EXIT_SUCCESS;
        }

        // SEED must come before any other randomly-generated option
        if (opts[SEED]) {
            srandom((unsigned) strtol(opts[SEED].last()->arg, NULL, 10));
        } else {
            struct timeval tv;
            struct timezone tz;
            gettimeofday(&tv, &tz);
            srand((unsigned) tv.tv_usec);
        }

        // Parse remaining options
        if (opts[BURN  ]) burn   = strtol(opts[BURN  ].last()->arg, NULL, 10);
        if (opts[TFINAL]) tfinal = strtol(opts[TFINAL].last()->arg, NULL, 10);

        // Warn whenever burn >= tfinal
        if (burn >= tfinal) {
            cerr << "Warning: burn >= tfinal suppresses output ("
                 << burn << " >= " << tfinal << ") \n";
        }
    }

    // Initialize the AR process parameters
    std::vector<double> params;
    params.push_back(+3.1378);
    params.push_back(-3.9789);
    params.push_back(+2.6788);
    params.push_back(-1.0401);
    params.push_back(+0.2139);
    params.push_back(-0.0133);
    ar::predictor<double> p(params.begin(), params.end(), rand0point1);

    std::vector<double> initial;
    for (size_t i = 0; i < params.size(); ++i) {
        initial.push_back(rand0point1());
    }
    p.initial_conditions(initial.begin());

    // Discard 0 <= t < burn
    for (; t < burn; ++t, ++p) {
        // NOP
    }

    // Output (t, x, y, z) during burn < t <= tfinal
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos;
    for (; t < tfinal; ++t, ++p) {
        cout << t << '\t' << *p << '\n';
    }

    return EXIT_SUCCESS;
}

option::ArgStatus Arg::IntNonNeg(const option::Option& opt, bool msg)
{
    char *p = 0;
    if (opt.arg) {
        double v = strtol(opt.arg, &p, 10);
        if (p != opt.arg && !*p && v >= 0) return option::ARG_OK;
    }
    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                    << " requires a nonnegative integer argument\n";
    }
    return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Double(const option::Option& opt, bool msg)
{
    char *p = 0;
    if (opt.arg) {
        double v = strtod(opt.arg, &p); (void) v; // Discard
        if (p != opt.arg && !*p) return option::ARG_OK;
    }
    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                   << " requires a floating point argument\n";
    }
    return option::ARG_ILLEGAL;
}
