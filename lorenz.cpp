// Copyright (C) 2012 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Generate (t, x, y, z) data from the Lorenz attractor on standard output.
 * Intended to be piped into another program, e.g. arsel.cpp, for analysis.
 */

#include "optionparser.h"

#include <sys/time.h>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

// Declarations for argument checking logic based upon optionparser.h examples
struct Arg : public option::Arg
{
    static option::ArgStatus Int      (const option::Option& opt, bool msg);
    static option::ArgStatus IntNonNeg(const option::Option& opt, bool msg);
    static option::ArgStatus Double   (const option::Option& opt, bool msg);
    static option::ArgStatus DoublePos(const option::Option& opt, bool msg);
};

// Command line argument declarations for optionparser.h usage
enum OptionIndex {
    UNKNOWN, BETA, BURN, DT, EVERY, MORE, SCHEME, RHO, SEED, SIGMA,
    TFINAL, INITX, INITY, INITZ, HELP
};
enum SchemeType {
    RK1, RK2, RK3
};
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "", Arg::None,
     "Usage: lorenz [OPTION]...\n"
     "Output tab-separated (t, x, y, z) data from the Lorenz equations:\n"
     "  d/dt {x,y,z} = {sigma*(y - x), x*(rho - z) - y, x*y - beta*z}\n"
     "Advance per Gottlieb and Shu 1998 (doi:10.1090/S0025-5718-98-00913-2).\n"
     "Individual columns may be extracted by piping to cut(1) utility.\n"
     "\n"
     "Options:" },
    {0,0,"","",Arg::None,0}, // table break
    {BETA,     0, "B", "beta",     Arg::Double,
     "  -B \t--beta=BETA  \t Beta coefficient defaulting to 8/3"              },
    {BURN,     0, "b", "burn",     Arg::DoublePos,
     "  -b \t--burn=BURN  \t \"Burn-in\" for 0 <= t < BURN defaulting to 500" },
    {DT,       0, "d", "dt",       Arg::DoublePos,
     "  -d \t--dt=DT      \t Fixed time step size defaulting to 0.01"         },
    {EVERY,    0, "e", "every",    Arg::IntNonNeg,
     "  -e \t--every=N    \t Output every Nth step defaulting to 1"           },
    {MORE,      0, "m", "more",    Arg::None,
     "  -m \t--more       \t Output more columns (xx, xy, xz, yy, yz, zz)"    },
    {RHO,      0, "R", "rho",      Arg::Double,
     "  -R \t--rho=RHO    \t Rho coefficient defaulting to 28"                },
    {SEED,     0, "s", "seed",     Arg::Double,
     "  -s \t--seed=SEED  \t Random seed defaulting to gettimeofday tv_usec"  },
    {SIGMA,    0, "S", "sigma",    Arg::Double,
     "  -S \t--sigma=SIGMA\t Sigma coefficient defaulting to 10"              },
    {TFINAL,   0, "t", "tfinal",   Arg::DoublePos,
     "  -t \t--tfinal=T   \t Advance time until t >= T defaulting to 3000"    },
    {INITX,    0, "x", "initx",    Arg::Double,
     "  -x \t--initx=POS  \t Initial x at t = 0 defaulting to rand()/RAND_MAX"},
    {INITY,    0, "y", "inity",    Arg::Double,
     "  -y \t--inity=POS  \t Initial y at t = 0 defaulting to rand()/RAND_MAX"},
    {INITZ,    0, "z", "initz",    Arg::Double,
     "  -z \t--initz=POS  \t Initial z at t = 0 defaulting to rand()/RAND_MAX"},
    {SCHEME, RK1, "1", "euler",    Arg::None,
     "  -1 \t--euler      \t Advance with 1st order Forward Euler scheme"     },
    {SCHEME, RK2, "2", "rk2",      Arg::None,
     "  -2 \t--rk2        \t Advance with 2nd order TVD Runge--Kutta scheme"  },
    {SCHEME, RK3, "3", "rk3",      Arg::None,
     "  -3 \t--rk3        \t Advance with default 3rd order TVD Runge--Kutta" },
    {HELP,     0, "h", "help",     Arg::None,
     "  -h \t--help       \t Display this help message and immediately exit"  },
    {0,0,0,0,0,0}
};

// Enable strict IEEE behavior for reproducibility
// Inlining should help offset at least some of the performance loss
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) > 40305
# pragma GCC push_options
#elif _MSC_VER > 1400
# pragma float_control(push)
# pragma float_control(precise, on)
#else
# warning Pushing precise floating point behavior unimplemented for compiler
#endif

/** Computation of the Lorenz equation right hand side. */
static inline void lorenz(
    const double  beta, const double rho,   const double sigma,
    const double  x,    const double y,     const double z,
            double &dxdt,       double &dydt,       double &dzdt)
{
     dxdt = sigma*(y - x);
     dydt = x*(rho - z) - y;
     dzdt = x*y - beta*z;
}

/** Advance by \c dt using one step of 1st order Forward Euler. */
static void euler(
        const double dt,
        const double beta, const double rho, const double sigma,
              double &t,
              double &x,         double &y,        double &z)
{
    double dxdt, dydt, dzdt;
    lorenz(beta, rho, sigma, x, y, z, dxdt, dydt, dzdt);
    x += dt*dxdt;
    y += dt*dydt;
    z += dt*dzdt;
    t += dt;
}

/**
 * Advance by \c dt using one step of 2nd order TVD Runge--Kutta.
 * This optimal scheme appears in Proposition 3.1 of Gottlieb and Shu 1998.
 */
static void tvd_rk2(
        const double dt,
        const double beta, const double rho, const double sigma,
              double &t,
              double &x,         double &y,        double &z)
{
    double dxdt, dydt, dzdt;

    lorenz(beta, rho, sigma, x, y, z, dxdt, dydt, dzdt);
    double u1x = x + dt*dxdt;
    double u1y = y + dt*dydt;
    double u1z = z + dt*dzdt;

    lorenz(beta, rho, sigma, u1x, u1y, u1z, dxdt, dydt, dzdt);
    x = (x + u1x + dt*dxdt)/2;
    y = (y + u1y + dt*dydt)/2;
    z = (z + u1z + dt*dzdt)/2;
    t += dt;
}

/**
 * Advance by \c dt using one step of 3rd order TVD Runge--Kutta.
 * The optimal scheme appears in Proposition 3.2 of Gottlieb and Shu 1998.
 */
static void tvd_rk3(
        const double dt,
        const double beta, const double rho, const double sigma,
              double &t,
              double &x,         double &y,        double &z)
{
    double dxdt, dydt, dzdt;

    lorenz(beta, rho, sigma, x, y, z, dxdt, dydt, dzdt);
    double u1x = x + dt*dxdt;
    double u1y = y + dt*dydt;
    double u1z = z + dt*dzdt;

    lorenz(beta, rho, sigma, u1x, u1y, u1z, dxdt, dydt, dzdt);
    double u2x = (3*x + u1x + dt*dxdt)/4;
    double u2y = (3*y + u1y + dt*dydt)/4;
    double u2z = (3*z + u1z + dt*dzdt)/4;

    lorenz(beta, rho, sigma, u2x, u2y, u2z, dxdt, dydt, dzdt);
    x = (x + 2*u2x + 2*dt*dxdt)/3;
    y = (y + 2*u2y + 2*dt*dydt)/3;
    z = (z + 2*u2z + 2*dt*dzdt)/3;
    t += dt;
}

// Disable strict IEEE behavior for remainder of file
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) > 40305
# pragma GCC pop_options
#elif _MSC_VER > 1400
# pragma float_control(pop)
#else
# warning Popping precise floating point behavior unimplemented for compiler
#endif

// TODO Simulation time t drifts slowly due to accumulation error

int main(int argc, char *argv[])
{
    using namespace std;

    double     beta   = 8./3;
    double     burn   = 500;
    double     dt     = 0.01;
    long       every  = 1;
    double     rho    = 28;
    SchemeType scheme = RK3;
    double     sigma  = 10;
    double     t      = 0;
    double     tfinal = 3000;
    bool       more   = false;
    double     x;
    double     y;
    double     z;
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

        // SEED must come before any other randomly-generated option, e.g. INITX
        if (opts[SEED]) {
            srandom((unsigned) strtol(opts[SEED].last()->arg, NULL, 10));
        } else {
            struct timeval tv;
            struct timezone tz;
            gettimeofday(&tv, &tz);
            srand((unsigned) tv.tv_usec);
        }

        // INITX, INITY, INITZ either specified or randomized
        x = opts[INITX] ? strtod(opts[INITX].last()->arg, NULL)
                        : (double)rand()/(double)RAND_MAX;
        y = opts[INITY] ? strtod(opts[INITY].last()->arg, NULL)
                        : (double)rand()/(double)RAND_MAX;
        z = opts[INITZ] ? strtod(opts[INITZ].last()->arg, NULL)
                        : (double)rand()/(double)RAND_MAX;

        // Parse remaining options
        if (opts[BETA  ]) beta   = strtod(opts[BETA  ].last()->arg, NULL);
        if (opts[BURN  ]) burn   = strtod(opts[BURN  ].last()->arg, NULL);
        if (opts[DT    ]) dt     = strtod(opts[DT    ].last()->arg, NULL);
        if (opts[EVERY ]) every  = strtol(opts[EVERY ].last()->arg, NULL, 10);
        if (opts[RHO   ]) rho    = strtod(opts[RHO   ].last()->arg, NULL);
        if (opts[SIGMA ]) sigma  = strtod(opts[SIGMA ].last()->arg, NULL);
        if (opts[TFINAL]) tfinal = strtod(opts[TFINAL].last()->arg, NULL);
        if (opts[MORE  ]) more   = true;

        // Warn whenever burn >= tfinal
        if (burn >= tfinal) {
            cerr << "Warning: burn >= tfinal suppresses output ("
                 << burn << " >= " << tfinal << ") \n";
        }
    }

    // Discard 0 <= t < burn
    switch (scheme) {  // Switch designed to avoid spurious jumps
        default:  cerr << "Sanity failure: unknown scheme at "
                       << __FILE__ << ":" << __LINE__ << '\n';
                  return EXIT_FAILURE;
        case RK1: while (t < burn) euler  (dt, beta, rho, sigma, t, x, y, z);
                  break;
        case RK2: while (t < burn) tvd_rk2(dt, beta, rho, sigma, t, x, y, z);
                  break;
        case RK3: while (t < burn) tvd_rk3(dt, beta, rho, sigma, t, x, y, z);
                  break;
    }

    // Output (t, x, y, z) during burn < t <= tfinal at periodic intervals
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos;
    do {

        cout << t << '\t' << x << '\t' << y << '\t' << z;
        if (more) {
            cout << '\t' << x*x << '\t' << x*y << '\t' << x*z
                                << '\t' << y*y << '\t' << y*z
                                               << '\t' << z*z;
        }
        cout << '\n';

        switch (scheme) {
        default:
            cerr << "Sanity failure: unknown scheme at "
                 << __FILE__ << ":" << __LINE__ << '\n';
            return EXIT_FAILURE;
        case RK1:
            for (long i = 0; i < every; ++i) {
                euler  (dt, beta, rho, sigma, t, x, y, z);
            }
            break;
        case RK2:
            for (long i = 0; i < every; ++i) {
                tvd_rk2(dt, beta, rho, sigma, t, x, y, z);
            }
            break;
        case RK3:
            for (long i = 0; i < every; ++i) {
                tvd_rk3(dt, beta, rho, sigma, t, x, y, z);
            }
            break;
        }

    } while (t < tfinal);

    return EXIT_SUCCESS;
}

option::ArgStatus Arg::Int(const option::Option& opt, bool msg)
{
    char *p = 0;
    if (opt.arg) {
        long v = strtol(opt.arg, &p, 10); (void) v; // Discard
        if (p != opt.arg && !*p) return option::ARG_OK;
    }
    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                    << " requires an integer argument\n";
    }
    return option::ARG_ILLEGAL;
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

option::ArgStatus Arg::DoublePos(const option::Option& opt, bool msg)
{
    char *p = 0;
    if (opt.arg) {
        double v = strtod(opt.arg, &p);
        if (p != opt.arg && !*p && v > 0) return option::ARG_OK;
    }
    if (msg) {
        (std::cerr << "Option ").write(opt.name, opt.namelen)
                   << " requires a positive floating point argument\n";
    }
    return option::ARG_ILLEGAL;
}
