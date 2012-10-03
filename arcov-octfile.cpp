// Copyright (C) 2012 Rhys Ulerich
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "ar.hpp"

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <list>
#include <vector>

#include <octave/oct.h>
#include <octave/oct-map.h>
#include <octave/ov-struct.h>
#include <octave/Cell.h>

/** @file
 * A GNU Octave function finding correlation matrices from arsel(...) output.
 */

// Compile-time defaults in the code also appearing in the help message
#define DEFAULT_ABSRHO true
#define STRINGIFY(x) STRINGIFY_HELPER(x)
#define STRINGIFY_HELPER(x) #x

DEFUN_DLD(
    arcov, args, nargout,
    "\tC = arcov (arsel1, arsel2, absrho)\n"
    "\tC = arcov (arsel1,         absrho)\n"
    "\tFind a correlation matrix for signal collections processed by arsel.\n"
    "\t\n"
    "\tUse raw samples and information computed by arsel(...) along with\n"
    "\tar::decorrelation_time to build an effective covariance matrix for\n"
    "\ttwo vectors of signals.  Input arsel1 should have been produced\n"
    "\tby 'arsel1 = arsel(d1, ...)' and similarly for input arsel2.\n"
    "\tThe signal sources d1 and d2 must have had the same dimensions.\n"
    "\tA structure is returned which each field contains a matrix indexed by\n"
    "\tthe corresponding row in d1 and d2 or a single, globally descriptive\n"
    "\tscalar.\n"
    "\t\n"
    "\tThe number of samples in d1 and d2 (i.e. the number of rows) is\n"
    "\treturned in field 'N'.  Given the observed autocorrelation structure,\n"
    "\ta decorrelation time 'T0' is computed by ar::decorrelation_time\n"
    "\tand used to estimate the effective signal covariance 'eff_cov'.\n"
    "\tThe number of effectively independent samples is returned in 'eff_N'.\n"
    "\tThe absolute value of the autocorrelation function will be used in\n"
    "\tcomputing the decorrelation times whenever absrho is true.\n"
    "\t\n"
    "\tWhen only d1 and arsel1 are provided, the autocovarance is found.\n"
    "\tWhen omitted, absrho defaults to " STRINGIFY(DEFAULT_ABSRHO) ".\n"
)
{
    // Process incoming positional arguments
    bool argsok = false;
    bool absrho = DEFAULT_ABSRHO;
    Octave_map arsel1, arsel2;
    switch (args.length())
    {
        case 3:  absrho = args(2).bool_value();
                 arsel2 = args(1).map_value();
                 arsel1 = args(0).map_value();
                 argsok = !error_state;
                 break;
        case 2:  if (!args(1).is_bool_scalar())  // Disambiguate cases
                 {
                     arsel2 = args(1).map_value();
                 }
                 else
                 {
                     absrho = args(1).bool_value();
                     arsel2 = args(0).map_value();
                 }
                 arsel1 = args(0).map_value();
                 argsok = !error_state;
                 break;
        case 1:  arsel2 = args(0).map_value();
                 arsel1 = args(0).map_value();
                 argsok = !error_state;
                 break;
        default: argsok = !error_state;
    }

    if (!argsok)
    {
        error("Invalid call to arcov.  Correct usage is: ");
        print_usage();
        return octave_value();
    }

    // Unpack required fields from arsel1 into typesafe instances
    const Cell         AR1      = arsel1.contents("AR"     )(0).cell_value();
    const Matrix       d1       = arsel1.contents("data"   )(0).matrix_value();
    const Cell         autocor1 = arsel1.contents("autocor")(0).cell_value();
    const ColumnVector gain1    = arsel1.contents("gain"   )(0).column_vector_value();
    const bool         submean1 = arsel1.contents("submean")(0).bool_value();

    // Unpack required fields from arsel1 into typesafe instances
    const Cell         AR2      = arsel2.contents("AR"     )(0).cell_value();
    const Matrix       d2       = arsel2.contents("data"   )(0).matrix_value();
    const Cell         autocor2 = arsel2.contents("autocor")(0).cell_value();
    const ColumnVector gain2    = arsel2.contents("gain"   )(0).column_vector_value();
    const bool         submean2 = arsel2.contents("submean")(0).bool_value();

    // Check that the unpack logic worked as expected
    if (error_state)
    {
        error("arcov: arsel1 or arsel2 lacked fields provided by arsel(...)");
        return octave_value();
    }

    // Ensure data sets d1 and d2 are conformant
    if (!(d1.dims() == d2.dims()))
    {
        error("arcov: arsel1.data and arsel2.data must have same dimensions");
        return octave_value();
    }
    const octave_idx_type M = d1.rows();  // Number of signals
    const octave_idx_type N = d1.cols();  // Samples per signal

    // Prepare storage to be returned to the caller
    Matrix _T0     (dim_vector(M, M));
    Matrix _eff_N  (dim_vector(M, M));
    Matrix _eff_cov(dim_vector(M, M));

    // Compute upper triangular portion of _T0, _eff_N, and _eff_cov
    // and store the result into both the upper and lower triangular storage
    for (octave_idx_type j = 0; j < M; ++j)
    {
        // Prepare iterators into the raw data for signal d1(j)
        ar::strided_adaptor<const double*> s1_begin(d1.data() + j + 0*M, M);
        ar::strided_adaptor<const double*> s1_end  (d1.data() + j + N*M, M);

        // Prepare an iterator over the autocorrelation function for d1(j)
        const RowVector AR1j = AR1(j).row_vector_value();
        ar::predictor<double> p1 = ar::autocorrelation(
                AR1j.fortran_vec() + 1, AR1j.fortran_vec() + AR1j.length(),
                gain1(j), autocor1(j).row_vector_value().fortran_vec());

        for (octave_idx_type i = j; i < M; ++i)
        {
            // Prepare an iterator over the autocorrelation function for d2(i)
            const RowVector AR2i = AR2 (i).row_vector_value();
            ar::predictor<double> p2 = ar::autocorrelation(
                    AR2i.fortran_vec() + 1, AR2i.fortran_vec() + AR2i.length(),
                    gain2(i), autocor2(i).row_vector_value().fortran_vec());

            // Compute the decorrelation time of d1(j) against d2(i)
            const double T0 = ar::decorrelation_time(N, p1, p2, absrho);

            // Compute the effective covariance given the decorrelation time
            ar::strided_adaptor<const double*> s2_begin(d2.data() + i + 0*M, M);
            double mu1, mu2, ncovar;
            welford_ncovariance(s1_begin, s1_end, s2_begin, mu1, mu2, ncovar);
            if (!submean1 && !submean2) ncovar += N*mu1*mu2;
            const double eff_cov = ncovar / (N - T0);

            // Save the findings into the (symmetric) result storage
            _T0     (i, j) = _T0     (j, i) = T0;
            _eff_N  (i, j) = _eff_N  (j, i) = N / T0;
            _eff_cov(i, j) = _eff_cov(j, i) = eff_cov;

            // Permit user to interrupt the computations at this time
            OCTAVE_QUIT;
        }
    }

    // Provide no results whenever an error was detected
    if (error_state)
    {
        warning("arcov: error detected; no results returned");
        return octave_value_list();
    }

    // Build map containing return fields
    Octave_map retval;
    retval.assign("absrho",    octave_value(absrho));
    retval.assign("eff_cov",   _eff_cov);
    retval.assign("eff_N",     _eff_N);
    retval.assign("T0",        _T0);

    return octave_value_list(retval);
}
