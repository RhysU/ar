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
    "\tC = arcov (d1, arsel1, d2, arsel2, absrho)\n"
    "\tC = arcov (d1, arsel1,             absrho)\n"
    "\tFind a correlation matrix for signal collections processed by arsel.\n"
    "\t\n"
    "\tUse raw samples and information computed by arsel(...) along with\n"
    "\tar::decorrelation_time to build an effective covariance matrix\n"
    "\tfor two vectors of signals.  Input d1 should raw samples and\n"
    "\tarsel1 should have been produced by 'arsel1 = arsel(d1, ...)'.\n"
    "\tSimilarly for d2 and arsel2.  Inputs d1 and d2 must have the same\n"
    "\tdimensions.  A structure is returned which each field contains a\n"
    "\tmatrix indexed by the corresponding row in d1 and d2 or a single,\n"
    "\tglobally descriptive scalar.\n"
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
    Matrix d1, d2;
    Octave_map arsel1, arsel2;
    switch (args.length())
    {
        case 5:  absrho = args(4).bool_value();    // Falls through
        case 4:  arsel2 = args(3).map_value();
                 d2     = args(2).matrix_value();
                 arsel1 = args(1).map_value();
                 d1     = args(0).matrix_value();
                 argsok = !error_state;
                 break;
        case 3:  absrho = args(2).bool_value();    // Falls through
        case 2:  arsel2 = args(1).map_value();
                 d2     = args(0).matrix_value();
                 arsel1 = args(1).map_value();
                 d1     = args(0).matrix_value();
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

    // Determine problem size based on d1 and d2
    if (!(d1.dims() == d2.dims()))
    {
        error("arcov: d1 and d2 must have same dimensions");
        return octave_value();
    }
    const octave_idx_type M = d1.rows();  // Number of signals
    const octave_idx_type N = d1.cols();  // Samples per signal

    // Unpack required fields from arsel1 into typesafe instances
    Cell&        AR1      = arsel1.contents("AR");
    Cell&        autocor1 = arsel1.contents("autocor");
    ColumnVector gain1    = arsel1.contents("gain")(0).column_vector_value();
    bool         submean1 = arsel1.contents("submean")(0).bool_value();

    // Unpack required fields from arsel1 into typesafe instances
    Cell&        AR2      = arsel2.contents("AR");
    Cell&        autocor2 = arsel2.contents("autocor");
    ColumnVector gain2    = arsel2.contents("gain")(0).column_vector_value();
    bool         submean2 = arsel2.contents("submean")(0).bool_value();

    // Check that the unpack logic worked as expected
    if (error_state)
    {
        error("arcov: arsel1 or arsel2 lacked fields provided by arsel(...)");
        return octave_value();
    }

    // Prepare storage to be returned to the caller
    Matrix _T0     (dim_vector(M, M));
    Matrix _eff_N  (dim_vector(M, M));
    Matrix _eff_cov(dim_vector(M, M));

    // Compute upper triangular portion of _T0, _eff_N, and _eff_cov
    // and store the result into both the upper and lower triangular storage
    for (octave_idx_type j = 0; j < M; ++j)
    {
        // Prepare iterators into the raw data for signal d1(j)
        ar::strided_adaptor<double*> s1_begin(&d1(j,0), M);
        ar::strided_adaptor<double*> s1_end  (&d1(j,N), M);

        // Prepare an iterator over the autocorrelation function for d1(j)
        RowVector AR1j = AR1(j).row_vector_value();
        ar::predictor<double> p1 = ar::autocorrelation(
                AR1j.fortran_vec(), AR1j.fortran_vec() + AR1j.length(),
                gain1(j), autocor1(j).row_vector_value().fortran_vec());

        for (octave_idx_type i = j; i < M; ++i)
        {
            // Prepare an iterator over the autocorrelation function for d2(i)
            RowVector AR2i = AR2(i).row_vector_value();
            ar::predictor<double> p2 = ar::autocorrelation(
                    AR2i.fortran_vec(), AR2i.fortran_vec() + AR2i.length(),
                    gain2(i), autocor2(i).row_vector_value().fortran_vec());

            // Compute the decorrelation time of d1(j) against d2(i)
            const double T0 = ar::decorrelation_time(N, p1, p2, absrho);

            // Compute the effective covariance given the decorrelation time
            ar::strided_adaptor<double*> s2_begin(&d2(i,0), M);
            double mu1, mu2, ncovar;
            welford_ncovariance(s1_begin, s1_end, s2_begin, mu1, mu2, ncovar);
            if (!submean1 && !submean2) ncovar += N*mu1*mu2;
            const double eff_cov = ncovar / (N - T0);

            // Save the findings into the (symmetric) result storage
            _T0     (i, j) = _T0     (j, i) = T0;
            _eff_N  (i, j) = _eff_N  (j, i) = N / T0;
            _eff_cov(i, j) = _eff_cov(j, i) = eff_cov;
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
