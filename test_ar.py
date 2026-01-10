#!/usr/bin/env python
# Copyright (C) 2025 Rhys Ulerich
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
Test the AR Python module using the ar.arsel function.
Uses only standard library for test infrastructure.
"""

import sys
import os

def main():
    """Test ar.arsel() against known test data."""

    # Try importing required modules
    try:
        import ar
    except ImportError:
        print("ar module not available, skipping Python AR test")
        return 0

    try:
        import numpy as np
    except ImportError:
        print("numpy not available, skipping Python AR test")
        return 0

    # Test data files
    coeff_file = "test0.coeff"
    data_file = "test0.dat"

    # Check if test files exist
    if not os.path.exists(coeff_file):
        print("Test file {} not found".format(coeff_file))
        return 1
    if not os.path.exists(data_file):
        print("Test file {} not found".format(data_file))
        return 1

    # Load expected coefficients
    expected_coeffs = []
    with open(coeff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                expected_coeffs.append(float(line))

    # Load time series data
    data = []
    with open(data_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                data.append(float(line))

    # Convert to numpy array (row vector)
    data_array = np.array([data])

    # Determine model order from number of coefficients
    max_order = len(expected_coeffs)

    # Fit AR model using ar.arsel
    # Use positional parameters matching the C++ test
    # arsel(data, submean, absrho, criterion, minorder, maxorder)
    result = ar.arsel(data_array, False, True, 'CIC', max_order, max_order)

    # Extract estimated coefficients (skip leading 1.0)
    estimated_coeffs = result.AR[0][1:]

    # Verify we got the expected number of coefficients
    if len(estimated_coeffs) != len(expected_coeffs):
        print("ERROR: Expected {} coefficients, got {}".format(
            len(expected_coeffs), len(estimated_coeffs)))
        return 1

    # Compare estimated vs expected coefficients
    print("{:>22s} {:>22s} {:>22s}".format(
        "Nominal Value", "arsel Result", "Percent Diff"))
    print("{:>22s} {:>22s} {:>22s}".format(
        "-------------", "------------", "------------"))

    max_error = 0.0
    for i, (expected, estimated) in enumerate(zip(expected_coeffs, estimated_coeffs)):
        if expected != 0:
            percent_diff = (estimated - expected) / expected * 100.0
        else:
            percent_diff = 0.0 if estimated == 0.0 else float('inf')

        max_error = max(max_error, abs(percent_diff))
        print("{:22.14g} {:22.14g} {:22.14g}".format(
            expected, estimated, percent_diff))

    print()
    print("{:>22s} {:22.14g}".format("mean of data:", result.mu[0]))
    print("{:>22s} {:22.14g}".format("\\sigma^2_\\epsilon:", result.sigma2eps[0]))
    print("{:>22s} {:22.14g}".format("signal gain:", result.gain[0]))
    print("{:>22s} {:22.14g}".format("\\sigma^2_x:", result.sigma2x[0]))

    # Accept up to 0.01% error as reasonable for numerical methods
    tolerance = 0.01
    if max_error > tolerance:
        print("\nERROR: Maximum error {:.6g}% exceeds tolerance {:.6g}%".format(
            max_error, tolerance))
        return 1

    print("\nPython AR test PASSED (max error: {:.6g}%)".format(max_error))
    return 0

if __name__ == '__main__':
    sys.exit(main())
