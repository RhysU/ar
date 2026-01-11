# Copyright (C) 2025 Rhys Ulerich
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""Test the AR Python module using the ar.arsel function."""

import sys
import ar
import numpy as np


def main():
    """Test ar.arsel() against known test data."""

    # Load expected coefficients and test data
    expected = [float(x) for x in open("test0.coeff").read().split()]
    data = [float(x) for x in open("test0.dat").read().split()]

    # Fit AR model using ar.arsel
    result = ar.arsel(
        data,
        False,
        True,
        "CIC",
        len(expected),
        len(expected),
    )
    estimated = result.AR[0][1:]

    if len(estimated) != len(expected):
        raise ValueError(
            "Expected {} coefficients, got {}".format(
                len(expected), len(estimated)
            )
        )

    # Compare coefficients
    print(
        "{:>22s} {:>22s} {:>22s}".format(
            "Nominal Value", "arsel Result", "Percent Diff"
        )
    )
    print(
        "{:>22s} {:>22s} {:>22s}".format(
            "-------------", "------------", "------------"
        )
    )

    max_error = 0.0
    for exp, est in zip(expected, estimated):
        diff = (
            (est - exp) / exp * 100.0
            if exp != 0
            else 0.0 if est == 0.0 else float("inf")
        )
        max_error = max(max_error, abs(diff))
        print("{:22.14g} {:22.14g} {:22.14g}".format(exp, est, diff))

    print()
    print("{:>22s} {:22.14g}".format("mean of data:", result.mu[0]))
    print(
        "{:>22s} {:22.14g}".format(
            "\\sigma^2_\\epsilon:", result.sigma2eps[0]
        )
    )
    print("{:>22s} {:22.14g}".format("signal gain:", result.gain[0]))
    print("{:>22s} {:22.14g}".format("\\sigma^2_x:", result.sigma2x[0]))

    if max_error > 0.01:
        raise ValueError(
            "Maximum error {:.6g}% exceeds tolerance 0.01%".format(max_error)
        )

    print("\nPython AR test PASSED (max error: {:.6g}%)".format(max_error))


if __name__ == "__main__":
    main()
