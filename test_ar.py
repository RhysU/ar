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

    # Fit AR model using ar.arsel with POSITIONAL arguments
    result_positional = ar.arsel(
        data,
        False,
        True,
        "CIC",
        len(expected),
        len(expected),
    )

    # Fit AR model using ar.arsel with KEYWORD arguments
    result_keyword = ar.arsel(
        data=data,
        submean=False,
        absrho=True,
        criterion="CIC",
        minorder=len(expected),
        maxorder=len(expected),
    )

    # Verify that positional and keyword calls produce identical results
    print("Verifying positional vs keyword arguments produce identical results...")

    def values_equal(a, b):
        """Check if two values are equal, handling NaN and inf correctly."""
        # If both are NaN, consider them equal
        if np.isnan(a) and np.isnan(b):
            return True
        # If both are inf with same sign, consider them equal
        if np.isinf(a) and np.isinf(b) and (a > 0) == (b > 0):
            return True
        # Otherwise use normal equality
        return a == b

    assert len(result_positional.AR[0]) == len(result_keyword.AR[0]), \
        "AR length mismatch between positional and keyword calls"
    for i, (pos_val, kw_val) in enumerate(zip(result_positional.AR[0], result_keyword.AR[0])):
        assert values_equal(pos_val, kw_val), \
            f"AR[0][{i}] mismatch: positional={pos_val}, keyword={kw_val}"
    assert values_equal(result_positional.mu[0], result_keyword.mu[0]), \
        f"mu mismatch: positional={result_positional.mu[0]}, keyword={result_keyword.mu[0]}"
    assert values_equal(result_positional.sigma2eps[0], result_keyword.sigma2eps[0]), \
        f"sigma2eps mismatch: positional={result_positional.sigma2eps[0]}, keyword={result_keyword.sigma2eps[0]}"
    assert values_equal(result_positional.gain[0], result_keyword.gain[0]), \
        f"gain mismatch: positional={result_positional.gain[0]}, keyword={result_keyword.gain[0]}"
    assert values_equal(result_positional.sigma2x[0], result_keyword.sigma2x[0]), \
        f"sigma2x mismatch: positional={result_positional.sigma2x[0]}, keyword={result_keyword.sigma2x[0]}"
    print("✓ Positional and keyword arguments produce IDENTICAL results\n")

    # Use positional result for the rest of the test
    result = result_positional
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
