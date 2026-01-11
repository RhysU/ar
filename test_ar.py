# Copyright (C) 2026 Rhys Ulerich
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""Test the AR Python module using the ar.arsel function."""

import numpy as np
import numpy.testing as npt

import ar


def main():
    """Test ar.arsel() against known test data."""
    # Load required inputs
    with open("rhoe.dat") as file:
        data = [float(line.strip()) for line in file.readlines()]
    with open("rhoe.coeff") as file:
        desired = [float(line.strip()) for line in file.readlines()]

    # Fit AR model using both positional and keyword arguments
    kwargs = {
        "data": data,
        "submean": True,
        "absrho": True,
        "criterion": "CIC",
        "minorder": len(desired),
        "maxorder": len(desired),
    }
    result_positional = ar.arsel(*kwargs.values())
    result_keyword = ar.arsel(**kwargs)

    # Confirm results match regardless of calling convention
    npt.assert_equal(result_positional.N, result_keyword.N)
    npt.assert_equal(result_positional.AR, result_keyword.AR)
    npt.assert_equal(result_positional.mu, result_keyword.mu)
    npt.assert_equal(result_positional.sigma2eps, result_keyword.sigma2eps)
    npt.assert_equal(result_positional.gain, result_keyword.gain)
    npt.assert_equal(result_positional.sigma2x, result_keyword.sigma2x)
    npt.assert_equal(result_positional.autocor, result_keyword.autocor)
    npt.assert_equal(result_positional.data, result_keyword.data)

    # Confirm coefficients sufficiently close to expected results
    actual = result_positional.AR[0][1:]
    npt.assert_allclose(actual=actual, desired=desired)


if __name__ == "__main__":
    main()
