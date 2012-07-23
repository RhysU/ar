Burg's method in header-only C++
================================

Description
-----------

This is a precision-agnostic implementation of Burg's recursive method for
estimating autoregressive model parameters based on the presentation in
[Collomb2009].   Many usability-related extensions have been added to permit
simply obtaining autocorrelation information from the resulting estimated
model.  The expressions differ slightly from those presented in [Press2007],
whose source cannot be distributed due to licensing restrictions.  However, the
results match Numerical Recipes code to floating point error against benchmark
data by [Bourke1998]::

	./test test1.coeff test1.dat

	         Coefficient   NumericalRecipes       burgs_method           PercentDiff
	         -----------   ----------------       ------------           -----------
	                 1.4    1.3722815374162    1.3722815374161  -7.4431168445089e-12
	                -0.7  -0.69066992983051   -0.6906699298305  -1.9289498094656e-12
	                0.04  0.034482949286973  0.034482949287216   7.0461579987592e-10
	                 0.7   0.74323368158906   0.74323368158868  -5.1385819821096e-11
	                -0.5  -0.46165133140387  -0.46165133140364  -4.9240227164672e-11

	  Mean^2 Discrepancy    1.0600415436226    1.0600415436232   6.0934192590181e-11

	./test test2.coeff test2.dat

	         Coefficient   NumericalRecipes       burgs_method           PercentDiff
	         -----------   ----------------       ------------           -----------
	               0.677   0.65876631854655   0.65876631854655  -1.1797144076101e-13
	               0.175   0.20659300731471   0.20659300731471  -4.2991698085266e-13
	               0.297   0.28423342541501   0.28423342541501  -3.5154230038347e-13
	               0.006  0.046456925303432  0.046456925303431  -3.2859614554535e-13
	              -0.114  -0.12126498977465  -0.12126498977465   1.2588601720054e-13
	              -0.083 -0.088294332237173 -0.088294332237173  -2.9863521436483e-13
	              -0.025  -0.05546648074175  -0.05546648074175  -3.2526174202682e-13

	  Mean^2 Discrepancy   0.97858693598551   0.97858693598555   3.8006302717997e-12

	./test test3.coeff test3.dat

	         Coefficient   NumericalRecipes       burgs_method           PercentDiff
	         -----------   ----------------       ------------           -----------
	                1.02    1.0251761581124    1.0251761581124   3.4654665451271e-13
	               -0.53  -0.52577027224279  -0.52577027224279   4.8567085121517e-13

	  Mean^2 Discrepancy    1.0577040559129    1.0577040559129  -6.5078532262362e-13

The implementation also permits extracting a sequence of AR(p) models for p
from one up to some maximum order::

	Order       RMS/N Coefficients
	-----       ----- ------------
	    0    0.001172     0.9995
	    1  5.0336e-06     1.9969   -0.99785
	    2  2.7649e-08      2.992    -2.9892    0.99725
	    3  6.0467e-11     3.9881    -5.9751     3.9859   -0.99891
	    4  6.6639e-14     4.9865    -9.9589     9.9578    -4.9848    0.99945
	    5  1.2883e-16      5.985    -14.939     19.906    -14.934     5.9811   -0.99903
	    6  4.8869e-17      6.772    -19.651     31.672    -30.617      17.75    -5.7142    0.78783
	    7  4.8602e-17     6.7138    -19.229      30.36    -28.354      15.41    -4.2621    0.28742   0.073894

Also included is a Toeplitz linear equation solver for a single right hand side
using O(3m^2) operations.  The algorithm is [Zohar1974]'s improvement of
[Trench1967]'s work::

	Topmost row of Toeplitz matrix is:
		1 2 3 5 7 11 13 17
	Leftmost column of Toeplitz matrix is:
		1 2 4 8 16 32 64 128
	Right hand size data is:
		1 2 3 4 5 6 7 8
	Expected solution is:
		-0.62963 0.148148 3.55556 -1.66667 0 -2 -1 2
	Solution computed by zohar_linear_solve is:
		-0.62963 0.148148 3.55556 -1.66667 7.10543e-15 -2 -1 2
	Term-by-term errors are:
		5.55112e-16 1.04361e-14 -2.70894e-14 9.99201e-15 -7.10543e-15 4.44089e-15 1.26565e-14 -9.32587e-15
	Sum of the absolute errors is:
		8.16014e-14
See [Bunch1985] for a discussion of the stability of Trench-like algorithms and
for faster, albeit much more complicated, variants.

Contents
--------

*burg.hpp*
  Standalone header implementing the Burg recursion and the Zohar Toeplitz solver.

*Makefile*
   Try ``make`` followed by ``./example`` or ``make check``.

*example.cpp*
   Extract a sequence of AR(p) models for a sample by [Collomb2009].

*zohar.cpp*
   Solve a nonsymmetric, real-valued Toeplitz set of linear equations.

*test.cpp*
   A test driver for testing ``burg.hpp`` against benchmarks by [Bourke1998].

*test\*.coeff*
*test\*.dat*
   Sample data and exact coefficients from [Bourke1998] used for ``make check``.

*rhoe.coeff*
*rhoe.dat*
   Sample turbulent total energy RMS fluctuation data and optimal coefficients
   found by automatically by ARMASA [Broersen2002].

*WuleYalker.tex*
   A derivation of some equations closely connected with the Yule--Walker
   system.  Solving these permits recovering autocorrelations from process
   parameters.

*Collomb_Burg.pdf*
   For posterity, a copy of [Collomb2009].

TODO
----

1. Implement the ``AIC``, ``AICc``, and ``AKICc`` model selection criteria
   following [Seghouane2004].  Add these as standalone routines but also
   incorporate them into the class.

2. Add a class to encapsulate the sequence of AR(p) models produced.  Include
   prediction both with and without noise and prediction error computations
   against known data.

3. Use the AR polynomial (e.g. [Broersen2006] equation 4.36) to obtain the
   autocorrelation for arbitrary lags ([Broersen2006] equation 4.52).

4. To find the lag 1, ..., p-1 autocorrelation boundary conditions given only
   process parameters, implement a ``Wule-Yalker`` solver based on the
   WuleYalker.tex write up using the Toeplitz-plus-Hankel solver approach due
   to [Merchant1982] which employs [Akaike1973].  The double Levinson recursion
   discussed by [Broersen2006] section 5.4 appears to be too numerically
   unstable to use in practice without requiring O(n^2) memory.

5. Implement the Ibrahim Optimum Tapered Burg as described by [Campbell1993]
   based on work in [Ibrahim1987a], [Ibrahim1987b], and [Ibrahim1989].

References
----------

-- [Akaike1973]    Akaike, Hirotugu. "Block Toeplitz Matrix Inversion." SIAM Journal on Applied Mathematics 24 (March 1973): 234-241. http://dx.doi.org/10.1137/0124024

-- [Bourke1998]    Bourke, Paul. AutoRegression Analysis, November 1998. http://paulbourke.net/miscellaneous/ar/

-- [Box2008]       Box, George E. P., Gwilym M. Jenkins, and Gregory C. Reinsel. Time Series Analysis : Forecasting and Control. 4 edition. John Wiley, June 2008.

-- [Broersen2002]  Broersen, P. M. T. "Automatic spectral analysis with time series models." IEEE Transactions on Instrumentation and Measurement 51 (April 2002): 211-216. http://dx.doi.org/10.1109/19.997814

-- [Broersen2006]  Broersen, P. M. T. Automatic autocorrelation and spectral analysis. Springer, 2006. http://dx.doi.org/10.1007/1-84628-329-9

-- [Bunch1985]     Bunch, James R. "Stability of Methods for Solving Toeplitz Systems of Equations." SIAM Journal on Scientific and Statistical Computing 6 (1985): 349-364. http://dx.doi.org/10.1137/0906025

-- [Campbell1993]  Campbell, W. and D. N. Swingler. "Frequency estimation performance of several weighted Burg algorithms." IEEE Transactions on Signal Processing 41 (March 1993): 1237-1247. http://dx.doi.org/10.1109/78.205726

-- [Collomb2009]   Cedrick Collomb. Burg's method, algorithm, and recursion, November 2009. http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf

-- [Ibrahim1987a]  Ibrahim, M. K. "Improvement in the speed of the data-adaptive weighted Burg technique." IEEE Transactions on Acoustics, Speech, and Signal Processing 35 (October 1987): 1474–1476. http://dx.doi.org/10.1109/TASSP.1987.1165046

-- [Ibrahim1987b]  Ibrahim, M. K. "On line splitting in the optimum tapered Burg algorithm." IEEE Transactions on Acoustics, Speech, and Signal Processing 35 (October 1987): 1476–1479. http://dx.doi.org/10.1109/TASSP.1987.1165047

-- [Ibrahim1989]   Ibrahim, M. K. "Correction to 'Improvement in the speed of the data-adaptive weighted Burg technique'." IEEE Transactions on Acoustics, Speech, and Signal Processing 37 (1989): 128. http://dx.doi.org/10.1109/29.17511

-- [Merchant1982]  Merchant, G. and T. Parks. "Efficient solution of a Toeplitz-plus-Hankel coefficient matrix system of equations." IEEE Transactions on Acoustics, Speech, and Signal Processing 30 (February 1982): 40-44. http://dx.doi.org/10.1109/TASSP.1982.1163845

-- [Press2007]     Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. Numerical recipes : The Art of Scientific Computing. Third edition. Cambridge University Press, September 2007.

-- [Seghouane2004] Seghouane, A. K. and M. Bekara. "A Small Sample Model Selection Criterion Based on Kullback's Symmetric Divergence." IEEE Transactions on Signal Processing 52 (December 2004): 3314-3323. http://dx.doi.org/10.1109/TSP.2004.837416

-- [Trench1967]    Trench, William F. Weighting coefficients for the prediction of stationary time series from the finite past. SIAM J. Appl. Math. 15, 6 (Nov. 1967), 1502-1510. http://www.jstor.org/stable/2099503

-- [Zohar1974]     Zohar, Shalhav. "The Solution of a Toeplitz Set of Linear Equations." J. ACM 21 (April 1974): 272-276. http://dx.doi.org/10.1145/321812.321822
