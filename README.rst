Burg's method and autoregressive model selection
================================================

Description
-----------

This is a precision-agnostic, header-only, C++ implementation of Burg's
recursive method for estimating autoregressive model parameters based on the
presentation in [Collomb2009].   Many usability-related extensions have been
added to permit simply obtaining autocorrelation information from the resulting
estimated model.  The expressions differ slightly from those presented in
[Press2007] (whose source cannot be distributed due to licensing restrictions).
Still, the results match Numerical Recipes code to floating point error against
benchmark data by [Bourke1998]::

	./test test1.coeff test1.dat

	         Coefficient   NumericalRecipes        burg_method           PercentDiff
	         -----------   ----------------        -----------           -----------
	                 1.4    1.3722815374162    1.3722815374161  -7.4431168445089e-12
	                -0.7  -0.69066992983051   -0.6906699298305  -1.9289498094656e-12
	                0.04  0.034482949286973  0.034482949287216   7.0461579987592e-10
	                 0.7   0.74323368158906   0.74323368158868  -5.1385819821096e-11
	                -0.5  -0.46165133140387  -0.46165133140364  -4.9240227164672e-11

	  Mean^2 Discrepancy    1.0600415436226    1.0600415436232   6.0934192590181e-11

	./test test2.coeff test2.dat

	         Coefficient   NumericalRecipes        burg_method           PercentDiff
	         -----------   ----------------        -----------           -----------
	               0.677   0.65876631854655   0.65876631854655  -1.1797144076101e-13
	               0.175   0.20659300731471   0.20659300731471  -4.2991698085266e-13
	               0.297   0.28423342541501   0.28423342541501  -3.5154230038347e-13
	               0.006  0.046456925303432  0.046456925303431  -3.2859614554535e-13
	              -0.114  -0.12126498977465  -0.12126498977465   1.2588601720054e-13
	              -0.083 -0.088294332237173 -0.088294332237173  -2.9863521436483e-13
	              -0.025  -0.05546648074175  -0.05546648074175  -3.2526174202682e-13

	  Mean^2 Discrepancy   0.97858693598551   0.97858693598555   3.8006302717997e-12

	./test test3.coeff test3.dat

	         Coefficient   NumericalRecipes        burg_method           PercentDiff
	         -----------   ----------------        -----------           -----------
	                1.02    1.0251761581124    1.0251761581124   3.4654665451271e-13
	               -0.53  -0.52577027224279  -0.52577027224279   4.8567085121517e-13

	  Mean^2 Discrepancy    1.0577040559129    1.0577040559129  -6.5078532262362e-13

The implementation also permits extracting a sequence of AR(p) models for p
from one up to some maximum order::

	AR      RMS/N      Gain Filter Coefficients
	--      -----      ---- -------------------
	 0   5.91e+00  1.00e+00 [ 1  ]
	 1   6.02e-04  9.82e+03 [ 1   -0.9999 ]
	 2   1.91e-05  3.09e+05 [ 1    -1.984    0.984 ]
	 3   3.33e-07  1.78e+07 [ 1    -2.959    2.951  -0.9913 ]
	 4   3.54e-08  1.67e+08 [ 1    -3.896     5.74   -3.789   0.9453 ]
	 5   2.82e-09  2.10e+09 [ 1    -4.803    9.375   -9.295    4.683  -0.9594 ]
	 6   6.26e-10  9.44e+09 [ 1    -5.649     13.5   -17.49    12.95   -5.195   0.8819 ]
	 7   1.15e-10  5.15e+10 [ 1    -6.446     18.2    -29.2    28.76    -17.4    5.987  -0.9037 ]
	
	AIC  selects model order 7 as best
	AICC selects model order 6 as best
	CIC  selects model order 6 as best

A variety of finite sample model selection criteria are implemented following
[Broersen2000].  In particular, the

* generalized information criterion (GIC),
* Akaike information criterion (AIC),
* consistent criterion BIC,
* minimally consistent criterion (MCC),
* asymptotically-corrected Akaike information criterion (AICC),
* finite information criterion (FIC),
* finite sample information criterion (FSIC), and
* combined information criterion (CIC)

are all implemented.  An included sample program called ``arsel`` uses CIC to
select the best model order given data from standard input.  For example,
``arsel --subtract-mean < rhoe.dat`` reproduces results from ARMASA
[Broersen2002] on a turbulence signal::

	# N                   1753
	# AR(p)               6
	# Mean                +0.2095528795620023
	# \sigma^2_\epsilon   +8.3374933107465524e-09
	# Gain                +4249.4034177677795
	# \sigma^2_x          +3.5429372570302398e-05
	-2.6990358158025196
	+2.8771725720036141
	-1.7247889420225018
	+0.75024991289684761
	-0.26867206160585461
	+0.067007468591674405

Also included is a Toeplitz linear equation solver for a single right hand side
using O(3m^2) operations.  This solver is useful for investigating the
correctness and numerical stability of estimated process parameters and
autocorrelation information.  The algorithm is [Zohar1974]'s improvement of
[Trench1967]'s work See [Bunch1985] for a discussion of the stability of
Trench-like algorithms and for faster, albeit much more complicated, variants.

::

	Topmost row of Toeplitz matrix is:
		1 2 3 5 7 11 13 17
	Leftmost column of Toeplitz matrix is:
		1 2 4 8 16 32 64 128
	Right hand side data is:
		1 2 3 4 5 6 7 8
	Expected solution is:
		-0.62963 0.148148 3.55556 -1.66667 0 -2 -1 2
	Solution computed by zohar_linear_solve is:
		-0.62963 0.148148 3.55556 -1.66667 7.10543e-15 -2 -1 2
	Term-by-term errors are:
		5.55112e-16 1.04361e-14 -2.70894e-14 9.99201e-15 -7.10543e-15 4.44089e-15 1.26565e-14 -9.32587e-15
	Sum of the absolute errors is:
		8.16014e-14


Contents
--------

*burg.hpp*
  Standalone header implementing Burg's recursive method, the Zohar Toeplitz
  solver, a variety of finite sample model selection criteria, and algorithmic
  helper routines.

*Makefile*
   Try ``make`` followed by ``make check``.  On Linux, try ``make stress`` to
   examine the implementation's performance when piping in plain text data.

*arsel.cpp*
   Given data on standard input, use Burg's method to compute a hierarchy
   of candidate models and select the best one using CIC.

*example.cpp*
   A test driver extracting a hierarchy of AR(p) models for a sample given by
   [Collomb2009].

*zohar.cpp*
   A test driver solving a nonsymmetric, real-valued Toeplitz set of linear
   equations.

*test.cpp*
   A test driver for testing ``burg.hpp`` against benchmarks by [Bourke1998].

*test\*.coeff*, *test\*.dat*
   Sample data and exact parameters from [Bourke1998] used for ``make check``.

*rhoe.coeff*, *rhoe.dat*
   Sample turbulent total energy RMS fluctuation data and optimal parameters
   found by automatically by ARMASA [Broersen2002].

*WuleYalker.tex*
   A derivation of some equations closely connected with the Yule--Walker
   system.  Solving these permits recovering autocorrelations from process
   parameters.

*FiniteSampleCriteria.tex*
   A catalog of the autoregressive model selection criteria implemented.

*Collomb_Burg.pdf*
   For posterity, a copy of [Collomb2009].

Todo
----

1. Add a class to encapsulate a single AR(p) model.  Include prediction both
   with and without noise and prediction error computations against known data.

2. To find the lag 1, ..., p-1 autocorrelation boundary conditions given only
   process parameters, implement a ``Wule-Yalker`` solver based on the
   WuleYalker.tex write up using the Toeplitz-plus-Hankel solver approach due
   to [Merchant1982] which employs [Akaike1973].  The double Levinson recursion
   discussed by [Broersen2006] section 5.4 appears to be too numerically
   unstable to use in practice without requiring O(n^2) memory.

3. Implement the Ibrahim Optimum Tapered Burg as described by [Campbell1993]
   based on work in [Ibrahim1987a], [Ibrahim1987b], and [Ibrahim1989].  This
   should reduce the sensitivity to phase shifted input signals when working
   with small data sets.

References
----------

-- [Akaike1973]      Akaike, Hirotugu. "Block Toeplitz Matrix Inversion." SIAM Journal on Applied Mathematics 24 (March 1973): 234-241. http://dx.doi.org/10.1137/0124024

-- [Bernardo1976]    Bernardo, J. M.  "Algorithm AS 103: Psi (digamma) function." Journal of the Royal Statistical Society.  Series C (Applied Statistics) 25 (1976). http://www.jstor.org/stable/2347257

-- [Bourke1998]      Bourke, Paul. AutoRegression Analysis, November 1998. http://paulbourke.net/miscellaneous/ar/

-- [Box2008]         Box, George E. P., Gwilym M. Jenkins, and Gregory C. Reinsel. Time Series Analysis : Forecasting and Control. 4 edition. John Wiley, June 2008.

-- [Broersen2000]    Broersen, P. M. T. "Finite sample criteria for autoregressive order selection." IEEE Transactions on Signal Processing 48 (December 2000): 3550-3558. http://dx.doi.org/10.1109/78.887047

-- [Broersen2002]    Broersen, P. M. T. "Automatic spectral analysis with time series models." IEEE Transactions on Instrumentation and Measurement 51 (April 2002): 211-216. http://dx.doi.org/10.1109/19.997814

-- [Broersen2006]    Broersen, P. M. T. Automatic autocorrelation and spectral analysis. Springer, 2006. http://dx.doi.org/10.1007/1-84628-329-9

-- [Bunch1985]       Bunch, James R. "Stability of Methods for Solving Toeplitz Systems of Equations." SIAM Journal on Scientific and Statistical Computing 6 (1985): 349-364. http://dx.doi.org/10.1137/0906025

-- [Campbell1993]    Campbell, W. and D. N. Swingler. "Frequency estimation performance of several weighted Burg algorithms." IEEE Transactions on Signal Processing 41 (March 1993): 1237-1247. http://dx.doi.org/10.1109/78.205726

-- [Collomb2009]     Cedrick Collomb. Burg's method, algorithm, and recursion, November 2009. http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf

-- [GalassiGSL]      M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078. \url{http://www.gnu.org/software/gsl/}

-- [Hurvich1989]     Hurvich, Clifford M. and Chih-Ling Tsai. "Regression and time series model selection in small samples." Biometrika 76 (June 1989): 297-307. http://dx.doi.org/10.1093/biomet/76.2.297

-- [Ibrahim1987a]    Ibrahim, M. K. "Improvement in the speed of the data-adaptive weighted Burg technique." IEEE Transactions on Acoustics, Speech, and Signal Processing 35 (October 1987): 1474–1476. http://dx.doi.org/10.1109/TASSP.1987.1165046

-- [Ibrahim1987b]    Ibrahim, M. K. "On line splitting in the optimum tapered Burg algorithm." IEEE Transactions on Acoustics, Speech, and Signal Processing 35 (October 1987): 1476–1479. http://dx.doi.org/10.1109/TASSP.1987.1165047

-- [Ibrahim1989]     Ibrahim, M. K. "Correction to 'Improvement in the speed of the data-adaptive weighted Burg technique'." IEEE Transactions on Acoustics, Speech, and Signal Processing 37 (1989): 128. http://dx.doi.org/10.1109/29.17511

-- [Merchant1982]    Merchant, G. and T. Parks. "Efficient solution of a Toeplitz-plus-Hankel coefficient matrix system of equations." IEEE Transactions on Acoustics, Speech, and Signal Processing 30 (February 1982): 40-44. http://dx.doi.org/10.1109/TASSP.1982.1163845

-- [Press2007]       Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. Numerical recipes : The Art of Scientific Computing. Third edition. Cambridge University Press, September 2007.

-- [Seghouane2004]   Seghouane, A. K. and M. Bekara. "A Small Sample Model Selection Criterion Based on Kullback's Symmetric Divergence." IEEE Transactions on Signal Processing 52 (December 2004): 3314-3323. http://dx.doi.org/10.1109/TSP.2004.837416

-- [Trenberth1984]   Trenberth, K. E. "Some effects of finite sample size and persistence on meteorological statistics. Part I: Autocorrelations." Monthly Weather Review 112 (1984). http://dx.doi.org/10.1175/1520-0493(1984)112%3C2359:SEOFSS%3E2.0.CO;2

-- [Trench1967]      Trench, William F. Weighting coefficients for the prediction of stationary time series from the finite past. SIAM J. Appl. Math. 15, 6 (Nov. 1967), 1502-1510. http://www.jstor.org/stable/2099503

-- [Vandevender1982] Vandevender, W. H. and K. H. Haskell. "The SLATEC mathematical subroutine library." ACM SIGNUM Newsletter 17 (September 1982): 16-21.  http://dx.doi.org/10.1145/1057594.1057595

-- [Welford1962]     Welford, B. P. "Note on a Method for Calculating Corrected Sums of Squares and Products." Technometrics 4 (1962). http://www.jstor.org/stable/1266577

-- [Zohar1974]       Zohar, Shalhav. "The Solution of a Toeplitz Set of Linear Equations." J. ACM 21 (April 1974): 272-276. http://dx.doi.org/10.1145/321812.321822
