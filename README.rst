Autoregressive modeling tools in header-only C++
================================================

Overview
--------

This package contains a precision-agnostic, header-only, C++ implementation of
Burg's recursive method for estimating autoregressive model parameters.  Many
usability-related extensions, in particular Octave- and Python-friendly
functions, have been added to permit simply obtaining autocorrelation
information from the resulting estimated model.

The implementation permits extracting a sequence of AR(p) models for p from one
up to some maximum order::

	Estimating at most an AR(7) model using 10 samples
	
	AR      RMS/N      Gain Filter Coefficients
	--      -----      ---- -------------------
	 0   5.91e+00  1.00e+00 [ 1  ]
	 1   1.71e-03  3.45e+03 [ 1   -0.9999 ]
	 2   2.01e-05  2.94e+05 [ 1    -1.994   0.9941 ]
	 3   2.55e-08  2.32e+08 [ 1    -2.987    2.987  -0.9994 ]
	 4   1.01e-09  5.83e+09 [ 1    -3.967    5.913   -3.927   0.9799 ]
	 5   2.61e-11  2.26e+11 [ 1    -4.934    9.789   -9.763    4.895   -0.987 ]
	 6   3.42e-12  1.73e+12 [ 1    -5.854    14.35   -18.86    14.02   -5.586   0.9322 ]
	 7   3.65e-13  1.62e+13 [ 1    -6.735    19.63   -32.12    31.85   -19.15    6.465  -0.9452 ]
	
	AIC  selects model order 7 as best
	AICC selects model order 6 as best
	CIC  selects model order 7 as best

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
select the best model order given data from standard input.  It also estimates
the effective sample size and corresponding variance using ideas from
[Trenberth1984], [Thiebaux1984], and [vonStorch2001].  For example, ``arsel
--subtract-mean < rhoe.dat`` reproduces results from ARMASA [Broersen2002] on a
turbulence signal::

	# absrho    true
	# criterion CIC
	# eff_N     28.18777014115533
	# eff_var   3.6732508957963943e-05
	# gain      4249.4040527950729
	# maxorder  512
	# minorder  0
	# mu        0.20955287956200269
	# mu_sigma  0.0011415499935005066
	# N         1753
	# AR(p)     6
	# sigma2eps 8.3374920647988362e-09
	# sigma2x   3.5429372570302933e-05
	# submean   true
	# T0        62.190091348891279
	# window_T0 1
	+1
	-2.6990334396411866
	+2.8771681702855281
	-1.7247852051789097
	+0.75024605955486146
	-0.26866837869957461
	+0.06700587276734557


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

*Makefile*
   Try ``make`` followed by ``make check``.  On Linux, try ``make stress`` to
   examine the implementation's performance when piping in plain text data.
   Octave and/or Python functionality also will be built in-place when possible.

*ar.hpp*
  The standalone header implementing all algorithms.  Complete API
  documentation is available at http://rhysu.github.com/ar.

*arsel.cpp*
   Given data on standard input, use Burg's method to compute a hierarchy of
   candidate models and select the best one using CIC.  Try ``arsel --help`` to
   see available options.  This is perhaps the most useful standalone utility.

*arsel-octfile.cpp*, *arcov-octfile.cpp*
   Provides arsel.cpp-like capabilities for GNU Octave.  This is perhaps the
   most feature-rich way to start using these AR tools.  See appendix A
   ("Dynamically Linked Functions") within [Octave] for implementation details.
   Also demonstrates how working storage may be reused across multiple
   invocations to reduce the number of allocations for processing data sets.

*ar-python.cpp*, *setup.py*
   Provides some functionality as a Python extension module called 'ar'.
   This is modeled after the Octave wrapper but is not yet as robust.

*test.cpp*
   A test driver for testing ``ar.hpp`` against benchmarks by [Bourke1998].

*example.cpp*
   A test driver extracting a hierarchy of AR(p) models for a sample given by
   [Collomb2009].

*zohar.cpp*
   A test driver solving a nonsymmetric, real-valued Toeplitz set of linear
   equations.

*collomb2009.cpp*, *faber1986.cpp*
   For implementation testing and comparison purposes, a nearly verbatim copy
   of the recursive denominator algorithmic variant presented in
   [Kay1981,Faber1986] and [Collomb2009].  See comments at *issue3.dat*
   regarding numerical stability.

*test\*.coeff*, *test\*.dat*
   Sample data and exact parameters from [Bourke1998] used for ``make check``.

*rhoe.coeff*, *rhoe.dat*
   Sample turbulent total energy RMS fluctuation data and optimal parameters
   found by automatically by ARMASA [Broersen2002].

*issue3.dat*
   A large dataset from Nicholas Malaya generated by the Lorenz attractor.  For
   AR(4) and higher order models, this data tickles an instability present in
   [Andersen1978]'s recursive denominator variant of Burg's algorithm.  Namely,
   this variant will return a non-stationary process with complex poles outside
   the unit circle.  See https://github.com/RhysU/ar/issues/3 for details.

*WuleYalker.tex*
   A derivation of some equations closely connected with the Yule--Walker
   system.  Solving these permits recovering autocorrelations from process
   parameters.

*FiniteSampleCriteria.tex*
   A catalog of all implemented autoregressive model selection criteria.

*optionparser.h*
   The Lean Mean C++ Option Parser from http://optionparser.sourceforge.net
   which is used to parse command line arguments within sample applications.

References
----------

-- [Akaike1973]      Akaike, Hirotugu. "Block Toeplitz Matrix Inversion." SIAM Journal on Applied Mathematics 24 (March 1973): 234-241. http://dx.doi.org/10.1137/0124024

-- [Andersen1978]    Andersen, N. "Comments on the performance of maximum entropy algorithms." Proceedings of the IEEE 66 (November 1978): 1581-1582. http://dx.doi.org/10.1109/PROC.1978.11160

-- [Bernardo1976]    Bernardo, J. M.  "Algorithm AS 103: Psi (digamma) function." Journal of the Royal Statistical Society.  Series C (Applied Statistics) 25 (1976). http://www.jstor.org/stable/2347257

-- [Bourke1998]      Bourke, Paul. AutoRegression Analysis, November 1998. http://paulbourke.net/miscellaneous/ar/

-- [Box2008]         Box, George E. P., Gwilym M. Jenkins, and Gregory C. Reinsel. Time Series Analysis : Forecasting and Control. 4 edition. John Wiley, June 2008.

-- [Broersen2000]    Broersen, P. M. T. "Finite sample criteria for autoregressive order selection." IEEE Transactions on Signal Processing 48 (December 2000): 3550-3558. http://dx.doi.org/10.1109/78.887047

-- [Broersen2002]    Broersen, P. M. T. "Automatic spectral analysis with time series models." IEEE Transactions on Instrumentation and Measurement 51 (April 2002): 211-216. http://dx.doi.org/10.1109/19.997814

-- [Broersen2006]    Broersen, P. M. T. Automatic autocorrelation and spectral analysis. Springer, 2006. http://dx.doi.org/10.1007/1-84628-329-9

-- [Bunch1985]       Bunch, James R. "Stability of Methods for Solving Toeplitz Systems of Equations." SIAM Journal on Scientific and Statistical Computing 6 (1985): 349-364. http://dx.doi.org/10.1137/0906025

-- [Campbell1993]    Campbell, W. and D. N. Swingler. "Frequency estimation performance of several weighted Burg algorithms." IEEE Transactions on Signal Processing 41 (March 1993): 1237-1247. http://dx.doi.org/10.1109/78.205726

-- [Collomb2009]     Cedrick Collomb. "Burg's method, algorithm, and recursion", November 2009. http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf

-- [Faber1986]       Faber, L. J. "Commentary on the denominator recursion for Burg's block algorithm." Proceedings of the IEEE 74 (July 1986): 1046-1047. http://dx.doi.org/10.1109/PROC.1986.13584

-- [GalassiGSL]      M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078. \url{http://www.gnu.org/software/gsl/}

-- [Hurvich1989]     Hurvich, Clifford M. and Chih-Ling Tsai. "Regression and time series model selection in small samples." Biometrika 76 (June 1989): 297-307. http://dx.doi.org/10.1093/biomet/76.2.297

-- [Ibrahim1987a]    Ibrahim, M. K. "Improvement in the speed of the data-adaptive weighted Burg technique." IEEE Transactions on Acoustics, Speech, and Signal Processing 35 (October 1987): 1474–1476. http://dx.doi.org/10.1109/TASSP.1987.1165046

-- [Ibrahim1987b]    Ibrahim, M. K. "On line splitting in the optimum tapered Burg algorithm." IEEE Transactions on Acoustics, Speech, and Signal Processing 35 (October 1987): 1476–1479. http://dx.doi.org/10.1109/TASSP.1987.1165047

-- [Ibrahim1989]     Ibrahim, M. K. "Correction to 'Improvement in the speed of the data-adaptive weighted Burg technique'." IEEE Transactions on Acoustics, Speech, and Signal Processing 37 (1989): 128. http://dx.doi.org/10.1109/29.17511

-- [Kay1981]         Kay, S. M. and S. L. Marple. "Spectrum analysis- A modern perspective." Proceedings of the IEEE 69 (November 1981): 1380-1419. http://dx.doi.org/10.1109/PROC.1981.12184

-- [Merchant1982]    Merchant, G. and T. Parks. "Efficient solution of a Toeplitz-plus-Hankel coefficient matrix system of equations." IEEE Transactions on Acoustics, Speech, and Signal Processing 30 (February 1982): 40-44. http://dx.doi.org/10.1109/TASSP.1982.1163845

-- [Octave]          Eaton, John W., David Bateman, and Søren Hauberg. GNU Octave Manual Version 3. Network Theory Limited, 2008. http://www.octave.org/

-- [Press2007]       Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. Numerical recipes : The Art of Scientific Computing. Third edition. Cambridge University Press, September 2007.

-- [Seghouane2004]   Seghouane, A. K. and M. Bekara. "A Small Sample Model Selection Criterion Based on Kullback's Symmetric Divergence." IEEE Transactions on Signal Processing 52 (December 2004): 3314-3323. http://dx.doi.org/10.1109/TSP.2004.837416

-- [vonStorch2001]   Hans von Storch and Francis W. Zwiers. Statistical analysis in climate research. Cambridge University Press, March 2001. ISBN 978-0521012300.

-- [Thiebaux1984]    Thiébaux, H. J. and F. W. Zwiers. "The Interpretation and Estimation of Effective Sample Size." J. Climate Appl. Meteor. 23 (May 1984): 800-811. http://dx.doi.org/10.1175/1520-0450(1984)023%253C0800:TIAEOE%253E2.0.CO;2

-- [Trenberth1984]   Trenberth, K. E. "Some effects of finite sample size and persistence on meteorological statistics. Part I: Autocorrelations." Monthly Weather Review 112 (1984). http://dx.doi.org/10.1175/1520-0493(1984)112%3C2359:SEOFSS%3E2.0.CO;2

-- [Trench1967]      Trench, William F. Weighting coefficients for the prediction of stationary time series from the finite past. SIAM J. Appl. Math. 15, 6 (Nov. 1967), 1502-1510. http://www.jstor.org/stable/2099503

-- [Vandevender1982] Vandevender, W. H. and K. H. Haskell. "The SLATEC mathematical subroutine library." ACM SIGNUM Newsletter 17 (September 1982): 16-21.  http://dx.doi.org/10.1145/1057594.1057595

-- [Welford1962]     Welford, B. P. "Note on a Method for Calculating Corrected Sums of Squares and Products." Technometrics 4 (1962). http://www.jstor.org/stable/1266577

-- [Zohar1974]       Zohar, Shalhav. "The Solution of a Toeplitz Set of Linear Equations." J. ACM 21 (April 1974): 272-276. http://dx.doi.org/10.1145/321812.321822
