// Except for any way in which it interferes with Cedrick Collomb's 2009
// copyright assertion in the article "Burg’s Method, Algorithm and Recursion":
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BURG_HPP
#define BURG_HPP

// The selection criteria routines might be sped up for floating point
// arguments given an appropriate digamma (psi) or Pochhammer symbol
// implementation.  To do so with the GNU Scientific Library, e.g., try
//
//     #include <gsl/gsl_sf_psi.h>
//     #define BURG_DIGAMMA(x) gsl_sf_psi(x)
//
//     #include <gsl/gsl_sf_gamma.h>
//     #define BURG_POCHHAMMER(a,x) gsl_sf_poch(a,x)
//
// before including this header.  Be sure to link the GSL with your binary.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <vector>

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Algorithms for autoregressive parameter estimation and manipulation.
 *
 * @{
 */

/**
 * Fit an autoregressive model to stationary time series data using
 * Burg's method.  That is, assuming a zero-mean model
 * \f{align}{
 *     x_n + a_1 x_{n - 1} + \dots + a_p x_{n - p} &= \epsilon_n
 *     &
 *     \epsilon_n &\sim{} N\left(0, \sigma^2_epsilon\right)
 *     \\
 *     \rho_k + a_1 \rho_{k-1} + \dots + a_p \rho_{k-p} &= 0
 *     &
 *     k &\geq{} p
 * \f}
 * find coefficients \f$a_i\f$ such that the sum of the squared errors in the
 * forward predictions \f$x_n = -a_1 x_{n-1} - \dots - a_p x_{n-p}\f$ and
 * backward predictions \f$x_n = -a_1 x_{n+1} - \dots - a_p x_{n+p}\f$ are both
 * minimized.  Either a single model of given order or a hierarchy of models up
 * to and including a maximum order may fit.
 *
 * The input data \f$\vec{x}\f$ are read from <tt>[data_first,data_last)</tt>
 * in a single pass.  The mean is computed using pairwise summation,
 * returned in \c mean, and \e removed from further consideration whenever
 * \c subtract_mean is true.  The estimated model parameters \f$a_i\f$
 * are output using \c params_first with the behavior determined by
 * the amount of data read, <tt>maxorder</tt>, and the \c hierarchy flag:
 * <ul>
 *     <li>If \c hierarchy is \c false, only the \f$a_1, \dots,
 *         a_\text{maxorder}\f$ parameters for an AR(<tt>maxorder</tt>) process
 *         are output.</li>
 *     <li>If \c hierarchy is true, the <tt>maxorder*(maxorder+1)/2</tt>
 *         coefficients for models AR(1), AR(2), ..., AR(maxorder) are output.
 *         </li>
 * </ul>
 * Note that the latter case is \e always computed;  The \c hierarchy flag
 * merely controls what is output.  In both cases, the maximum order is limited
 * by the number of data samples provided and is output to \c maxorder.
 *
 * One mean squared discrepancy \f$\sigma^2_\epsilon\f$, also called the
 * innovation variance, and gain \f$\sigma^2_x / \sigma^2_\epsilon\f$ are
 * output for each model using \c sigma2e_first and \c gain_first.  The
 * autocorrelations for lags <tt>[1,k]</tt> are output using \c autocor_first.
 * When \c hierarchy is true, only lags <tt>[1,m]</tt> should be applied for
 * some AR(<tt>m</tt>) model.  The lag zero autocorrelation is always equal to
 * one and is therefore never output.  Outputting the lag \c k autocorrelation
 * is technically redundant as it may be computed from \f$a_i\f$ and lags
 * <tt>1, ..., k-1</tt>.  Autocovariances may be computed by multiplying the
 * autocorrelations by \f$\text{gain} \sigma^2_\epsilon\f$.
 *
 * The implementation has been refactored
 * heavily from Cedrick Collomb's 2009 article <a
 * href="http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf">"Burg’s
 * Method, Algorithm and Recursion"</a>.  In particular, iterators
 * are employed, the working precision is selectable using \c mean,
 * the mean squared discrepancy calculation has been added, some loop
 * index transformations have been performed, and all lower order models
 * may be output during the recursion using \c hierarchy.  Gain and
 * autocorrelation calculations have been added based on sections 5.2
 * and 5.3 of Broersen, P.  M.  T. Automatic autocorrelation and spectral
 * analysis. Springer, 2006.  http://dx.doi.org/10.1007/1-84628-329-9.
 *
 * @param[in]     data_first    Beginning of the input data range.
 * @param[in]     data_last     Exclusive end of the input data range.
 * @param[out]    mean          Mean of data computed using pairwise summation.
 * @param[in,out] maxorder      On input, the maximum model order desired.
 *                              On output, the maximum model order computed.
 * @param[out]    params_first  Model parameters for a single model or
 *                              for an entire hierarchy of models
 *                              depending upon \c hierarchy.
 * @param[out]    sigma2e_first Mean squared discrepancy for a single model or
 *                              for an entire hierarchy of models
 *                              depending upon \c hierarchy.
 * @param[out]    gain_first    Model gain for a single model
 *                              for an entire hierarchy of models
 *                              depending upon \c hierarchy.
 * @param[out]    autocor_first Lag one through lag maxorder autocorrelations.
 * @param[in]     subtract_mean Should \c mean be subtracted from the data?
 * @param[in]     hierarchy     Should the entire hierarchy of model fits
 *                              be output?
 *
 * @returns the number data values processed within [data_first, data_last).
 */
template <class InputIterator,
          class Value,
          class OutputIterator1,
          class OutputIterator2,
          class OutputIterator3,
          class OutputIterator4>
std::size_t burg_method(InputIterator     data_first,
                        InputIterator     data_last,
                        Value&            mean,
                        std::size_t&      maxorder,
                        OutputIterator1   params_first,
                        OutputIterator2   sigma2e_first,
                        OutputIterator3   gain_first,
                        OutputIterator4   autocor_first,
                        const bool        subtract_mean = false,
                        const bool        hierarchy = false)
{
    using std::bind2nd;
    using std::copy;
    using std::distance;
    using std::fill;
    using std::inner_product;
    using std::min;
    using std::minus;

    typedef typename std::vector<Value> vector;
    typedef typename vector::size_type size;

    // Initialize f from [data_first, data_last) and fix number of samples
    vector f(data_first, data_last);
    const size N = f.size();

    // Compute the mean of f using pairwise summation and output it.
    // Pairwise chosen instead of of Kahan for speed trade off and to avoid
    // algorithmic nonsense when working precision is exact (e.g. rational).
    vector b(N, 0);
    {
        // First pass copies f into b reducing by up to a factor of 2
        // Requires b has been adequately sized and filled with zeros
        for (size i = 0; i < N; ++i)
            b[i/2] += f[i];

        // Initialize i to the next power of 2 lower than N
        size t = N, i = !!N;
        while (t /= 2)
            i *= 2;

        // Recurse on the now power-of-two, smaller problem size
        while (i /= 2)
            for (size_t j = 0; j < i; ++j)
                b[j] = b[2*j] + b[2*j+1];
    }
    mean = b[0] / N;

    // At most maxorder N-1 can be fit from N samples.  Beware N may be zero.
    maxorder = min(static_cast<size>(maxorder)+1, N)-1;

    // Short circuit if no work was requested or is possible.
    if (maxorder == 0) return N;

    // Subtract the mean of the data if requested
    if (subtract_mean)
        transform(f.begin(), f.end(), f.begin(), bind2nd(minus<Value>(), mean));

    // Initialize mean squared discrepancy sigma2e and Dk
    Value sigma2e = inner_product(f.begin(), f.end(), f.begin(), Value(0));
    Value Dk = - f[0]*f[0] - f[N - 1]*f[N - 1] + 2*sigma2e;
    sigma2e /= N;

    // Initialize recursion
    copy(f.begin(), f.end(), b.begin());
    vector Ak(maxorder + 1, Value(0));
    Ak[0] = 1;
    Value gain = 1;
    vector autocor;
    autocor.reserve(maxorder);

    // Perform Burg recursion
    for (size kp1 = 1; kp1 <= maxorder; ++kp1)
    {
        // Compute mu from f, b, and Dk and then update sigma2e and Ak using mu
        // Afterwards, Ak[1:kp1] contains AR(k) coefficients by the recurrence
        // By the recurrence, Ak[kp1] will also be the reflection coefficient
        const Value mu = 2/Dk*inner_product(f.begin() + kp1, f.end(),
                                            b.begin(), Value(0));
        sigma2e *= (1 - mu*mu);
        for (size n = 0; n <= kp1/2; ++n)
        {
            Value t1 = Ak[n] - mu*Ak[kp1 - n];
            Value t2 = Ak[kp1 - n] - mu*Ak[n];
            Ak[n] = t1;
            Ak[kp1 - n] = t2;
        }

        // Update the gain per Broersen 2006 equation (5.25)
        gain *= 1 / (1 - Ak[kp1]*Ak[kp1]);

        // Compute and output the next autocorrelation coefficient
        // See Broersen 2006 equations (5.28) and (5.31) for details
        autocor.push_back(-inner_product(autocor.rbegin(), autocor.rend(),
                                         Ak.begin() + 1, Ak[kp1]));

        // Output parameters and the input and output variances when requested
        if (hierarchy || kp1 == maxorder)
        {
            params_first = copy(Ak.begin() + 1, Ak.begin() + kp1 + 1,
                                params_first);
            *sigma2e_first++ = sigma2e;
            *gain_first++    = gain;
        }

        // Update f, b, and then Dk for the next iteration if another remains
        if (kp1 < maxorder)
        {
            for (size n = 0; n < N - kp1; ++n)
            {
                Value t1 = f[n + kp1] - mu*b[n];
                Value t2 = b[n] - mu*f[n + kp1];
                f[n + kp1] = t1;
                b[n] = t2;
            }
            Dk = (1 - mu*mu)*Dk - f[kp1]*f[kp1] - b[N - kp1 - 1]*b[N - kp1 - 1];
        }
    }

    // Output the lag [1,maxorder] autocorrelation coefficients in single pass
    copy(autocor.begin(), autocor.end(), autocor_first);

    // Return the number of values processed in [data_first, data_last)
    return N;
}

/**
 * Solve a Toeplitz set of linear equations.  That is, find \f$s_{n+1}\f$
 * satisfying
 * \f[
 *      L_{n+1} s_{n+1} = d_{n+1}
 *      \mbox{ where }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    r_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{a}\f$, \f$\vec{r}\f$, and \f$\vec{d}\f$.  The dimension of the
 * problem is fixed by <tt>n = distance(a_first, a_last)</tt>.  A symmetric
 * Toeplitz solve can be performed by having \f$\vec{a}\f$ and \f$vec{r}\f$
 * iterate over the same data.  The Hermitian case requires two buffers with
 * \f$vec{r}\f$ being the conjugate of \f$\vec{a}\f$.  The working precision
 * is fixed by the \c value_type of \c d_first.
 *
 * The algorithm is from Zohar, Shalhav. "The Solution of a Toeplitz Set of
 * Linear Equations." J. ACM 21 (April 1974): 272-276.
 * http://dx.doi.org/10.1145/321812.321822.  It has complexity like
 * <tt>O(2*(n+1)^2)</tt>.  Zohar improved upon earlier work from Page 1504 from
 * Trench, William F. "Weighting Coefficients for the Prediction of Stationary
 * Time Series from the Finite Past." SIAM Journal on Applied Mathematics 15
 * (November 1967): 1502-1510.  http://www.jstor.org/stable/2099503.  See
 * Bunch, James R. "Stability of Methods for Solving Toeplitz Systems of
 * Equations." SIAM Journal on Scientific and Statistical Computing 6 (1985):
 * 349-364. http://dx.doi.org/10.1137/0906025 for a discussion of the
 * algorithms stability characteristics.
 *
 * @param[in]  a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]  a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in]  r_first Beginning of the range containing \f$\vec{r}\f$.
 * @param[in]  d_first Beginning of the range containing \f$\vec{d}\f$
 *                     which should have <tt>n+1</tt> entries available.
 * @param[out] s_first Beginning of the output range to which
 *                     <tt>n+1</tt> entries will be written.
 */
template<class RandomAccessIterator,
         class InputIterator,
         class OutputIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        RandomAccessIterator r_first,
                        InputIterator        d_first,
                        OutputIterator       s_first)
{
    using std::copy;
    using std::distance;
    using std::inner_product;
    using std::reverse_iterator;

    // Tildes indicate transposes while hats indicate reversed vectors.

    // InputIterator::value_type determines the working precision
    typedef typename std::iterator_traits<InputIterator>::value_type value;
    typedef typename std::vector<value> vector;
    typedef typename vector::size_type size;

    // Determine problem size using [a_first,a_last) and ensure nontrivial
    typename std::iterator_traits<RandomAccessIterator>::difference_type dist
            = distance(a_first, a_last);
    if (dist < 1) throw std::invalid_argument("distance(a_first, a_last) < 1");
    const size n = static_cast<size>(dist);

    // Allocate working storage and set initial values for recursion:
    vector s;    s   .reserve(n+1); s   .push_back( *d_first);
    vector ehat; ehat.reserve(n  ); ehat.push_back(-a_first[0]);
    vector g;    g   .reserve(n  ); g   .push_back(-r_first[0]);
    value lambda  = 1 - a_first[0]*r_first[0];

    // Though recursive updates to s and g can be done in-place, updates to
    // ehat seemingly require one additional vector for storage:
    //
    // "This sequence of computations is economical of storage.  It is only
    // necessary to retain quantities computed at level m - 1 until the
    // computations at level m are complete." [Trench1967, page 1504]
    vector next_ehat; next_ehat.reserve(n);

    // Recursion for i = {1, 2, ..., n - 1}:
    for (size i = 1; i < n; ++i)
    {

        reverse_iterator<RandomAccessIterator> rhat_first(r_first + i);

        // \theta_i =  \delta_{i+1}  - \tilde{s}_i \hat{r}_i
        const value neg_theta = inner_product(
                s.begin(), s.end(), rhat_first, value(-(*++d_first)));

        // \eta_i   = -\rho_{-(i+1)} - \tilde{a}_i \hat{e}_i
        const value neg_eta   = inner_product(
                ehat.begin(), ehat.end(), a_first, value(a_first[i]));

        // \gamma_i = -\rho_{i+1}    - \tilde{g}_i \hat{r}_i
        const value neg_gamma = inner_product(
                g.begin(), g.end(), rhat_first, value(r_first[i]));

        /*
         * s_{i+1} = \bigl(\begin{smallmatrix}
         *              s_i + (\theta_i/\lambda_i) \hat{e}_i \\
         *              \theta_i/\lambda_i
         *          \end{smallmatrix}\bigr)
         *
         * \hat{e}_{i+1} = \bigl(\begin{smallmatrix}
         *                     \eta_i/\lambda_i \\
         *                     \hat{e}_i + (\eta_i/\lambda_i) g_i
         *                 \end{smallmatrix}\bigr)
         *
         * g_{i+1} = \bigl(\begin{smallmatrix}
         *               g_i + (\gamma_i/\lambda_i) \hat{e}_i \\
         *               \gamma_i/\lambda_i
         *           \end{smallmatrix}\bigr)
         */
        const value theta_by_lambda = -neg_theta/lambda;
        const value   eta_by_lambda = -neg_eta  /lambda;
        const value gamma_by_lambda = -neg_gamma/lambda;
        next_ehat.clear();
        next_ehat.push_back(eta_by_lambda);
        for (size j = 0; j < i; ++j)
        {
            s[j] += theta_by_lambda*ehat[j];
            next_ehat.push_back(ehat[j] + eta_by_lambda*g[j]);
            g[j] += gamma_by_lambda*ehat[j];
        }
        s.push_back(theta_by_lambda);
        g.push_back(gamma_by_lambda);
        ehat.swap(next_ehat);

        // \lambda_{i+1} = \lambda_i - \eta_i \gamma_i / \lambda_i
        lambda -= neg_eta*neg_gamma/lambda;
    }

    // Recursion for i = n differs slightly per Zohar's "Last Computed Values"
    // Computing g_n above was unnecessary but the incremental expense is small
    {
        reverse_iterator<RandomAccessIterator> rhat_first(r_first + n);

        // \theta_n =  \delta_{n+1}  - \tilde{s}_n \hat{r}_n
        const value neg_theta = inner_product(
                s.begin(), s.end(), rhat_first, value(-(*++d_first)));

        /*
         * s_{n+1} = \bigl(\begin{smallmatrix}
         *              s_n + (\theta_n/\lambda_n) \hat{e}_n \\
         *              \theta_n/\lambda_n
         *          \end{smallmatrix}\bigr)
         */
        const value theta_by_lambda = -neg_theta/lambda;
        for (size j = 0; j < n; ++j)
        {
            s[j] += theta_by_lambda*ehat[j];
        }
        s.push_back(theta_by_lambda);
    }

    // Output solution
    copy(s.begin(), s.end(), s_first);
}

/**
 * Solve a Toeplitz set of linear equations in-place.  That is, compute
 * \f[
 *      L_{n+1}^{-1} d_{n+1}
 *      \mbox{ for }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    r_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{a}\f$, \f$\vec{r}\f$, and \f$\vec{d}\f$.  The dimension of the
 * problem is fixed by <tt>n = distance(a_first, a_last)</tt>.  A symmetric
 * Toeplitz solve can be performed by having \f$\vec{a}\f$ and \f$vec{r}\f$
 * iterate over the same data.  The Hermitian case requires two buffers with
 * \f$vec{r}\f$ being the conjugate of \f$\vec{a}\f$.  The working precision
 * is fixed by the \c value_type of \c d_first.
 *
 * @param[in]     a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]     a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in]     r_first Beginning of the range containing \f$\vec{r}\f$.
 * @param[in,out] d_first Beginning of the range containing \f$\vec{d}\f$.
 *                        Also the beginning of the output range to which
 *                        <strong><tt>n+1</tt></strong> entries will be
 *                        written.
 */
template<class RandomAccessIterator,
         class ForwardIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        RandomAccessIterator r_first,
                        ForwardIterator      d_first)
{
    return zohar_linear_solve(a_first, a_last, r_first, d_first, d_first);
}


/**
 * Solve a real-valued, symmetric Toeplitz set of linear equations in-place.
 * That is, compute
 * \f[
 *      L_{n+1}^{-1} d_{n+1}
 *      \mbox{ for }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    a_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{a}\f$ and \f$\vec{d}\f$.  The dimension of the problem is
 * fixed by <tt>n = distance(a_first, a_last)</tt>.  The working precision is
 * fixed by the \c value_type of \c d_first.
 *
 * @param[in]     a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]     a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in,out] d_first Beginning of the range containing \f$\vec{d}\f$.
 *                        Also the beginning of the output range to which
 *                        <strong><tt>n+1</tt></strong> entries will be
 *                        written.
 */
template<class RandomAccessIterator,
         class ForwardIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        ForwardIterator      d_first)
{
    return zohar_linear_solve(a_first, a_last, a_first, d_first);
}

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Method-specific estimation variance routines following Broersen.
 *
 * For details see either the FiniteSampleCriteria.tex write up or Broersen,
 * P. M. T. "Finite sample criteria for autoregressive order selection." IEEE
 * Transactions on Signal Processing 48 (December 2000): 3550-3558.
 * http://dx.doi.org/10.1109/78.887047.
 *
 * @{
 */

/** Denotes the sample mean was subtracted from a signal before estimation. */
struct mean_subtracted
{
    /** Computes the empirical variance estimate for order zero. */
    template <typename Result, typename Integer>
    static Result empirical_variance_zero(Integer N)
    {
        assert(N >= 1);

        Result den = N;
        return 1 / den;
    }
};

/** Denotes the sample mean was retained in a signal during estimation. */
struct mean_retained
{
    /** Computes the empirical variance estimate for order zero. */
    template <typename Result, typename Integer>
    static Result empirical_variance_zero(Integer)
    {
        return 0;
    }
};

/**
 * A parent type for autoregressive process parameter estimation techniques.
 *
 * Each subclass should have an <tt>empirical_variance(N, i)</tt> method
 * following Broersen, P. M. and H. E. Wensink. "On Finite Sample Theory for
 * Autoregressive Model Order Selection." IEEE Transactions on Signal
 * Processing 41 (January 1993): 194+.
 * http://dx.doi.org/10.1109/TSP.1993.193138.
 */
struct estimation_method {};

/** Represents estimation by solving the Yule--Walker equations. */
template <class MeanHandling>
class YuleWalker : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(i >= 0);
        assert(i <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        Result num = N - i, den = N*(N + 2);
        return num / den;
    }
};

/** Represents estimation using Burg's recursive method. */
template <class MeanHandling>
class Burg : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(i >= 0);
        assert(i <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        Result den = N + 1 - i;
        return 1 / den;
    }
};

/** Represents forward and backward prediction least squares minimization. */
template <class MeanHandling>
class LSFB : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(i >= 0);
        assert(i <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        // Factorizing expression will cause problems in unsigned arithmetic
        Result den = N + Result(3)/2 - Result(3)/2 * i;
        return 1 / den;
    }
};

/** Represents forward prediction least squares minimization. */
template <class MeanHandling>
class LSF : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(i >= 0);
        assert(i <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        // Factorizing expression will cause problems in unsigned arithmetic
        Result den = N + 2 - 2*i;
        return 1 / den;
    }
};

/** An STL-ready binary_function for a given method's empirical variance. */
template <class EstimationMethod,
          typename Result,
          typename Integer1 = std::size_t,
          typename Integer2 = Integer1>
struct empirical_variance_function
    : public std::binary_function<Result,Integer1,Integer2>
{
    Result operator() (Integer1 N, Integer2 i) const
    {
        return EstimationMethod::template empirical_variance<Result>(N, i);
    }
};

/**
 * An STL AdaptableGenerator for a given method's empirical variance.
 *
 * On each <tt>operator()</tt> invocation the method's empirical variance is
 * returned for the current model order.  The first invocation returns the
 * result for model order zero.
 */
template <class EstimationMethod,
          typename Result,
          typename Integer1 = std::size_t,
          typename Integer2 = Integer1>
class empirical_variance_generator
{
private:

    Integer1 N;
    Integer2 i;

public:

    typedef Result result_type;

    empirical_variance_generator(Integer1 N) : N(N), i(0) {}

    Result operator() ()
    {
        return EstimationMethod::template empirical_variance<Result>(N, i++);
    }
};

/**
 * An immutable RandomAccessIterator over a method's empirical variance
 * sequence.
 *
 * Facilitates using algorithms like <tt>std::copy</tt>,
 * <tt>std::accumulate</tt>, and <tt>std::partial_sum</tt> when comparing a
 * hierarchy of models during model order selection.
 *
 * The (N+1)-length sequence of orders 0, 1, ..., N is iterated given sample
 * size N.  Default constructed instances represent past-end iterators.
 */
template <class EstimationMethod,
          typename Result,
          typename Integer1 = std::size_t,
          typename Integer2 = Integer1>
class empirical_variance_iterator
    : public std::iterator<std::random_access_iterator_tag, Result>
{
private:
    typedef std::iterator<std::random_access_iterator_tag, Result> base;

    empirical_variance_iterator(Integer1 N, Integer2 i) : N(N), i(i) {}

    Integer1 N;
    Integer2 i;

public:
    typedef typename base::difference_type   difference_type;
    typedef typename base::iterator_category iterator_category;
    typedef typename base::pointer           pointer;
    typedef typename base::reference         reference;
    typedef typename base::value_type        value_type;

    /** Construct a past-end iterator */
    empirical_variance_iterator() : N(0), i(1) {}  // Null

    /** Construct an iterator over sequence order 0, 1, ..., N (inclusive). */
    empirical_variance_iterator(Integer1 N) : N(N), i(0)
        { assert(N >= 1); }

    // Forward traversal

    empirical_variance_iterator& operator++()
        { assert(N >= 1); ++i; return *this; }

    empirical_variance_iterator operator++(int)
        { assert(N >= 1); return empirical_variance_iterator(N, i++); }

    empirical_variance_iterator operator+(const difference_type& k)
        { assert(N >= 1); return empirical_variance_iterator(N, i + k); }

    empirical_variance_iterator& operator+=(const difference_type& k)
        { assert(N >= 1); i += k; return *this; }

    // Backward traversal

    empirical_variance_iterator& operator--()
        { assert(N >= 1); --i; return *this; }

    empirical_variance_iterator operator--(int)
        { assert(N >= 1); return empirical_variance_iterator(N, i--); }

    empirical_variance_iterator operator-(const difference_type& k)
        { assert(N >= 1); return empirical_variance_iterator(N, i - k); }

    empirical_variance_iterator& operator-=(const difference_type& k)
        { assert(N >= 1); i -= k; return *this; }

    // Distance support

    difference_type operator-(const empirical_variance_iterator& other)
    {
        if (!this->N) {
            return 1 + other.N - other.i;
        } else if (!other.N) {
            return -static_cast<difference_type>(1 + this->N - this->i);
        } else {
            assert(this->N == other.N);
            return this->i - other.i;
        }
    }

    // EqualityComparable

    bool operator==(const empirical_variance_iterator& other)
    {
        if (!this->N) {
            return other.i >= static_cast<Integer2>(other.N + 1);
        } else if (!other.N) {
            return this->i >= static_cast<Integer2>(this->N + 1);
        } else {
            return this->N == other.N && this->i == other.i;
        }
    }

    bool operator!=(const empirical_variance_iterator& other)
        { return !(*this == other); }

    // LessThanComparable will trigger assertion on nonsense N cases

    bool operator<(const empirical_variance_iterator& other)
        { return (*this - other) < 0; }

    bool operator<=(const empirical_variance_iterator& other)
        { return (*this - other) <= 0; }

    bool operator>(const empirical_variance_iterator& other)
        { return (*this - other) > 0; }

    bool operator>=(const empirical_variance_iterator& other)
        { return (*this - other) >= 0; }

    // Dereference operations

    value_type operator*() const
    {
        assert(i <= N);

        return EstimationMethod::template empirical_variance<Result>(N, i);
    }

    value_type operator[](const difference_type &k) const
    {
        assert(i + k <= N);

        return EstimationMethod::template empirical_variance<Result>(N, i + k);
    }
};

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Criteria for autoregressive model order selection following Broersen.
 *
 * For details see either the FiniteSampleCriteria.tex write up or Broersen,
 * P. M. T. "Finite sample criteria for autoregressive order selection." IEEE
 * Transactions on Signal Processing 48 (December 2000): 3550-3558.
 * http://dx.doi.org/10.1109/78.887047.
 *
 * @{
 */

/**
 * A parent type for autoregressive model selection criterion.
 *
 * Each subclass should have an <tt>overfit_penalty(N, p)</tt> method
 * following Broersen, P. M. and H. E. Wensink. "On Finite Sample Theory for
 * Autoregressive Model Order Selection." IEEE Transactions on Signal
 * Processing 41 (January 1993): 194+.
 * http://dx.doi.org/10.1109/TSP.1993.193138.
 */
struct criterion
{
    /** Compute the underfit penalty given \f$\sigma^2_\epsilon. */
    template <typename Result, typename Input>
    static Result underfit_penalty(Input sigma2e)
    {
        using std::log;
        return log(Result(sigma2e));
    }
};

/**
 * Represents the generalized information criterion (GIC).  The penalty factor
 * \f$\alpha\f$ is controlled by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <int AlphaNumerator = 3, int AlphaDenominator = 1>
struct GIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        return Result(AlphaNumerator) * p / (N * AlphaDenominator);
    }
};

/** Represents the Akaike information criterion (AIC). */
struct AIC : public GIC<2> {};

/** Represents the consistent criterion BIC. */
struct BIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        using std::log;
        return log(Result(N)) * p / N;
    }
};

/** Represents the minimally consistent criterion (MCC). */
struct MCC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        using std::log;
        return 2*log(log(Result(N))) * p / N;
    }
};

/**
 * Represents the asymptotically-corrected Akaike information criterion (AICC).
 */
struct AICC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        return 2 * Result(p) / (N - p - 1);
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to a particular
 * \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled by
 * <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class EstimationMethod,
          int AlphaNumerator = 3,
          int AlphaDenominator = 1 >
struct FIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        // This accumulate invocation is inefficient but always correct
        using std::accumulate;
        typedef empirical_variance_iterator<
                EstimationMethod, Result, Integer1, Integer2
            > evi_type;
        Result sum = accumulate(evi_type(N), evi_type(N) + p, Result(0));

        return AlphaNumerator * sum / AlphaDenominator;
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * YuleWalker \ref estimation_method.  The penalty factor \f$\alpha\f$ is
 * controlled by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<YuleWalker<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = YuleWalker<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));

        Result num = (2*N - p - 1)*p;  // Avoids non-positive values
        Result den = 2*N*(N + 2);

        return AlphaNumerator * (v0 + num/den) / AlphaDenominator;
    }
};

// Specializations of the FIC for efficiency when digamma is available.
#ifdef BURG_DIGAMMA

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * Burg \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled
 * by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<Burg<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result t  = Burg<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        t -= BURG_DIGAMMA(N + 1);
        t += BURG_DIGAMMA(N + 1 - p);
        return AlphaNumerator * t / AlphaDenominator;
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * LSFB \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled
 * by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<LSFB<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result t  = LSFB<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = BURG_DIGAMMA((3 + 2*N) / Result(3)    );
        Result b = BURG_DIGAMMA((3 + 2*N) / Result(3) - p);
        return AlphaNumerator * (t - 2*(a - b)/3) / AlphaDenominator;
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * LSF \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled
 * by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<LSF<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result t  = LSF<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = BURG_DIGAMMA((2 + N) / Result(2)    );
        Result b = BURG_DIGAMMA((2 + N) / Result(2) - p);
        return AlphaNumerator * (t - (a - b)/2) / AlphaDenominator;
    }
};

#endif /* BURG_DIGAMMA */

/**
 * Represents the finite sample information criterion (FSIC) as applied to a
 * particular \ref estimation_method.
 */
template <class EstimationMethod>
struct FSIC : public criterion
{

private:

    /** A helper to compute \f$\frac{1+v}{1-v}\f$ for the method in use. */
    template <typename Result, typename Integer1, typename Integer2>
    class product_iterator
        : public empirical_variance_iterator<
                EstimationMethod, Result, Integer1, Integer2
          >
    {
    private:
        typedef empirical_variance_iterator<
                EstimationMethod, Result, Integer1, Integer2
            > base;

    public:
        typedef typename base::difference_type   difference_type;
        typedef typename base::iterator_category iterator_category;
        typedef typename base::pointer           pointer;
        typedef typename base::reference         reference;
        typedef typename base::value_type        value_type;

        product_iterator(Integer1 N) : base(N) {}

        value_type operator*() const
        {
            const value_type v = this->base::operator*();
            return (1 + v) / (1 - v);
        }

        value_type operator[](const difference_type &k) const
        {
            const value_type v = this->base::operator[](k);
            return (1 + v) / (1 - v);
        }

    };

public:

    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        // This accumulate invocation is inefficient but always correct
        using std::multiplies;
        using std::accumulate;
        product_iterator<Result, Integer1, Integer2> first(N), last(N);
        last += p; // Avoids overloading requirements for (last + p)

        return accumulate(first, last, Result(1), multiplies<Result>()) - 1;
    }
};

// Specializations of the FIC for efficiency when Pochhammer is available.
#ifdef BURG_POCHHAMMER

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * YuleWalker \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<YuleWalker<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = YuleWalker<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = BURG_POCHHAMMER(N*(N+3) - p, p);
        Result b = BURG_POCHHAMMER(Result(1 + N) - N*N);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * Burg \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<Burg<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = Burg<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = BURG_POCHHAMMER(-Result(1) - N, p);
        Result b = BURG_POCHHAMMER( Result(1) - N, p);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * LSFB \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<LSFB<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = LSFB<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = BURG_POCHHAMMER((-Result(2)-2*N)/3, p);
        Result b = BURG_POCHHAMMER(( Result(2)-2*N)/3, p);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * LSF \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<LSF<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = LSF<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = BURG_POCHHAMMER((-Result(1)-N)/2, p);
        Result b = BURG_POCHHAMMER(( Result(1)-N)/2, p);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

#endif /* BURG_POCHHAMMER */


/**
 * Represents the combined information criterion (CIC) as applied to a
 * particular \ref estimation_method.
 */
template <class EstimationMethod>
struct CIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        using std::max;
        return max(
                FSIC<EstimationMethod>::template overfit_penalty<Result>(N, p),
                FIC <EstimationMethod>::template overfit_penalty<Result>(N, p)
            );
    }
};

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Algorithmic helpers for autoregressive model order selection.
 *
 * @{
 */

/**
 * Find the best model from a hierarchy of candidates using a \ref criterion.
 *
 * @param[in]  N        Sample count used to compute \f$\sigma^2_\epsilon\f$.
 * @param[in]  ordfirst The model order corresponding to \c first.
 *                      When \f$sigma^2_epsilon\f$ is produced entirely by
 *                      \ref burg_method, this should be \c 1u.
 * @param[in]  first    Beginning of the range holding \f$\sigma^2_\epsilon\f$
 * @param[in]  last     Exclusive end of input range.
 * @param[out] crit     Value assigned to each model by the criterion.
 *
 * @return The distance from \c first to the best model.
 */
template <class    Criterion,
          typename Integer1,
          typename Integer2,
          class    InputIterator,
          class    OutputIterator>
typename std::iterator_traits<InputIterator>::difference_type
select_model(Integer1       N,
             Integer2       ordfirst,
             InputIterator  first,
             InputIterator  last,
             OutputIterator crit)
{
    using std::iterator_traits;
    using std::numeric_limits;

    typedef InputIterator iterator;
    typedef typename iterator_traits<iterator>::difference_type difference;
    typedef typename iterator_traits<iterator>::value_type      value;

    value best_val = numeric_limits<value>::max();
    difference best_pos = -1, dist = -1;

    while (first != last) {

        value underfit  = Criterion::template underfit_penalty<value>(
                                *first++);
        value overfit   = Criterion::template overfit_penalty<value>(
                                N, ++dist + ordfirst);
        value candidate = underfit + overfit;
        *crit++ = candidate;

        if (candidate < best_val) {
            best_val = candidate;
            best_pos = dist;
        }

    }

    return best_pos;;
}

namespace { // anonymous

// Used to discard output per http://stackoverflow.com/questions/335930/
struct null_output : std::iterator< std::output_iterator_tag, null_output >
{
    template <typename T> void operator=(const T&) {}

    null_output& operator++()
        { return *this; }

    null_output operator++(int)
        { null_output it(*this); ++*this; return it; }

    null_output& operator*()
        { return *this; }
};

}

/**
 * Find the best model from a hierarchy of candidates using a \ref criterion.
 *
 * @param[in]  N        Sample count used to compute \f$\sigma^2_\epsilon\f$.
 * @param[in]  ordfirst The model order corresponding to \c first.
 *                      When \f$sigma^2_epsilon\f$ is produced entirely by
 *                      \ref burg_method, this should be \c 1u.
 * @param[in]  first    Beginning of the range holding \f$\sigma^2_\epsilon\f$
 * @param[in]  last     Exclusive end of input range.
 *
 * @return The distance from \c first to the best model.
 */
template <class    Criterion,
          typename Integer1,
          typename Integer2,
          class    InputIterator>
typename std::iterator_traits<InputIterator>::difference_type
select_model(Integer1      N,
             Integer2      ordfirst,
             InputIterator first,
             InputIterator last)
{
    return select_model<Criterion>(N, ordfirst, first, last, null_output());
}

/**
 * @}
 */

#endif /* BURG_HPP */
