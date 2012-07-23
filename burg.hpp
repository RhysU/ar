// Except for any way in which it interferes with Cedrick Collomb's 2009
// copyright assertion in the article "Burg’s Method, Algorithm and Recursion":
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BURG_HPP
#define BURG_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <iterator>
#include <stdexcept>
#include <vector>

/**
 * Fit an autoregressive model to stationary, zero-mean time series data using
 * Burg's algorithm.  That is, assuming the model
 * \f{align}{
 *     x_n + a_1 x_{n - 1} + \dots + a_p x_{n - p} &= \epsilon_n
 *     &
 *     \epsilon_n &\sim{} N\left(0, \sigma^2_epsilon\right)
 *     \\
 *     \rho_k + a_1 \rho_{k-1} + \dots + a_p \rho_{k-p} &= 0
 *     &
 *     k &\geq{} p
 * \f}
 * find coefficients \f$a_i\f$ such that the sum of the squared errors in both
 * the forward prediction \f$x_n = \sum_{i=1}^m - a_i x_{n-i}\f$ and backward
 * prediction \f$x_n = \sum_{i=1}^m - a_i x_{n+i}\f$ are minimized.  Either a
 * single model of given order or a hierarchy of models up to and including a
 * maximum order may returned.
 *
 * The input data \f$\vec{x}\f$, which should have zero mean, is read from
 * <tt>[data_first,data_last)</tt> in a single pass.  The estimated model
 * parameters \f$a_i\f$ are output into <tt>[params_first, params_last)</tt>
 * with the behavior determined by both <tt>k = distance(params_first,
 * params_last)</tt> and the \c hierarchy flag:
 * <ul>
 *     <li>If \c hierarchy is \c false, only the \f$a_1, \dots, a_k\f$
 *         parameters for an AR(<tt>k</tt>) process are output.</li>
 *     <li>If \c hierarchy is true, the <tt>m*(m+1)/2</tt> coefficients for
 *         models AR(1), AR(2), ..., AR(m) up to and including order <tt>m =
 *         floor(sqrt(2*k))</tt> are output.</li>
 * </ul>
 * Note that the latter case is \e always computed;  The \c hierarchy flag
 * merely controls what is output.
 *
 * One mean squared discrepancy \f$\sigma^2_\epsilon\f$, also called the
 * innovation variance, and gain \f$\sigma^2_x / \sigma^2_\epsilon\f$ are
 * output for each model using \c sigma2e_first and \c gain_first.  The
 * autocorrelations for lags up to the maximum model order are output using \c
 * cor_first.  When \c hierarchy is true, Only lags <tt>[1,m]</tt> should be
 * applied for a model AR(<tt>m</tt>) when \c hierarchy is true.  The lag zero
 * autocorrelation is always one and is never output.  Autocovariances may be
 * computed by multiplying the autocorrelations by \f$\text{gain}
 * \sigma^2_\epsilon\f$.
 *
 * The implementation has been refactored heavily from Cedrick Collomb's 2009
 * article <a
 * href="http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf">"Burg’s
 * Method, Algorithm and Recursion"</a>.  In particular, iterators are
 * employed, the working precision depends on the output parameter precision,
 * the mean squared discrepancy calculation has been added, some loop index
 * transformations have been performed, and all lower order models may be
 * output during the recursion using \c hierarchy.  Gain and autocorrelation
 * calculations have been added based on sections 5.2 and 5.3 of Broersen, P.
 * M.  T. Automatic autocorrelation and spectral analysis. Springer, 2006.
 * http://dx.doi.org/10.1007/1-84628-329-9.
 *
 * @returns the number data values processed within [data_first, data_last).
 */
template <class InputIterator,
          class ForwardIterator1,
          class OutputIterator1,
          class OutputIterator2,
          class ForwardIterator2>
std::size_t burg_algorithm(InputIterator    data_first,
                           InputIterator    data_last,
                           ForwardIterator1 params_first,
                           ForwardIterator1 params_last,
                           OutputIterator1  sigma2e_first,
                           OutputIterator2  gain_first,
                           ForwardIterator2 cor_first,
                           const bool hierarchy = false)
{
    using std::copy;
    using std::distance;
    using std::fill;
    using std::inner_product;
    using std::numeric_limits;
    using std::sqrt;

    // ForwardIterator1::value_type determines the working precision
    typedef typename std::iterator_traits<ForwardIterator1>::value_type value;
    typedef typename std::vector<value> vector;
    typedef typename vector::size_type size;

    // Initialize f from [data_first, data_last), b by copying f, and data size
    // TODO Use pairwise summation here to stably compute and subtract mean
    vector f(data_first, data_last), b(f);
    const size N = f.size();

    // Get the order or maximum order of autoregressive model(s) desired.
    // When hierarchy is true, the maximum order is the index of the largest
    // triangular number that will fit within [params_first, params_last).
    const size m = hierarchy ? sqrt(2*distance(params_first, params_last))
                             :        distance(params_first, params_last);

    // Initialize mean squared discrepancy sigma2e and Dk
    value sigma2e = inner_product(f.begin(), f.end(), f.begin(), value(0));
    value Dk = - f[0]*f[0] - f[N - 1]*f[N - 1] + 2*sigma2e;
    sigma2e /= N;

    // Initialize recursion
    vector Ak(m + 1, value(0));
    Ak[0] = 1;
    value gain = 1;
    ForwardIterator2 cor(cor_first);

    // Perform Burg recursion
    for (size kp1 = 1; kp1 <= m; ++kp1)
    {
        // Compute mu from f, b, and Dk and then update sigma2e and Ak using mu
        // Afterwards, Ak[1:kp1] contains AR(k) coefficients by the recurrence
        // By the recurrence, Ak[kp1] will also be the reflection coefficient
        const value mu = 2/Dk*inner_product(f.begin() + kp1, f.end(),
                                            b.begin(), value(0));
        sigma2e *= (1 - mu*mu);
        for (size n = 0; n <= kp1/2; ++n)
        {
            const value t1 = Ak[n] - mu*Ak[kp1 - n];
            const value t2 = Ak[kp1 - n] - mu*Ak[n];
            Ak[n] = t1;
            Ak[kp1 - n] = t2;
        }

        // Update the gain per Broersen 2006 equation (5.25)
        gain *= 1 / (1 - Ak[kp1]*Ak[kp1]);

        // Compute and output the next autocorrelation coefficient
        // See Broersen 2006 equations (5.28) and (5.31) for details
        *cor++ = -inner_product(Ak.rend()-kp1, Ak.rend()-1, cor_first, Ak[kp1]);

        // Output parameters and the input and output variances when requested
        if (hierarchy || kp1 == m)
        {
            params_first = copy(Ak.begin() + 1, Ak.begin() + kp1 + 1,
                                params_first);
            *sigma2e_first++ = sigma2e;
            *gain_first++    = gain;
        }

        // Update f, b, and then Dk for the next iteration if another remains
        if (kp1 < m)
        {
            for (size n = 0; n < N - kp1; ++n)
            {
                const value t1 = f[n + kp1] - mu*b[n];
                const value t2 = b[n] - mu*f[n + kp1];
                f[n + kp1] = t1;
                b[n] = t2;
            }
            Dk = (1 - mu*mu)*Dk - f[kp1]*f[kp1] - b[N - kp1 - 1]*b[N - kp1 - 1];
        }
    }

    // Defensively NaN any unspecified locations within the coeffs range
    if (numeric_limits<value>::has_quiet_NaN)
        fill(params_first, params_last, numeric_limits<value>::quiet_NaN());

    // Return the number of values processed in [data_first, data_last)
    return N;
}

/**
 * Convert AR(p) process parameters \f$a_i\f$ into reflection coefficients
 * \f$k_i\f$.  The parameters are defined by
 * \f[
 *   x_n + a_0 x_{n-1} + \dots + a_{p-1} x_{n - (p + 1)} = \epsilon_n
 * \f]
 * while the reflection coefficients are the negative of the partial
 * autocorrelations as defined within the classical Levinson-Durbin recursion.
 * The model order is determined by <tt>p = distance(first, last)</tt>.
 *
 * The in-place conversion algorithm is from section 5.4.3 of Broersen, P. M.
 * T. Automatic autocorrelation and spectral analysis. Springer, 2006.
 * http://dx.doi.org/10.1007/1-84628-329-9.
 *
 * @param[in,out] first On input, the beginning of the range containing
 *                      the parameters.  On output, the range contains
 *                      the reflection coefficients.
 * @param[in,out] last  The exclusive end of the range.
 *
 * @return the variance of the process normalized by the variance of the
 *         input noise \f$\sigma_epsilon\f$ where \f$\epsilon \sim{}
 *         N(0,\sigma_epsilon)\f$.  That is, \f$\prod_{i=0}^{p-1}
 *         \left(1-k_i^2\right)^{-1}\f$ is returned.
 */
template <class BidirectionalIterator>
typename std::iterator_traits<BidirectionalIterator>::value_type
reflection_coefficients(BidirectionalIterator first,
                        BidirectionalIterator last)
{
    using std::distance;
    using std::reverse_iterator;

    typedef BidirectionalIterator iterator;
    typedef reverse_iterator<iterator> reverse;
    typedef typename std::iterator_traits<iterator>::value_type value;
    typedef typename std::iterator_traits<iterator>::difference_type difference;

    // Determine problem size using [first,last) and ensure nontrivial
    const difference n = distance(first, last);
    if (n < 1) throw std::invalid_argument("distance(first, last) < 1");

    // Initialize and recurse.  Output is *rfirst on every iteration.
    value gain = 1;
    reverse rfirst(last);
    for (difference i = n; i --> 1 ;) {

        // Preserve the final parameter as the current reflection coefficient
        const value k = *rfirst++;
        const value mu = 1 / (1 - k*k);
        gain *= mu;

        // Compute the recursive inputs by traversing from both ends
        // Front write occurs second so it "wins" in odd-length iterations
        iterator front(first);
        reverse  back(rfirst);
        for (difference j = 0; j <= (i-1)/2; ++j) {
            const value tf = *front;
            const value tb = *back;
            *back++  = mu*(tb - k*tf);
            *front++ = mu*(tf - k*tb);
        }

    }

    return gain;
}

/**
 * Convert AR(p) reflection coefficients \f$k_i\f$ into process parameters
 * \f$a_i\f$.  The parameters are defined by
 * \f[
 *   x_n + a_0 x_{n-1} + \dots + a_{p-1} x_{n - (p + 1)} = \epsilon_n
 * \f]
 * while the reflection coefficients are the negative of the partial
 * autocorrelations as defined within the classical Levinson-Durbin recursion.
 * The model order is determined by <tt>p = distance(first, last)</tt>.
 *
 * The in-place conversion algorithm is from section 5.4.2 of Broersen, P. M.
 * T. Automatic autocorrelation and spectral analysis. Springer, 2006.
 * http://dx.doi.org/10.1007/1-84628-329-9.
 *
 * @param[in,out] first On input, the beginning of the range containing
 *                      the reflection coefficients.  On output, the range
 *                      contains the parameters.
 * @param[in,out] last  The exclusive end of the range.
 */
template <class BidirectionalIterator>
void parameters(BidirectionalIterator first,
                BidirectionalIterator last)
{
    using std::distance;
    using std::reverse_iterator;

    typedef BidirectionalIterator iterator;
    typedef reverse_iterator<iterator> reverse;
    typedef typename std::iterator_traits<iterator>::value_type value;
    typedef typename std::iterator_traits<iterator>::difference_type difference;

    // Determine problem size using [first,last) and ensure nontrivial
    const difference n = distance(first, last);
    if (n < 1) throw std::invalid_argument("distance(first, last) < 1");

    // Initialize, hard coding result in simple p = 1 case, and then recurse.
    iterator k(first);
    ++k;
    for (difference i = 1; i < n; ++i, ++k) {

        // Compute the recursive inputs by traversing from both ends
        // Front write occurs second so it "wins" in odd-length iterations
        iterator front(first);
        reverse  back(k);
        for (difference j = 0; j <= (i-1)/2; ++j) {
            const value tf = *front;
            const value tb = *back;
            *back++  = tb + (*k)*tf;
            *front++ = tf + (*k)*tb;
        }

    }
}

/**
 * Compute AR(p) lag \f$k < p\f$ autocorrelations \f$\rho_1,\dots,\rho_{p-1}\f$
 * given reflection coefficients \f$k_i\f$.  The reflection coefficients are
 * the negative of the partial autocorrelations as defined within the classical
 * Levinson-Durbin recursion.  The model order is determined by <tt>p =
 * distance(k_first, k_last)</tt>.  Autocorrelations for \f$k >= p\f$ may be
 * computed directly from \f$a_i\f$ using that
 * \f[
 *   \rho_k = - a_0 \rho_{k-1} + \dots + a_{p-1} x_{k - (p + 1)}
 * \f]
 * where \f$a_i\f$ are the corresponding process parameters.
 *
 * The algorithm is from section 5.4.4 of Broersen, P. M.  T. Automatic
 * autocorrelation and spectral analysis. Springer, 2006.
 * http://dx.doi.org/10.1007/1-84628-329-9.  Because the algorithm operates
 * in-place on <tt>[k_first, k_last)</tt> by incrementally building process
 * parameters of increasingly higher orders, it is highly subject to numerical
 * stability issues.
 *
 * @param[in,out] k_first   The beginning of the range containing reflection
 *                          coefficients from, e.g., reflection_coefficients().
 *                          Values <i>will be overwritten</i> by parameters.
 * @param[in,out] k_last    The exclusive end of the reflection coefficients.
 * @param[out]    rho_first The beginning of the computed autocorrelations.
 *                          The zeroth output is the lag 1 autocorrelation.
 */
template <class BidirectionalIterator,
          class ForwardIterator>
void autocorrelations(BidirectionalIterator k_first,
                      BidirectionalIterator k_last,
                      ForwardIterator       rho_first)
{
    using std::distance;
    using std::inner_product;
    using std::reverse_iterator;

    typedef BidirectionalIterator iterator;
    typedef reverse_iterator<iterator> reverse;
    typedef typename std::iterator_traits<iterator>::value_type value;
    typedef typename std::iterator_traits<iterator>::difference_type difference;

    // Determine problem size using [first,last) and ensure nontrivial
    const difference n = distance(k_first, k_last);
    if (n < 1) throw std::invalid_argument("distance(k_first, k_last) < 1");

    // Initialize, hard coding result in simple p = 1 case, and then recurse.
    ForwardIterator rho(rho_first);
    iterator k(k_first);
    *rho++ = -*(k++);
    for (difference i = 1; i < n; ++i, ++k, ++rho) {

        // Compute the recursive inputs by traversing from both ends
        // Front write occurs second so it "wins" in odd-length iterations
        iterator front(k_first);
        reverse  back(k);
        for (difference j = 0; j <= (i-1)/2; ++j) {
            const value tf = *front;
            const value tb = *back;
            *back++  = (*k)*tb + tf;
            *front++ = (*k)*tf + tb;
        }

        // Compute the k-th autocorrelation using Broersen's equation (5.31)
        *rho  = inner_product(rho_first, rho, k_first, *k);
        *rho *= -1;
    }
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
 * \f$vec{r}\f$ being the conjugate of \f$\vec{a}\f$.
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
 * @param[in]  d_first Beginning of the range containing \f$\vec{d}\f$.
 * @param[out] s_first Beginning of the output range to which
 *                     <strong><tt>n+1</tt></strong> entries will be
 *                     written.
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

    // OutputIterator::value_type determines the working precision
    typedef typename std::iterator_traits<OutputIterator>::value_type value;
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
    for (size i = 1; i < n; ++i) {

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
        for (size j = 0; j < i; ++j) {
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
        for (size j = 0; j < n; ++j) {
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
 * given \f$\vec{a}\f$, \f$\vec{r}\f$, and \f$\vec{d}\f$.  The dimension of
 * the problem is fixed by <tt>n = distance(a_first, a_last)</tt>.  A symmetric
 * Toeplitz solve can be performed by having \f$\vec{a}\f$ and \f$vec{r}\f$
 * iterate over the same data.  The Hermitian case requires two buffers with
 * \f$vec{r}\f$ being the conjugate of \f$\vec{a}\f$.
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

#endif /* BURG_HPP */
