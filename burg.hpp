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
#include <vector>

/**
 * Use Burg's recursion to find coefficients \f$a_i\f$ such that the sum
 * of the squared errors in both the forward linear prediction \f$x_n =
 * \sum_{i=1}^m - a_i x_{n-i}\f$ and backward linear prediction \f$x_n =
 * \sum_{i=1}^m - a_i x_{n+i}\f$ are minimized.  Input data \f$\vec{x}$
 * is taken from the range [data_first, data_last) in a single pass.
 *
 * Coefficients \f$\vec{a}\f$ are stored in [coeffs_first, coeffs_last)
 * with the model order determined by both <tt>k = distance(coeffs_first,
 * coeffs_last)</tt> and the \c hierarchy flag.  If \c hierarchy is
 * false, the coefficients for an AR(<tt>k</tt>) process are output.  If \c
 * hierarchy is true, the <tt>m*(m+1)/2</tt> coefficients for models AR(1),
 * AR(2), ..., AR(m) up to order <tt>m = floor(sqrt(2*k))</tt> are output.
 * The mean squared discrepancy value \f$\sigma_\epsilon^2\f$ is also output
 * for each model.  Each corresponding AR(p) prediction model has the form
 * \f[
 *   x_n + a_0 x_{n-1} + \dots + a_{p-1} x_{n - (p + 1)} = \epsilon_n
 * \f]
 * where \f$\epsilon_n\f$ has variance \f$\sigma_\epsilon^2\f$.
 *
 * The implementation has been refactored from of Cedrick Collomb's 2009
 * article <a
 * href="http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf">"Burg’s
 * Method, Algorithm and Recursion"</a>.  In particular, iterators are
 * employed, the working precision depends on the output coefficients
 * precision, a mean squared discrepancy calculation has been added, some loop
 * index transformations have been performed, and all lower order models may be
 * output during the recursion using \c hierarchy.
 *
 * @returns the number data values processed within [data_first, data_last).
 */
template <class InputIterator, class ForwardIterator, class OutputIterator>
std::size_t burg_algorithm(InputIterator   data_first,
                           InputIterator   data_last,
                           ForwardIterator coeffs_first,
                           ForwardIterator coeffs_last,
                           OutputIterator  msd_first,
                           const bool hierarchy = false)
{
    using std::copy;
    using std::distance;
    using std::fill;
    using std::inner_product;
    using std::iterator_traits;
    using std::numeric_limits;
    using std::sqrt;

    // ForwardIterator::value_type determines the working precision
    typedef typename iterator_traits<ForwardIterator>::value_type value;
    typedef typename std::vector<value> vector;
    typedef typename vector::size_type size;

    // Initialize f from [data_first, data_last), b by copying f, and data size
    vector f(data_first, data_last), b(f);
    const size N = f.size();

    // Get the order or maximum order of autoregressive model(s) desired.
    // When hierarchy is true, the maximum order is the index of the largest
    // triangular number that will fit within [coeffs_first, coeffs_last).
    const size m = hierarchy ? sqrt(2*distance(coeffs_first, coeffs_last))
                             :        distance(coeffs_first, coeffs_last);

    // Initialize mean squared discrepancy msd and Dk
    value msd = inner_product(f.begin(), f.end(), f.begin(), value(0));
    value Dk = - f[0]*f[0] - f[N - 1]*f[N - 1] + 2*msd;
    msd /= N;

    // Initialize Ak
    vector Ak(m + 1, value(0));
    Ak[0] = 1;

    // Burg recursion
    for (size kp1 = 1; kp1 <= m; ++kp1)
    {
        // Compute mu from f, b, and Dk and then update msd and Ak using mu
        const value mu = 2/Dk*inner_product(f.begin() + kp1, f.end(),
                                            b.begin(), value(0));
        msd *= (1 - mu*mu);
        for (size n = 0; n <= kp1/2; ++n)
        {
            const value t1 = Ak[n] - mu*Ak[kp1 - n];
            const value t2 = Ak[kp1 - n] - mu*Ak[n];
            Ak[n] = t1;
            Ak[kp1 - n] = t2;
        }

        // Here, Ak[1:kp1] contains AR(k) coefficients by the recurrence.
        if (hierarchy || kp1 == m)
        {
            // Output coefficients and the mean squared discrepancy
            coeffs_first = copy(Ak.begin() + 1, Ak.begin() + (kp1 + 1),
                                coeffs_first);
            *msd_first++ = msd;
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
    {
        fill(coeffs_first, coeffs_last, numeric_limits<value>::quiet_NaN());
    }

    // Return the number of values processed in [data_first, data_last)
    return N;
}


namespace { // anonymous

// Used to discard output per http://stackoverflow.com/questions/335930/
struct null_output_iterator
    : std::iterator< std::output_iterator_tag, null_output_iterator >
{

    template <typename T> void operator=(const T&) { }

    null_output_iterator& operator++() { return *this; }

    null_output_iterator operator++(int) {
        null_output_iterator it(*this);
        ++*this;
        return it;
    }

    null_output_iterator& operator*() { return *this; }
};

}

/**
 * Use Burg's recursion to find coefficients \f$a_i\f$ such that the sum
 * of the squared errors in both the forward linear prediction \f$x_n =
 * \sum_{i=1}^m a_i x_{n-i}\f$ and backward linear prediction \f$x_n =
 * \sum_{i=1}^m a_i x_{n+i}\f$ are minimized.  Input data \f$\vec{x}$
 * is taken from the range [data_first, data_last) in a single pass.
 *
 * Coefficients \f$\vec{a}\f$ are stored in [coeffs_first, coeffs_last) with
 * the model order determined by <tt>distance(coeffs_first, coeffs_last)</tt>.
 *
 * @returns the number data values processed within [data_first, data_last).
 */
template <class InputIterator, class ForwardIterator>
std::size_t burg_algorithm(InputIterator   data_first,
                           InputIterator   data_last,
                           ForwardIterator coeffs_first,
                           ForwardIterator coeffs_last)
{
    return burg_algorithm(data_first, data_last,
                          coeffs_first, coeffs_last,
                          null_output_iterator(), false);
}

/**
 * Find the solution to a Toeplitz set of linear equations.  That is,
 * satisfy the equation
 * \f[
 *      L_{n+1} s_{n+1} = d_{n+1}
 *      \mbox{ where }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    r_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{s}\f$, \f$\vec{r}\f$, and \f$\vec{d}\f$.  The dimension of the
 * problem is fixed by <tt>n = distance(a_first, a_last)</tt>.  Though less
 * efficient than a symmetric solution procedure, \f$\vec{a}\f$ and
 * \f$vec{r}\f$ may iterate over the same data.
 *
 * The algorithm is from Zohar, Shalhav. "The Solution of a
 * Toeplitz Set of Linear Equations." J. ACM 21 (April 1974):
 * 272-276. http://dx.doi.org/10.1145/321812.321822.  It has complexity
 * like <tt>O(2*(n+1)^2)</tt>.  See Bunch, James R. "Stability of Methods
 * for Solving Toeplitz Systems of Equations." SIAM Journal on Scientific and
 * Statistical Computing 6 (1985): 349-364. http://dx.doi.org/10.1137/0906025
 * for a discussion of the algorithms stability characteristics.
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
    // Tildes indicate transposes while hats indicate reversed vectors.

    using std::copy;
    using std::distance;
    using std::inner_product;
    using std::iterator_traits;
    using std::reverse_iterator;
    using std::transform;

    // OutputIterator::value_type determines the working precision
    typedef typename iterator_traits<OutputIterator>::value_type value;
    typedef typename std::vector<value> vector;
    typedef typename vector::size_type size;

    // Determine problem size using [a_first,a_last)
    const size n = distance(a_first, a_last);

    // Allocate working storage and set initial values for recursion:
    vector s;    s   .reserve(n+1); s   .push_back( *d_first);
    vector ehat; ehat.reserve(n+1); ehat.push_back(-a_first[0]);
    vector g;    g   .reserve(n+1); g   .push_back(-r_first[0]);
    value lambda  = 1 - a_first[0]*r_first[0];

    // Recursion for i = 1, 2, ..., n:
    for (size i = 1; i <= n; ++i) {

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
         * Using std::transform would have been nice but impossible
         * given the nature of the \hat{e} and g updates.
         *
         * s_{i+1} = \bigl(\begin{smallmatrix}
         *              s_i + (\theta_i/\lambda_i) \hat{e}_i \\
         *              \theta_i/\lambda_i
         *          \end{smallmatrix}\bigr)
         *
         * \hat{e}_{i+1} = \bigl(\begin{smallmatrix}
         *                     \eta_i/\lambda_i \\
         *                     \hat{e}_i + (\ega_i/\lambda_i) g_i
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
        for (size j = 0; j < i; ++j) {
            const value ejhat = ehat[j], gj = g[j];
            s[j]     += theta_by_lambda*ejhat;
            ehat[j]  +=   eta_by_lambda*gj;
            g[j]     += gamma_by_lambda*ejhat;
        }
        s   .push_back(theta_by_lambda);
        ehat.push_back(  eta_by_lambda);
        g   .push_back(gamma_by_lambda);

        // \lambda_{i+1} = \lambda_i - \eta_i \gamma_i / \lambda_i
        lambda -= neg_eta*neg_gamma/lambda;
    }

    // Output solution
    copy(s.begin(), s.end(), s_first);
}

namespace { // anonymous

// Provides an infinite input iterator returning zeros
template <typename T>
struct zero_input_iterator : std::iterator<
        std::input_iterator_tag, T, std::ptrdiff_t, const T*, const T&
    >
{
    T operator*() { return 0; }
    zero_input_iterator& operator++()    { return *this; }
    zero_input_iterator& operator++(int) { return *this; }
};

}

/**
 * Find the homogeneous solution to a Toeplitz set of linear equations.  That
 * is, satisfy the equation
 * \f[
 *      L_{n+1} s_{n+1} = 0
 *      \mbox{ where }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    r_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{s}\f$ and \f$\vec{r}\f$.  The dimension of the problem is
 * fixed by <tt>n = distance(a_first, a_last)</tt>.  Though less efficient than
 * a symmetric solution procedure, \f$\vec{a}\f$ and \f$vec{r}\f$ may iterate
 * over the same data.
 *
 * @param[in]  a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]  a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in]  r_first Beginning of the range containing \f$\vec{r}\f$.
 * @param[out] s_first Beginning of the output range to which
 *                     <strong><tt>n+1</tt></strong> entries will be
 *                     written.
 */
template<class RandomAccessIterator,
         class OutputIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        RandomAccessIterator r_first,
                        OutputIterator       s_first)
{
    // OutputIterator::value_type determines the working precision
    using std::iterator_traits;
    typedef typename iterator_traits<OutputIterator>::value_type value;
    return zohar_linear_solve(a_first, a_last, r_first,
                              zero_input_iterator<value>(), s_first);
}

#endif /* BURG_HPP */
