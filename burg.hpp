// Except for any way in which it interferes with Cedrick Collomb's 2009
// copyright assertion in his article refactored from of Cedrick Collomb's
// 2009 article "Burg’s Method, Algorithm and Recursion":
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BURG_ALGORITHM_HPP
#define BURG_ALGORITHM_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <iterator>
#include <vector>

/**
 * Use Burg's recursion to find coefficients \f$a_i\f$ such that the sum
 * of the squared errors in both the forward linear prediction \f$x_n =
 * \sum_{i=1}^m a_i x_{n-i}\f$ and backward linear prediction \f$x_n =
 * \sum_{i=1}^m a_i x_{n+i}\f$ are minimized.  Input data \f$\vec{x}$
 * is taken from the range [data_first, data_last) in a single pass.
 *
 * Coefficients \f$\vec{a}\f$ are stored in [coeffs_first, coeffs_last) with
 * the model order determined by both <tt>k = distance(coeffs_first,
 * coeffs_last)</tt> and the \c all_models flag.  If \c all_models is
 * false, the coefficients for an AR(<tt>k</tt>) process are output.  If \c
 * all_models is true, the <tt>m*(m+1)/2</tt> coefficients for models
 * AR(1), AR(2), ..., AR(m) up * to order <tt>m = floor(sqrt(2*k))</tt> are
 * output.  One mean squared discrepancy value is also output for each model.
 *
 * The implementation has been
 * refactored from of Cedrick Collomb's 2009 article <a
 * href="http://www.emptyloop.com/technotes/A%20tutorial%20on%20Burg's%20method,%20algorithm%20and%20recursion.pdf">"Burg’s
 * Method, Algorithm and Recursion"</a>.  In particular, iterators are
 * employed, the working precision depends on the output coefficients
 * precision, output \f$\vec{a}\f$ has the opposite sign from Collomb's
 * notation to better match other authors, a mean squared discrepancy
 * calculation has been added, some loop index transformations have been
 * performed, and all lower order models may be output during the recursion
 * using \c all_models.
 *
 * @returns the number data values processed within [data_first, data_last).
 */
template<class InputIterator, class ForwardIterator, class OutputIterator>
std::size_t burg_algorithm(InputIterator   data_first,
                           InputIterator   data_last,
                           ForwardIterator coeffs_first,
                           ForwardIterator coeffs_last,
                           OutputIterator  msd_first,
                           const bool all_models = false)
{
    using std::distance;
    using std::fill;
    using std::inner_product;
    using std::iterator_traits;
    using std::negate;
    using std::numeric_limits;
    using std::sqrt;
    using std::transform;

    // ForwardIterator::value_type determines the working precision
    typedef typename iterator_traits<ForwardIterator>::value_type value;
    typedef typename std::vector<value> vector;
    typedef typename vector::size_type size;

    // Initialize f from [data_first, data_last), b by copying f, and data size
    vector f(data_first, data_last), b(f);
    const size N = f.size();

    // Get the order or maximum order of autoregressive model(s) desired.
    // When all_models is tru, the maximum order is the index of the largest
    // triangular number that will fit within [coeffs_first, coeffs_last).
    const size m = all_models ? sqrt(2*distance(coeffs_first, coeffs_last))
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
        if (all_models || kp1 == m)
        {
            // Output negated coefficients and the mean squared discrepancy
            coeffs_first = transform(Ak.begin() + 1, Ak.begin() + (kp1 + 1),
                                     coeffs_first, negate<value>());
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

    // Defensively NaN any remaining locations within the coeffs range
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

    template<typename T> void operator=(const T&) { }

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
template<class InputIterator, class ForwardIterator>
std::size_t burg_algorithm(InputIterator   data_first,
                           InputIterator   data_last,
                           ForwardIterator coeffs_first,
                           ForwardIterator coeffs_last)
{
    return burg_algorithm(data_first, data_last,
                          coeffs_first, coeffs_last,
                          null_output_iterator(), false);
}

#endif /* BURG_ALGORITHM_HPP */
