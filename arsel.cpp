// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "burg.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

template <class Real, class In>
std::size_t process(In& in,
                    Real& mean,
                    std::size_t &order,
                    std::vector<Real>& params,
                    std::vector<Real>& sigma2e,
                    std::vector<Real>& gain,
                    std::vector<Real>& autocor,
                    bool subtract_mean)
{
    // Use burg_method to estimate a hierarchy of AR models from input data
    params .reserve(order*(order + 1)/2);
    sigma2e.reserve(order + 1);
    gain   .reserve(order + 1);
    autocor.reserve(order + 1);
    const std::size_t N = burg_method(std::istream_iterator<Real>(in),
                                      std::istream_iterator<Real>(),
                                      mean,
                                      order,
                                      std::back_inserter(params),
                                      std::back_inserter(sigma2e),
                                      std::back_inserter(gain),
                                      std::back_inserter(autocor),
                                      subtract_mean,
                                      /* output hierarchy? */ true);

    // Find the best model according to CIC accounting for subtract_mean.
    typename std::vector<Real>::difference_type best;
    if (subtract_mean) {
        best = select_model<CIC<Burg<mean_subtracted> > >(
                    N, 0u, sigma2e.begin(), sigma2e.end());
    } else {
        best = select_model<CIC<Burg<mean_retained> > >(
                    N, 0u, sigma2e.begin(), sigma2e.end());
    }

    // Trim away everything but the best model (ranges might overlap)
    std::copy_backward(params.begin() + best*(best+1)/2,
                       params.begin() + best*(best+1)/2 + best,
                       params.begin() + best);
    params.resize(best);
    sigma2e[0] = sigma2e[best]; sigma2e.resize(1);
    gain   [0] = gain   [best]; gain   .resize(1);
    autocor.resize(best + 1);

    return N;
}

// Provides nice formatting of real-valued quantities to maximum precision
template<class CharT, class Traits, class Number>
static std::basic_ostream<CharT,Traits>& append_real(
        std::basic_ostream<CharT,Traits>& os,
        Number value)
{
    // Compute the displayed precision:
    //     Magic 2 is the width of a sign and a decimal point
    //     Magic 3 is the width of a sign, leading zero, and decimal point
    using std::pow;
    static const int append_prec  = std::numeric_limits<Number>::digits10;
    static const int append_width = append_prec + 5;
    static const double fixedmax = pow(10.0,   append_width - append_prec - 2 );
    static const double fixedmin = pow(10.0, -(append_width - append_prec - 3));

    // Format in fixed or scientific form as appropriate in given width
    // Care taken to not perturb observable ostream state after function call
    std::ios::fmtflags savedflags;
    std::streamsize savedprec;
    if (value >= fixedmin && value <= fixedmax) {
        savedflags = os.setf(std::ios::fixed      | std::ios::right,
                             std::ios::floatfield | std::ios::adjustfield);
        savedprec = os.precision(append_prec);
    } else {
        savedflags = os.setf(std::ios::scientific | std::ios::right,
                             std::ios::floatfield | std::ios::adjustfield);
        savedprec = os.precision(append_width - 9);
    }
    os << std::setw(append_width) << value;
    os.precision(savedprec);
    os.setf(savedflags);

    return os;
}


// Test burg_method against synthetic data
int main(int argc, char *argv[])
{
    using namespace std;

    // Process a possible --subtract-mean flag shifting arguments if needed
    bool subtract_mean = false;
    if (argc > 1 && 0 == strcmp("--subtract-mean", argv[1])) {
        subtract_mean = true;
        argv[1] = argv[0];
        ++argv;
        --argc;
    }

    double mean;
    size_t order = 512;
    vector<double> params, sigma2e, gain, autocor;
    size_t N = process(cin, mean, order, params,
                       sigma2e, gain, autocor, subtract_mean);

    append_real(cout << "# N                   ", N                 ) << endl;
    append_real(cout << "# Mean                ", mean              ) << endl;
    append_real(cout << "# \\sigma^2_\\epsilon ", sigma2e[0]        ) << endl;
    append_real(cout << "# Gain                ", gain[0]           ) << endl;
    append_real(cout << "# \\sigma^2_x         ", gain[0]*sigma2e[0]) << endl;

    for (size_t i = 0; i < params.size(); ++i)
        append_real(cout, params[i]) << endl;

    return 0;
}
