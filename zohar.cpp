// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file
 * Test \ref ar::zohar_linear_solve routine against a reference solution.
 */

#include "ar.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>

int main()
{
    using namespace ar;
    using namespace std;

    vector<double> a;
    a.push_back( 2);
    a.push_back( 3);
    a.push_back( 5);
    a.push_back( 7);
    a.push_back(11);
    a.push_back(13);
    a.push_back(17);

    cout << "Topmost row of Toeplitz matrix is: \n\t1 ";
    copy(a.begin(), a.end(), ostream_iterator<double>(cout," "));
    cout << endl;

    vector<double> r;
    r.push_back(  2);
    r.push_back(  4);
    r.push_back(  8);
    r.push_back( 16);
    r.push_back( 32);
    r.push_back( 64);
    r.push_back(128);

    assert(a.size() == r.size());

    cout << "Leftmost column of Toeplitz matrix is: \n\t1 ";
    copy(r.begin(), r.end(), ostream_iterator<double>(cout," "));
    cout << endl;

    vector<double> d;
    d.push_back(1);
    d.push_back(2);
    d.push_back(3);
    d.push_back(4);
    d.push_back(5);
    d.push_back(6);
    d.push_back(7);
    d.push_back(8);

    assert(d.size() == a.size() + 1);

    cout << "Right hand side data is:\n\t";
    copy(d.begin(), d.end(), ostream_iterator<double>(cout," "));
    cout << endl;

    vector<double> exp;
    exp.push_back(-17./27);
    exp.push_back(  4./27);
    exp.push_back( 32./ 9);
    exp.push_back(- 5./ 3);
    exp.push_back(  0.   );
    exp.push_back(- 2.   );
    exp.push_back(- 1.   );
    exp.push_back(  2.   );

    assert(exp.size() == d.size());

    cout << "Expected solution is:\n\t";
    copy(exp.begin(), exp.end(), ostream_iterator<double>(cout," "));
    cout << endl;

    cout << "Solution computed by zohar_linear_solve is:\n\t";
    zohar_linear_solve(a.begin(), a.end(), r.begin(), d.begin());
    copy(d.begin(), d.end(), ostream_iterator<double>(cout," "));
    cout << endl;

    vector<double> err(exp.size());
    for (size_t i = 0; i < err.size(); ++i) {
        err[i] = exp[i] - d[i];
    }

    cout << "Term-by-term errors are:\n\t";
    copy(err.begin(), err.end(), ostream_iterator<double>(cout," "));
    cout << endl;

    double abserr = 0;
    for (size_t i = 0; i < exp.size(); ++i) abserr += std::abs(err[i]);
    cout << "Sum of the absolute errors is:\n\t" << abserr << endl;

    return 0;
}
