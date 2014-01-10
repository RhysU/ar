#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

#include "real.hpp"

using namespace std;

/** @file
 * A test harness for Cedrick Collomb's Burg algorithm variant.
 *
 * Taken from Cedrick Collomb. "Burg's method, algorithm, and recursion",
 * November 2009 available at http://www.emptyloop.com/technotes/.
 */

/**
 * Returns in vector coefficients calculated using Burg algorithm applied to
 * the input source data x
 */
static void BurgAlgorithm( vector<real>& coeffs, const vector<real>& x )
{
    // GET SIZE FROM INPUT VECTORS
    size_t N = x.size() - 1;
    size_t m = coeffs.size();

    // INITIALIZE Ak
    vector<real> Ak( m + 1, 0.0 );
    Ak[ 0 ] = 1.0;

    // INITIALIZE f and b
    vector<real> f( x );
    vector<real> b( x );

    // INITIALIZE Dk
    real Dk = 0.0;
    for ( size_t j = 0; j <= N; j++ )
    {
        Dk += 2.0 * f[ j ] * f[ j ];
    }
    Dk -= f[ 0 ] * f[ 0 ] + b[ N ] * b[ N ];

    // BURG RECURSION
    for ( size_t k = 0; k < m; k++ )
    {
        // COMPUTE MU
        real mu = 0.0;
        for ( size_t n = 0; n <= N - k - 1; n++ )
        {
            mu += f[ n + k + 1 ] * b[ n ];
        }
        mu *= -2.0 / Dk;

        // UPDATE Ak
        for ( size_t n = 0; n <= ( k + 1 ) / 2; n++ )
        {
            real t1 = Ak[ n ] + mu * Ak[ k + 1 - n ];
            real t2 = Ak[ k + 1 - n ] + mu * Ak[ n ];
            Ak[ n ] = t1;
            Ak[ k + 1 - n ] = t2;
        }

        // UPDATE f and b
        for ( size_t n = 0; n <= N - k - 1; n++ )
        {
            real t1 = f[ n + k + 1 ] + mu * b[ n ];
            real t2 = b[ n ] + mu * f[ n + k + 1 ];
            f[ n + k + 1 ] = t1;
            b[ n ] = t2;
        }

        // UPDATE Dk
        Dk = ( 1.0 - mu * mu ) * Dk - f[ k + 1 ] * f[ k + 1 ] - b[ N - k - 1 ] * b[ N - k - 1 ];
    }

    // ASSIGN COEFFICIENTS
    coeffs.assign( ++Ak.begin(), Ak.end() );
}

/** Fit data from standard input using \ref BurgAlgorithm. */
int main( int argc, char* argv[] )
{
    int order;
    if (argc < 2) {
        cerr << "Missing mandatory order argument" << endl;
        return 1;
    }
    order = atoi(argv[1]);

    // Load data from cin
    vector<real> data;
    copy(istream_iterator<real>(cin), istream_iterator<real>(), back_inserter(data));

    // Compute AR model of given order
    vector<real> coeffs( order );
    BurgAlgorithm( coeffs, data );

    // Output them, including a leading 1 not computed by BurgAlgorithm
    cout.precision(numeric_limits<real>::digits10 + 2);
    cout << showpos;
    cout << 1.0 << endl;
    copy(coeffs.begin(), coeffs.end(), ostream_iterator<real>(cout, "\n"));

    return 0;
}
