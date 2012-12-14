#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

using namespace std;

// Returns in vector coefficients calculated using Burg algorithm applied to the input source data x
// Taken directly from http://www.emptyloop.com/technotes/
void BurgAlgorithm( vector<double>& coeffs, const vector<double>& x )
{
    // GET SIZE FROM INPUT VECTORS
    size_t N = x.size() - 1;
    size_t m = coeffs.size();

    // INITIALIZE Ak
    vector<double> Ak( m + 1, 0.0 );
    Ak[ 0 ] = 1.0;

    // INITIALIZE f and b
    vector<double> f( x );
    vector<double> b( x );

    // INITIALIZE Dk
    double Dk = 0.0;
    for ( size_t j = 0; j <= N; j++ )
    {
        Dk += 2.0 * f[ j ] * f[ j ];
    }
    Dk -= f[ 0 ] * f[ 0 ] + b[ N ] * b[ N ];

    // BURG RECURSION
    for ( size_t k = 0; k < m; k++ )
    {
        // COMPUTE MU
        double mu = 0.0;
        for ( size_t n = 0; n <= N - k - 1; n++ )
        {
            mu += f[ n + k + 1 ] * b[ n ];
        }
        mu *= -2.0 / Dk;

        // UPDATE Ak
        for ( size_t n = 0; n <= ( k + 1 ) / 2; n++ )
        {
            double t1 = Ak[ n ] + mu * Ak[ k + 1 - n ];
            double t2 = Ak[ k + 1 - n ] + mu * Ak[ n ];
            Ak[ n ] = t1;
            Ak[ k + 1 - n ] = t2;
        }

        // UPDATE f and b
        for ( size_t n = 0; n <= N - k - 1; n++ )
        {
            double t1 = f[ n + k + 1 ] + mu * b[ n ];
            double t2 = b[ n ] + mu * f[ n + k + 1 ];
            f[ n + k + 1 ] = t1;
            b[ n ] = t2;
        }

        // UPDATE Dk
        Dk = ( 1.0 - mu * mu ) * Dk - f[ k + 1 ] * f[ k + 1 ] - b[ N - k - 1 ] * b[ N - k - 1 ];
    }

    // ASSIGN COEFFICIENTS
    coeffs.assign( ++Ak.begin(), Ak.end() );
}

// Example program using Burg’s algorithm
int main( int argc, char* argv[] )
{
    int order;
    if (argc < 2) {
        cerr << "Missing mandatory order argument";
        return 1;
    }
    order = atoi(argv[1]);

    // Load data from cin
    vector<double> original;
    copy(istream_iterator<double>(cin), istream_iterator<double>(), back_inserter(original));

    // Compute AR model of given order
    vector<double> coeffs( order, 0.0 );
    BurgAlgorithm( coeffs, original );

    // Output them
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos;
    cout << 1.0 << endl;
    copy(coeffs.begin(), coeffs.end(), ostream_iterator<double>(cout, "\n"));

    return 0;
}
