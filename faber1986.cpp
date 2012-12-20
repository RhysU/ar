/** @file
 * A test harness for Andersen's Burg algorithm variant as implemented by
 * Faber.
 */

/** Maximum number of input data points */
#define MAXSIZE (100000)

/** Maximum model order to fit */
#define MAXORD  (50)

/**
 * This function implements Burg's block data processing algorithm to compute
 * lattice [k] parameters.  The autoregressive (a) parameters are also computed
 * and returned (after the final iteration only).  This method minimizes the
 * sum of the forward and backward squared errors, subject to the constraint on
 * the AR parameters of the Levinson/Durbin algorithm.  This program is based
 * on equations in:
 *     S. M. Kay and S. L. Marple, "Spectrum Analysis--
 *     A Modern Perspective."  Proc. IEEE, Vol. 69,
 *     No. 11, Nov. 1981, pp. 1380-1419.
 * Equation numbers are noted in comments.  Note that eq. (2.74) for the
 * denominator recursion is wrong in that paper; it is corrected here.
 *
 * Adapted from L. J. Faber, "Commentary on the denominator recursion for
 * Burg's block algorithm," Proc.  IEEE, Vol. 74, No. 7, Jul. 1986, pp.
 * 1046-1047.
 *
 * Relative to the code appearing in [Faber1986]:
 * \li The maximum input size and order and been modified.
 * \li Double precision is used throughout.  ANSI C is used.
 * \li The method has been renamed to \ref faber1986.
 *
 * @param[in]  data input data to analyze
 * @param[in]  ndat number of data points
 * @param[out] k    reflection coefficients
 * @param[out] a    autoregressive coefficients
 * @param[out] err  AR prediction error energy
 * @param[in]  p    system order (# of coefficients)
 */
static void faber1986(double data[],
                      int    ndat,
                      double k[],
                      double a[],
                      double err[],
                      int p)
{
    int i;                    /* order index, i <= i <= p */
    int j;                    /* order sub-index, j < i */
    int n;                    /* time index, 0 <= n < ndat */
    double e[MAXSIZE];        /* forward prediction error[n] */
    double b[MAXSIZE];        /* backward pred. error[n] */
    double den;               /* denominator, eq. (2.74) */
    double num;               /* numerator, eq. (2.73) */
    double dk;                /* double prec. holder for k[i] */
    double ta[MAXORD];        /* temporary holder for a[i] */

    /******* Initialize *******/
    a[0] = 1;
    err[0] = k[0] = dk = 0;
    for (n = 0; n < ndat; n++) {
        e[n] = b[n] = data[n];
        err[0] += data[n] * data[n];
    }
    den = 2 * err[0];

    /******* Order Recursion *******/
    for (i = 1; i <= p; i++) {
        /**** Compute denom., corrected eq. (2.74) ****/
        den = den * (1 - dk * dk) - e[i-1] * e[i-1]
               - b[ndat - 1] * b[ndat - 1];
        /**** Compute k[i] eq. (2.73) ****/
        num = 0;
        for (n = i; n < ndat; n++)
            num += b[n-1] * e[n];
        num *= -2;
        k[i] = dk = num / den;
        /**** Update b[n], e[n], eq. (2.54), (2.52) ****/
        for (n=ndat - 1; n >= i; n--) {
            b[n] = b[n-1] + dk * e[n];
            e[n] = e[n] + dk * b[n-1];
        }
        /**** Update a[j], eq. (2.45) ****/
        for (j=1; j<i; j++)
            ta[j] = a[j] + dk*a[i-j];
        for (j=1; j<i; j++)
            a[j]=ta[j];
        a[i] = dk;
        /**** Update err[i], eq. (2.46) ****/
        err[i] = err[i-1] * (1 - dk * dk);
    }
}       /* ends faber1986 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

using namespace std;

/** Fit data from standard input using \ref faber1986. */
int main( int argc, char* argv[] )
{
    int p;
    if (argc < 2) {
        cerr << "Missing mandatory order argument" << endl;
        return 1;
    }
    p = atoi(argv[1]);
    if (p > MAXORD) {
        cerr << "Requested order exceeds limit MAXORD = " << MAXORD << endl;
        return 1;
    }

    // Load data from cin
    vector<double> data;
    copy(istream_iterator<double>(cin), istream_iterator<double>(), back_inserter(data));
    if (data.size() > MAXSIZE) {
        cerr << "Input data size exceeds limit MAXSIZE = " << MAXSIZE << endl;
        return 1;
    }

    // Compute AR model of given order
    vector<double> a(p+1), k(p+1), err(p+1);
    faber1986(&data[0], data.size(), &k[0], &a[0], &err[0], p);

    // Output them
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << showpos;
    copy(a.begin(), a.end(), ostream_iterator<double>(cout, "\n"));

    return 0;
}
