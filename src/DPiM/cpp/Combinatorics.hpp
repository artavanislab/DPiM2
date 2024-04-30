/* 
 * File:   Combinatorics.hpp
 * Author: jmintser
 *
 * Created on June 2, 2011, 6:06 PM
 */

#ifndef _COMBINATORICS_HPP
#define	_COMBINATORICS_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <numeric> // for accumulate
#include <functional> // for multiplies

using namespace std;

class Combinatorics {
public:

    // Count number of combinations C(n,k)
    static int nchoosek(int n, int k) {
        int c;
        if (k > n / 2)
            k = n - k;
        if (k <= 1) {
            c = (int) pow((double) n, (double) k);
        } else {
            vector<double> nums(k);
            for (unsigned int i = 0; i < nums.size(); i++) {
                nums[i] = (double) (n - k + 1 + i) / (i + 1);
            }
            double tmp = accumulate(nums.begin(), nums.end(), 1.0, multiplies<double>());
            c = static_cast<int> (tmp);
        }
        return c;
    }

    // Count log number of combinations C(n,k)
    // This can be used when n is large relative to k and the coefficient overflows

    static double ln_nchoosek(int n, int k) {
        double c;
        double ln_c;
        if (k > n / 2)
            k = n - k;
        if (k <= 1) {
            c = pow((double) n, (double) k);
            ln_c = log(c);
        } else {
            vector<double> nums(k);
            for (unsigned int i = 0; i < nums.size(); i++) {
                nums[i] = log((double) (n - k + 1 + i) / (i + 1));
            }
            ln_c = accumulate(nums.begin(), nums.end(), 0.0);
        }
        ////cout << "ln_choosek: " << n << " " << k << "\t" << ln_c << endl; ////
        return ln_c;
    }

    static inline int factorial(int a) {
        int result = 1;
        while (a > 1) {
            result *= a;
            a--;
        }
        return result;
    }

    /* Returns the probability of a run of at least r consecutive successes in n
     * independent trials where the probability of success at any one trial is p.
     * Need both functions below
     * Implemented following Mark B. Villarino, 'Probability of a Run', arxiv (2006)
     * NOTE: appears to break if compiles > -O1
     */
    static double run_prob(int r, int n, double p) {
        double Bnr = calc_beta_ln(r, n, p);
        double Bnrr = calc_beta_ln(r, (n - r), p);
        double zn = Bnr - pow(p, (double) r) * Bnrr;
        ////cout << "Bnr " << Bnr << "\t" << "Bnrr " << Bnrr << "\t" << "zn " << zn << "\t" << "yn " << yn << endl; ////
        if ((zn < 0) | (zn > 1))
            zn = 0;
        double yn = 1 - zn;
        return yn;
    }

    // special sub-function of &run_prob

    static double calc_beta_ln(int r, int n, double p) {
        double q = 1 - p;
        vector<double> Bnr_vec(floor((double) n / (r + 1)) + 1);
        for (unsigned int l = 0; l < Bnr_vec.size(); l++) {
            double tmp_ln = ln_nchoosek(n - l*r, l) + l * log(q * pow(p, (double) r));
            double tmp = exp(tmp_ln);
            Bnr_vec[l] = pow(-1.0, (double) l) * tmp;
        }
        /////cout << "Bnr_vec: "; copy (Bnr_vec.begin(), Bnr_vec.end(), ostream_iterator<double> (cout, " ")); cout << endl; ////
        double Bnr = accumulate(Bnr_vec.begin(), Bnr_vec.end(), 0.0);
        return Bnr;
    }

};

#endif	/* _COMBINATORICS_HPP */

