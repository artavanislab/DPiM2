/* 
 * File:   Histogram.hpp
 * Author: jmintser
 *
 * Created on August 24, 2010, 5:00 PM
 */

#ifndef _HISTOGRAM_HPP
#define	_HISTOGRAM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <gsl/gsl_histogram.h>

using namespace std;

class Histogram {
public:
    size_t firstBin;
    size_t lastBin;
    size_t nBins;
    gsl_histogram * h;
    gsl_histogram_pdf * p;

    Histogram();
    Histogram(map<int,int> &inputMap);
    Histogram(vector<int> &inputVec);
    size_t numBins();
    void print ();
    int sample ();

    virtual ~Histogram();
private:

};

#endif	/* _HISTOGRAM_HPP */

