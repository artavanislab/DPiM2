/* 
 * File:   Histogram.cpp
 * Author: jmintser
 * 
 * Created on August 24, 2010, 5:00 PM
 */

#include "Histogram.hpp"

Histogram::Histogram() {
    h = gsl_histogram_alloc(1);
}

Histogram::Histogram(map<int, int> &inputMap) {
    map<int, int>::reverse_iterator rMapIt = inputMap.rbegin();
    firstBin = 0;
    lastBin = (*rMapIt).first;
    nBins = lastBin + 1;
    std::vector<double> range_vec(nBins + 1);
    for (unsigned int i = 0; i < range_vec.size(); i++) {
        range_vec[i] = i;
    }
    double* range_arr = &range_vec[0];

    h = gsl_histogram_alloc(nBins);
    gsl_histogram_set_ranges(h, range_arr, nBins + 1);

    for (map<int, int>::iterator mapIt = inputMap.begin(); mapIt != inputMap.end(); mapIt++) {
        gsl_histogram_accumulate(h, (*mapIt).first, (*mapIt).second);
    }

    p = gsl_histogram_pdf_alloc(nBins);
    gsl_histogram_pdf_init(p, h);
}

Histogram::Histogram(vector<int> &inputVec) {
    firstBin = 0;
    lastBin = inputVec.size() - 1;
    nBins = inputVec.size();
    std::vector<double> range_vec(nBins + 1);
    for (unsigned int i = 0; i < range_vec.size(); i++) {
        range_vec[i] = i;
    }
    double* range_arr = &range_vec[0];

    h = gsl_histogram_alloc(nBins);
    gsl_histogram_set_ranges(h, range_arr, nBins + 1);

    for (unsigned int i = 0; i < inputVec.size(); i++) {
        gsl_histogram_accumulate(h, i, inputVec[i]);
    }

    p = gsl_histogram_pdf_alloc(nBins);
    gsl_histogram_pdf_init(p, h);
}

size_t Histogram::numBins() {
    return gsl_histogram_bins(h);
}

void Histogram::print() {
//    FILE * f = fopen ("hist.dat", "wb");
//    gsl_histogram_fprintf(f, h, "%g", "%g");
//    fclose (f);
    gsl_histogram_fprintf(stdout, h, "%g", "%g");
}

int Histogram::sample() {
    double randFraction = rand() / ((double) RAND_MAX);
    int x = gsl_histogram_pdf_sample(p, randFraction);
    return x;
}

Histogram::~Histogram() {
    gsl_histogram_pdf_free(p);
    gsl_histogram_free(h);
}
