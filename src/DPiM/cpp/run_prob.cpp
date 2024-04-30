// g++ run_prob.cpp -o run_prob
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>

#include "Combinatorics.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Returns the probability of a run of at least r consecutive successes in n" << endl;
    cout << "independent trials where the probability of success at any one trial is p." << endl;
    cout << "Usage: run_prob [r] [n] [p]\n";
    return 0;
  }
  int r = atoi(argv[1]);
  int n = atoi(argv[2]);
  double p = atof(argv[3]);

  //double rp = Combinatorics::run_prob(r, n, p);
  double rp = Combinatorics::run_prob(r, n, p);
  cout << rp << endl;


  return 0;
}
