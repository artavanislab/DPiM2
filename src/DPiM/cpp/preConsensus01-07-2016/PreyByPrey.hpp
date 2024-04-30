/* 
 * File:   PreyByPrey.hpp
 * Author: jmintser
 *
 * Created on July 21, 2010, 12:23 PM
 */

#ifndef _PREYBYPREY_HPP
#define	_PREYBYPREY_HPP

#include "BaitByPrey.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <numeric>
#include <algorithm>
#include <limits>
#include <math.h>

using namespace std;

class BaitByPrey;

class PreyByPrey {
public:
  map<string,int> preyRef2Idx;
  vector<string> uniqPreyReferenceVec;
  set<string> directBait_PreyPair;
  vector< vector<double> > preyByPrey;
  vector< vector<double> > preyByPreyHyge;
  vector<double> preyByPreySum; // sum of all TSC's for each protein
  
  int shuffled;
  multimap<double,string> scoreMMap; 
  //map<string, double> scoreMap;
  //  map<string, vector<int> > scoreMapIndices;
  //  map<string, vector<double> > extraScoreMap;

  //multimap<double,string> simulatedScoreMMap;
  //map<string, vector<double> > simulatedScoreMap;
    
  PreyByPrey(BaitByPrey &bbp_mat);
    
  void makePBPmatrix(BaitByPrey &bbp_mat);
  double sumSum();
  double rowSum(int row, int start, int end);
  double colSum(int col, int start, int end);
  void computeHyge();
  void outputScores();
  void json(string outFile);
  double vectorMean(vector<double> &vec);
  double vectorStd(vector<double> &vec);
  void print ();

  virtual ~PreyByPrey();
private:
  double maxScore; // how to handle a hypergeom test that returns zero
};

#endif	/* _PREYBYPREY_HPP */

