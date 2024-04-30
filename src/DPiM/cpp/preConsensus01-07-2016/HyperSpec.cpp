/* 
 * File:   HyperSpec.cpp
 * Author: jmintser
 *
 * Created on July 20, 2010, 3:04 PM
 */

#include "BaitByPrey.hpp"
#include "PreyByPrey.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <ctime>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void testJson(BaitByPrey sbp_mat, string outFile);

/*
 * 
 */
int main(int argc, char** argv) {
  /* Default arguments */
  unsigned int tim = static_cast<unsigned int>(std::time(0));
  string protlenDefault = "/home/glocke/DPiM/interfly_fbgn_avgTryp.out";
  po::options_description desc("Usage");
  desc.add_options()
    ("in", po::value<string>()->required(),
     "REQUIRED: APMS input data")
    ("out", po::value<string>()->required(),
     "REQUIRED: write HGScore here (JSON)")
    ("protlen", po::value<string>()->default_value(protlenDefault),
     "tab delimited length correction file for NSAF")
    ("sim", po::value<bool>()->default_value(false),
     "if true, scramble the input data for simulated HGScores")
    ("seed", po::value<unsigned int>()->default_value(tim),
     "random seed")
    ;

  po::variables_map opts;
  po::store(po::parse_command_line(argc, argv, desc), opts);
  try {
    po::notify(opts);
  } catch (std::exception& e) {
    std::cerr << desc << std::endl;
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  string inFile = opts["in"].as<string>();
  string outFile = opts["out"].as<string>();
  string lengthCorrectionFile = opts["protlen"].as<string>();
  bool simulateFlag = opts["sim"].as<bool>();
  srand( opts["seed"].as<unsigned int>() );

  BaitByPrey sbp_mat = BaitByPrey(inFile, lengthCorrectionFile);
  
  sbp_mat.reScale(BaitByPrey::NSAF);
  sbp_mat.reScale(BaitByPrey::Sqrt);
  sbp_mat.reScale(BaitByPrey::Round);

  //// SIMULATED HYPER-SCORE ////
  if (simulateFlag)
    sbp_mat.simulateDataset();

  /*
    cout << sbp_mat.sumSum() << " sum(sum( searchByPrey ) )" << endl;
    cout << sbp_mat.maxMax() << " max(max( searchByPrey ) )" << endl;
    cout << sbp_mat.sumSum() << " sum(sum( searchByPrey ) )" << endl;
    cout << sbp_mat.maxMax() << " max(max( searchByPrey ) )" << endl;
    cout << pbp_mat.sumSum() << " sum(sum( preyByPrey   ) )" << endl;
  */

  PreyByPrey pbp_mat = PreyByPrey(sbp_mat);

  //pbp_mat.print(); DEBUG
  //return 0;

  pbp_mat.computeHyge();
  pbp_mat.json(outFile);

  return (EXIT_SUCCESS);
}

void testJson(BaitByPrey sbp_mat, string outFile) {
  
  PreyByPrey pbp_mat = PreyByPrey(sbp_mat);

  size_t nProt = 10;
  vector<double> dummyRow;
  for (size_t i=0; i < nProt; ++i) {
    dummyRow.push_back( (i + 0.1) / 4.0);
  }
  
  pbp_mat.preyByPrey.resize(nProt, dummyRow);
  pbp_mat.preyByPreyHyge.resize(nProt, dummyRow);
  pbp_mat.json(outFile);
}
