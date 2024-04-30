/* 
 * File:   HyperSpec.cpp
 * Author: jmintser, glocke
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

void debugPrint(BaitByPrey sbp_mat, string sid, string prey);
void debugPrint2(BaitByPrey sbp_mat, string sidFile, string prey1, string prey2);
  
/*
 * 
 */
int main(int argc, char** argv) {
  /* Default arguments */
  unsigned int tim = static_cast<unsigned int>(std::time(0));
  //string protlenDefault = "/home/glocke/DPiM/nsaf/dmel-all-translation-r6.09.aveLen.tsv";
  string protlenDefault = "/home/glocke/DPiM/nsaf/REVdmel-all-translation-r6.07_TAGS_sorted_trEl_vir.aveLen.tsv";

  po::options_description desc("Usage");
  desc.add_options()
    ("in", po::value<string>()->required(),
    //("in", po::value<string>()->default_value(""),
     "REQUIRED: APMS input data")
    ("out", po::value<string>()->required(),
     //("out", po::value<string>()->default_value(""),
     "REQUIRED: write HGScore here (JSON)")
    ("protlen", po::value<string>()->default_value(protlenDefault),
     "tab delimited length correction file for NSAF")
    ("sim", po::value<bool>()->default_value(false),
     "if true, scramble the input data for simulated HGScores")
    ("seed", po::value<unsigned int>()->default_value(tim),
     "random seed")
    ("verbose", po::value<bool>()->default_value(false),
     "if true, print status updates")
    ("hgdump", po::value<bool>()->default_value(false),
     "if true, print full contingency tables for all non-zero interactions (LOTS)")
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
  bool verbose = opts["verbose"].as<bool>();
  bool dumpContingency = opts["hgdump"].as<bool>();
  srand( opts["seed"].as<unsigned int>() );

  
  BaitByPrey sbp_mat = BaitByPrey(inFile, lengthCorrectionFile);
  sbp_mat.verbose = verbose;
  
  sbp_mat.reScale(BaitByPrey::NSAF);
  sbp_mat.reScale(BaitByPrey::Sqrt); //commented this out for testing the HGScore without sqrt
  sbp_mat.reScale(BaitByPrey::Round);
  //debugPrint(sbp_mat, "244369", "FBgn0259139");
  //debugPrint(sbp_mat, "244369", "FBgn0022023");
  //debugPrint2(sbp_mat, "/home/glocke/DPiM/augRemap/nrBait_yRej_nBrand_11-16-2016/qdir/test1-12-19-2016/FBgn0022023-FBgn0034237.coappear.sids", "FBgn0022023", "FBgn0034237");
  //exit(0);
  
  //// SIMULATED HYPER-SCORE ////
  if (simulateFlag)
    sbp_mat.simulateDataset();

  PreyByPrey pbp_mat = PreyByPrey(sbp_mat);
  pbp_mat.dumpContingency = dumpContingency;
  
  if (verbose) {
    cout << sbp_mat.sumSum() << " sum(sum( searchByPrey ) )" << endl;
    cout << sbp_mat.maxMax() << " max(max( searchByPrey ) )" << endl;
    cout << sbp_mat.sumSum() << " sum(sum( searchByPrey ) )" << endl;
    cout << sbp_mat.maxMax() << " max(max( searchByPrey ) )" << endl;
    cout << pbp_mat.sumSum() << " sum(sum( preyByPrey   ) )" << endl;
  }

  
  //pbp_mat.print(); DEBUG
  //return 0;

  pbp_mat.computeHyge();
  if (dumpContingency) {
    return EXIT_SUCCESS;
  }
  pbp_mat.json(outFile);
  pbp_mat.outputScores();
  
  return (EXIT_SUCCESS);
}


void debugPrint(BaitByPrey sbp_mat, string sid, string prey)  {   // DEBUG 2016-12-19
  //for(auto const& it : searchId2Idx) {
  int sidIdx = sbp_mat.searchId2Idx.at(sid);
  int preyIdx = sbp_mat.preyRef2Idx.at(prey);
  
  cerr << "sbp[" << sid << "][" << prey << "] = "
       << sbp_mat.searchByPrey.at(sidIdx).at(preyIdx) << endl;    
}

void debugPrint2(BaitByPrey sbp_mat, string sidFile, string prey1, string prey2)  {   // DEBUG 2016-12-19

  std::ifstream is(sidFile);
  std::istream_iterator<string> start(is), end;
  std::vector<string> sids(start, end);
  std::cout << "Read " << sids.size() << " sids" << std::endl;
  
  int prey1Idx = sbp_mat.preyRef2Idx.at(prey1);
  int prey2Idx = sbp_mat.preyRef2Idx.at(prey2);

  int sum=0;
  
  for (std::vector<string>::iterator it = sids.begin(); it != sids.end(); ++it){
    int sidIdx = sbp_mat.searchId2Idx.at(*it);

    int tsc1 = sbp_mat.searchByPrey.at(sidIdx).at(prey1Idx);
    int tsc2 = sbp_mat.searchByPrey.at(sidIdx).at(prey2Idx);
    
    cerr << "sbp[" << *it << "][" << prey1 << "] = " << tsc1 << endl;
    cerr << "sbp[" << *it << "][" << prey2 << "] = " << tsc2 << endl;
    sum+=(tsc1<tsc2)?tsc1:tsc2;
  }
  cerr << "sum = " << sum << endl;
}
