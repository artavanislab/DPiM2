#include "Network.hpp"
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector< std::string > stringVector;
#endif
#ifndef SCORE_ROW
#define SCORE_ROW 1
typedef std::map<std::string, double> scoreRow;
#endif
#ifndef SCORE_MAT
#define SCORE_MAT 1
typedef std::map< std::string, scoreRow> scoreMat;
#endif
#ifndef PROT_PAIR
#define PROT_PAIR 1
typedef std::pair<std::string, std::string> protPair;
#endif

namespace po = boost::program_options;

po::variables_map getOptions(int argc, char* argv[]);
stringVector readList(const std::string inFile);
scoreMat getMean(const Network& n);
scoreMat getSD(const Network& n);
protPair findMax(const scoreMat& mat);

int main  (int argc, char* argv[]) {
  
  po::variables_map opts = getOptions(argc, argv);

  const std::string jsonList = opts["in"].as<std::string>();
  const std::string outFile = opts["out"].as<std::string>();
  const unsigned lookback = opts["lookback"].as<unsigned>();
  const double convCutoff = opts["converge"].as<double>();
  
  stringVector oldJson = readList(jsonList);
  if (oldJson.size() < lookback) {
    std::cerr << "json list has fewer elements than lookback" << std::endl;
    exit(1);
  }
  stringVector newJson;
  while(newJson.size() < lookback) {
    newJson.push_back(oldJson.back());
    oldJson.pop_back();
  }

  std::cout << "Ingesting JSON..." << std::endl;
  Network n(oldJson);

  std::cout << "Gather statistics..." << std::endl;;
  scoreMat oldMean = getMean(n);

  return 0;
  BOOST_FOREACH(std::string jsonFile, newJson) {
    n.addJson(jsonFile);
  }

  std::cout << "get first Mean" << std::endl;;
  scoreMat newMean = getMean(n);
  scoreMat sd = getMean(n);

  protPair maxPair = findMax(sd);
  double maxSD = sd.at(maxPair.first).at(maxPair.second);
  double diff = newMean.at(maxPair.first).at(maxPair.second) -
    oldMean.at(maxPair.first).at(maxPair.second);
  double testStat = std::abs(diff)/maxSD;
  if (convCutoff > testStat) {
    std::cout << "success" << std::endl;
  } else {
    std::cout << "test statistic for " << maxPair.first << "<->"
	      << maxPair.second << " = " << testStat << "which exceeds cutoff"
	      << std::endl;
  }
}

po::variables_map getOptions(int argc, char* argv[]) {
  std::string jsonListDefault = "/home/glocke/DPiM/dpim4/consTest0_01-11-2016/json.list";
  po::options_description desc("Usage");
  desc.add_options()
    ("out", po::value<std::string>()->required(),
     "REQUIRED: write stuff here (JSON)")
    ("in", po::value<std::string>()->default_value(jsonListDefault),
     "(soon to be) REQUIRED: json.list")
    ("lookback", po::value<unsigned>()->default_value(1),
     "compare the latest iteration to the one ~ ago")
    ("converge", po::value<double>()->default_value(0.001),
     "largest allowable relative change")
    ;
  po::variables_map opts;
  po::store(po::parse_command_line(argc, argv, desc), opts);
  try {
    po::notify(opts);
  } catch (std::exception& e) {
    std::cerr << desc << std::endl;
    std::cerr << "Error: " << e.what() << std::endl;
    exit(1);
  }
  return opts;
}

stringVector readList(const std::string inFile) {

  stringVector ret;
  
  std::string line;
  std::ifstream myFile;
  myFile.open(inFile.c_str());
  
  if (!myFile.is_open()) {
    std::cerr << "failed to open list file " << inFile << std::endl;
    return ret;
  }
  
  while (! myFile.eof() ) {
    getline (myFile,line);
    if (line[0]=='#') { continue; } // comments
    if (line.length()==0) { continue; } // eof

    //boost::algorithm::trim_right(line); chomp apparently unnecessary
    ret.push_back(line);
  }

  return ret;
}

// get mean for every protein pair
scoreMat getMean(const Network& n) {
  scoreMat ret;

  stringSet::const_iterator p1, p2;
  for (p1=n.allProteins.begin(); p1 != n.allProteins.end(); ++p1) {
    p2 = p1;
    ++p2;
    for (; p2 != n.allProteins.end(); ++p2) {
      ret[*p1][*p2] = n.mean(*p1, *p2);
      std::cout << "ret[" << *p1 << "][" << *p2 << "] = " << ret[*p1][*p2]
		<< std::endl;
    }
  }

  return ret;
}

scoreMat getSD(const Network& n) {
  scoreMat ret;
  return ret;
}
protPair findMax(const scoreMat& mat) {

  protPair ret;
  return ret;
}

