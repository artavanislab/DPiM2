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
#ifndef MAT_ELEM
#define MAT_ELEM 1
struct matElem {
  size_t i,j;
  double val;
};
#endif
#ifndef MAT_ELEM_VECTOR
#define MAT_ELEM_VECTOR 1
typedef std::vector<matElem> matElemVector;
#endif
#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector<double> doubleVector;
#endif

namespace po = boost::program_options;

po::variables_map getOptions(int argc, char* argv[]);
stringVector readList(const std::string inFile);
scoreMat getMean(const Network& n);
scoreMat getVar(const Network& n);
matElem findMax(const symMat & mat);

void report(const std::string outFile, const stringIBimap& allProteins,
	    const matElemVector& maxVar, const doubleVector& scoreAtMaxVar,
	    const matElemVector& maxDiff, const doubleVector& scoreAtMaxDiff); 

int main  (int argc, char* argv[]) {
  
  po::variables_map opts = getOptions(argc, argv);

  const std::string jsonList = opts["in"].as<std::string>();
  const std::string outFile = opts["out"].as<std::string>();
  const unsigned lookback = opts["lookback"].as<unsigned>();
  const double convCutoff = opts["converge"].as<double>();
  unsigned nProt = opts["nprot"].as<unsigned>();
  
  stringVector jsonFiles = readList(jsonList);
  
  Network n(nProt);

  std::vector<matElem> maxVar, maxDiff;
  std::vector<double> scoreAtMaxVar, scoreAtMaxDiff;
  unsigned minCnt = 2;
  {
    unsigned cnt = 0;
    symMat prevScore, diff;
    BOOST_FOREACH(std::string jsonFile, jsonFiles) {
      cnt++;
      std::cout << "reading " << jsonFile << std::endl;
      n.addJson(jsonFile);
      
      if (cnt > minCnt) {
	n.stats.calcVariance();
	maxVar.push_back(findMax(n.stats.var));
	scoreAtMaxVar.push_back(n.stats.newMean(maxVar.back().i,
						maxVar.back().j));
	
	diff = prevScore - n.stats.newMean;
	maxDiff.push_back(findMax(blas::element_prod(diff, diff)));
	scoreAtMaxVar.push_back(n.stats.newMean(maxVar.back().i,
						maxVar.back().j));
      }
      prevScore = n.stats.newMean;
    }
  }

  report(outFile, n.allProteins, maxVar, scoreAtMaxVar, maxDiff,
	 scoreAtMaxDiff);
	 
	 
  /*
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
  */
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
    ("nprot", po::value<unsigned>()->default_value(7016),
     "number of unique proteins expected in dataset")
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

  std::string row, col;
  stringIBimap::const_iterator p1, p2;
  for (p1=n.allProteins.begin(); p1 != n.allProteins.end(); ++p1) {
    p2 = p1;
    ++p2;
    for (; p2 != n.allProteins.end(); ++p2) {
      row = p1->left;
      col = p2->left;
      if (row.compare(col) > 0) {
	std::swap(row, col);
      }
      
      ret[row][col] = n.mean(row, col);
    }
  }

  return ret;
}

scoreMat getVar(const Network& n) {
  scoreMat ret;

  std::string row, col;
  stringIBimap::const_iterator p1, p2;
  for (p1=n.allProteins.begin(); p1 != n.allProteins.end(); ++p1) {
    p2 = p1;
    ++p2;
    for (; p2 != n.allProteins.end(); ++p2) {
      row = p1->left;
      col = p2->left;
      if (row.compare(col) > 0) {
	std::swap(row, col);
      }
      ret[row][col] = n.variance(row, col);
    }
  }
  return ret;
}

// ignore diagonal elements
matElem findMax(const symMat & mat) {
  size_t N = mat.size1();
  double max=0;
  matElem ret;
  for (size_t i = 1; i < N-1; ++i) {
    for (size_t j = i+1; j < N; ++j) {
      if (mat(i,j) > max) {
	max=mat(i,j);
	ret.i=i;
	ret.j=j;
	ret.val=max;
      }
    }
  }

  return val;
}
/*
protPair findMax(const scoreMat& mat) {

  protPair ret;
  return ret;
}

*/

void report(const std::string outFile, const stringIBimap& allProteins,
	    const matElemVector& maxVar, const doubleVector& scoreAtMaxVar,
	    const matElemVector& maxDiff, const doubleVector& scoreAtMaxDiff) {
  
}
