/*
 *  find and report a few statistics for exploratory purposes
 */

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
  size_t row,col;
  double val;
}
#ifndef STAT
#define STAT 1
struct stat {
  size_t row,col;
  double score, diff, var, snr;
  double precScore, prevDiff, prevVar, prevSnr;
};
#endif
#ifndef STAT_VECTOR
#define STAT_VECTOR 1
typedef std::vector<stat> statVector;
#endif
#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector<double> doubleVector;
#endif

namespace po = boost::program_options;
namespace blas = boost::numeric::ublas;

po::variables_map getOptions(int argc, char* argv[]);
stringVector readList(const std::string inFile);
scoreMat getMean(const Network& n);
scoreMat getVar(Network& n);
matElem findMax(const symMat& mat);
void findMax(double& max, size_t& row, size_t& col, const symMat & mat);
void findMax(double& max, size_t& row, size_t& col, const symMat& mat,
	     const symMat& refMat, float refMin);
void pushStat(statVector &s, double score, double diff, double var, double snr,
	      double prevScore, double prevDiff, double prevVar,
	      double prevSnr);
void reportStat(const std::string &outFile, unsigned firstStep,
		const stringIBimap& allProteins, const statVector& s)

/*
void report(const std::string &outFile, unsigned firstStep,
	    const stringIBimap& allProteins,
	    const matElemVector& maxVar, const doubleVector& scoreAtMaxVar,
	    const doubleVector& diffAtMaxVar,
	    const matElemVector& maxDiff, const doubleVector& scoreAtMaxDiff,
	    const doubleVector& varAtMaxDiff,
	    const matElemVector& maxSNR, const doubleVector& scoreAtMaxSNR,
	    const doubleVector& diffAtMaxSNR, const doubleVector& varAtMaxSNR);
*/
template<typename T> void dumpSize(T x);

int main  (int argc, char* argv[]) {
  
  po::variables_map opts = getOptions(argc, argv);

  const std::string jsonList = opts["in"].as<std::string>();
  const std::string outFile = opts["out"].as<std::string>();
  //const unsigned lookback = opts["lookback"].as<unsigned>();
  //const double convCutoff = opts["converge"].as<double>();
  unsigned nProt = opts["nprot"].as<unsigned>();
  
  stringVector jsonFiles = readList(jsonList);
  
  Network n(nProt);

  // var := variance of scores for each pair
  // diff := (meanAtN - meanAtN-1)^2
  // SNR := diff / var
  // where all operations are element-wise (not proper matrix operations)
  statVector maxVar, maxDiff, maxSNR;
  unsigned minCnt = 2;
  {
    unsigned cnt = 0;
    double minScore = 10;
    symMat prevScore, diff, prevDiff, SNR;
    BOOST_FOREACH(std::string jsonFile, jsonFiles) {
      cnt++;
      std::cout << "reading " << jsonFile << std::endl;
      n.addJson(jsonFile);
      
      if (cnt > 1) {
	std::cout << "\tcalc diff" << std::endl;
	diff = prevScore - n.stats.newMean;
	diff = blas::element_prod(diff, diff);
      }
      if (cnt > minCnt) {
	std::cout << "\tobserve statistics" << std::endl;
	n.stats.calcVariance();

	// find statistics at the pair with highest square difference
	maxDiff.push_back(findMax(diff));
	
	scoreAtMaxDiff.push_back(n.stats.newMean(maxDiff.back().row,
						 maxDiff.back().col));
	varAtMaxDiff.push_back(n.stats.var(maxDiff.back().row,
					   maxDiff.back().col));
	

	// find statistics at the pair with highest variance
	max = findMax(n.stats.var);
	maxVar.push_back(findMax(n.stats.var));
	scoreAtMaxVar.push_back(n.stats.newMean(maxVar.back().row,
						maxVar.back().col));
	diffAtMaxVar.push_back(diff(maxVar.back().row, maxVar.back().col));
	


	// find statistics at the pair with highest diff scaled by variance
	SNR = element_div(diff, n.stats.var);
	maxSNR.push_back(findMax(SNR, n.stats.newMean, minScore));
	scoreAtMaxSNR.push_back(n.stats.newMean(maxSNR.back().row,
						maxSNR.back().col));
	varAtMaxSNR.push_back(n.stats.var(maxSNR.back().row,
					  maxSNR.back().col));
	diffAtMaxSNR.push_back(diff(maxSNR.back().row, maxSNR.back().col));
      }
      prevScore = n.stats.newMean;
      prevDiff = diff;
    }
  }

  /*
  dumpSize(n.allProteins);
  dumpSize(maxVar);
  dumpSize(scoreAtMaxVar);
  dumpSize(maxDiff);
  dumpSize(scoreAtMaxDiff);
  dumpSize(maxDiffPrev);
  */
  report(outFile, minCnt, n.allProteins, maxVar, scoreAtMaxVar, diffAtMaxVar,
	 maxDiff, scoreAtMaxDiff, varAtMaxDiff, maxSNR, scoreAtMaxSNR,
	 diffAtMaxSNR, varAtMaxSNR);
	 
	 
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

scoreMat getVar(Network& n) {
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
      ret[row][col] = n.var(row, col);
    }
  }
  return ret;
}

// ignore diagonal elements
void findMax(double& max, size_t& row, size_t& col, const symMat & mat) {
  max=0;
  for (symMat::const_iterator1 it1 = mat.begin1(); it1 !=mat.end1(); ++it1) {
    for (symMat::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
      if (*it2 > max) {
	max = *it2;
	row = it2.index1();
	col = it2.index2();
      }
    }
  }

  return;
}

void findMax(double& max, size_t& row, size_t& col, const symMat& mat,
	     const symMat& refMat, float refMin)
{
  max=0;
  for (symMat::const_iterator1 it1 = mat.begin1(); it1 !=mat.end1(); ++it1) {
    for (symMat::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
      if (refMat(it2.index1(), it2.index2()) < refMin) {
	continue;
      }
      if (*it2 > max) {
	max = *it2;
	row = it2.index1();
	col = it2.index2();
      }
    }
  }

  return;
}


void report(const std::string outFile, unsigned firstStep,
	    const stringIBimap& allProteins,
	    const matElemVector& maxVar, const doubleVector& scoreAtMaxVar,
	    const doubleVector& diffAtMaxVar,
	    const matElemVector& maxDiff, const doubleVector& scoreAtMaxDiff,
	    const doubleVector& varAtMaxDiff,
	    const matElemVector& maxSNR, const doubleVector& scoreAtMaxSNR,
	    const doubleVector& diffAtMaxSNR, const doubleVector& varAtMaxSNR)
{
  std::ofstream out;
  out.open(outFile);
  if (out.fail()) {
    std::cerr << "failed to write to " << outFile << std::endl;
    return;
  }

  out << "step\tmaxVar\tMVscore\tMVdiff\tMV1\tMV2\t"
      << "maxDiff\tMDscore\tMDvar\tMD1\tMD2\t"
      << "maxSNR\tSNRscore\tSNRdiff\tSNRvar\tSNR1\tSNR2\t"
	    << std::endl;
  for (size_t i = 0; i < maxVar.size(); ++i) {
    out << firstStep+i << "\t";
    out << maxVar.at(i).val << "\t";
    out << scoreAtMaxVar.at(i) << "\t";
    out << diffAtMaxVar.at(i) << "\t";
    out << allProteins.right.at(maxVar.at(i).row) << "\t";
    out << allProteins.right.at(maxVar.at(i).col) << "\t";
    
    out << maxDiff.at(i).val << "\t";
    out << scoreAtMaxDiff.at(i) << "\t";
    out << varAtMaxDiff.at(i) << "\t";
    out << allProteins.right.at(maxDiff.at(i).row) << "\t";
    out << allProteins.right.at(maxDiff.at(i).col) << "\t";

    out << maxSNR.at(i).val << "\t";
    out << scoreAtMaxSNR.at(i) << "\t";
    out << diffAtMaxSNR.at(i) << "\t";
    out << varAtMaxSNR.at(i) << "\t";
    out << allProteins.right.at(maxSNR.at(i).row) << "\t";
    out << allProteins.right.at(maxSNR.at(i).col) << "\t";
    out << std::endl;
  }

  out.close();

  return;
}

template<typename T> void dumpSize(T x) {
  std::cout << x.size() << std::endl;
}
