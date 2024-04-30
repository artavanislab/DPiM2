/*
 *  find and report a few statistics for exploratory purposes
 */

#include "Network.hpp"
#include <boost/algorithm/string.hpp>
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
#ifndef PROT_PAIR
#define PROT_PAIR 1
typedef std::pair<std::string, std::string> protPair;
#endif
#ifndef STAT
#define STAT 1
struct stat {
  size_t row,col;
  double score, diff, var;
  double prevScore, prevDiff, prevVar;
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

/*
#ifndef SCORE_ROW
#define SCORE_ROW 1
typedef std::map<std::string, double> scoreRow;
#endif
#ifndef SCORE_MAT
#define SCORE_MAT 1
typedef std::map< std::string, scoreRow> scoreMat;
#endif

 */

namespace po = boost::program_options;
namespace blas = boost::numeric::ublas;

po::variables_map getOptions(int argc, char* argv[]);
stringVector readList(const std::string inFile);
void mainParse(const stringVector& jsonFiles, Network& n, unsigned minCnt,
	       statVector& maxVar, statVector& maxDiff, statVector& maxSNR,
	       statVector& foldChange);
//	       statVector& diffRat, statVector& SNRRat);
void findMax(double& max, size_t& row, size_t& col, const symMat & mat);
void findMax(double& max, size_t& row, size_t& col, const symMat& mat,
	     const symMat& refMat, float refMin);
void pushStat(statVector &sv, size_t row, size_t col, double score, double diff,
	      double var, double prevScore, double prevDiff, double prevVar);
std::string makeFileName(std::string base, std::string insert);
void reportStat(const std::string &outFile, unsigned firstStep,
		const stringIBimap& allProteins, const statVector& sv);
void reportEdges(const std::string &outFile, const std::string &colName,
		 const stringIBimap& allProteins, const symMat& mat,
		 double min);

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

  // find the maximum edge for each statistic
  // var := variance of scores for each pair
  // diff := (meanAtN - meanAtN-1)^2
  // SNR := diff / var
  // where all operations are element-wise (not proper matrix operations)
  unsigned minCnt = 2;
  statVector maxVar, maxDiff, maxSNR, foldChange;
  mainParse(jsonFiles, n, minCnt, maxVar, maxDiff, maxSNR, foldChange);
  //statVector maxVar, maxDiff, maxSNR, diffRat, SNRRat;
  //mainParse(jsonFiles, n, minCnt, maxVar, maxDiff, maxSNR, diffRat, SNRRat);

  /*
  dumpSize(n.allProteins);
  dumpSize(maxVar);
  dumpSize(maxDiff);
  dumpSize(maxDiffPrev);
  */
  
  reportStat(makeFileName(outFile,"maxVar"), minCnt+1, n.allProteins, maxVar);
  reportStat(makeFileName(outFile,"maxDiff"), minCnt+1, n.allProteins, maxDiff);
  reportStat(makeFileName(outFile,"maxSNR"), minCnt+1, n.allProteins, maxSNR);
  reportStat(makeFileName(outFile,"foldChange"), minCnt+1, n.allProteins,
	     foldChange);
  //reportStat(makeFileName(outFile,"diffRat"), minCnt+1, n.allProteins, diffRat);
  //reportStat(makeFileName(outFile,"SNRRat"), minCnt+1, n.allProteins, SNRRat);

  double minReport = 2.5;
  reportEdges(makeFileName(outFile,"aveScores"), "aveScore", n.allProteins,
	      n.stats.newMean, minReport);
  reportEdges(makeFileName(outFile,"scoreVariance"), "scoreVar", n.allProteins,
	      n.stats.variance(), minReport);
  
	 
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
    ("in", po::value<std::string>()->required(),
     "REQUIRED: json.list")
    ("out", po::value<std::string>()->required(),
     "REQUIRED: write stuff here (JSON)")
    ("lookback", po::value<unsigned>()->default_value(1),
     "compare the latest iteration to the one ~ ago")
    ("converge", po::value<double>()->default_value(0.001),
     "largest allowable relative change")
    ("nprot", po::value<unsigned>()->default_value(7210),
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

void mainParse(const stringVector& jsonFiles, Network& n, unsigned minCnt,
	       statVector& maxVar, statVector& maxDiff, statVector& maxSNR,
	       statVector& foldChange)
//statVector& diffRat, statVector& SNRRat)
{
  unsigned cnt = 0;
  double minScore = 10;
  symMat prevScore, prevVar, diff, prevDiff, SNR;
  double max;
  size_t row, col;
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
      SNR = element_div(diff, n.stats.var);

      // find statistics at the pair with highest square difference
      findMax(max, row, col, diff);
      pushStat(maxDiff, row, col, n.stats.newMean(row, col), diff(row, col),
	       n.stats.var(row, col), prevScore(row, col), prevDiff(row, col),
	       prevVar(row, col));
	
      // find statistics at the pair with highest variance
      findMax(max, row, col, n.stats.var);
      pushStat(maxVar, row, col, n.stats.newMean(row, col), diff(row, col),
	       n.stats.var(row, col), prevScore(row, col), prevDiff(row, col),
	       prevVar(row, col));

      // find statistics at the pair with highest diff scaled by variance
      findMax(max, row, col, SNR, n.stats.newMean,
	      minScore);
      pushStat(maxSNR, row, col, n.stats.newMean(row, col), diff(row, col),
	       n.stats.var(row, col), prevScore(row, col), prevDiff(row, col),
	       prevVar(row, col));


      // find statistics at the pair with highest diff_n / score_n 
      findMax(max, row, col, element_div(diff, n.stats.newMean),
	      n.stats.newMean, minScore);
      pushStat(foldChange, row, col, n.stats.newMean(row, col), diff(row, col),
	       n.stats.var(row, col), prevScore(row, col), prevDiff(row, col),
	       prevVar(row, col));

      /*
      // find statistics at the pair with highest diff_n / diff_n-1 
      findMax(max, row, col, element_div(diff, prevDiff), n.stats.newMean,
	      minScore);
      pushStat(diffRat, row, col, n.stats.newMean(row, col), diff(row, col),
	       n.stats.var(row, col), prevScore(row, col), prevDiff(row, col),
	       prevVar(row, col));

      // find statistics at the pair with highest snr_n / snr_n-1 
      findMax(max, row, col, element_div(SNR, element_div(prevDiff, prevVar)),
	      n.stats.newMean, minScore);
      pushStat(SNRRat, row, col, n.stats.newMean(row, col), diff(row, col),
	       n.stats.var(row, col), prevScore(row, col), prevDiff(row, col),
	       prevVar(row, col));
      */
    }

    prevScore = n.stats.newMean;
    if (cnt > 1) {
      prevDiff = diff;
      n.stats.calcVariance();
      prevVar = n.stats.var;
    }
  }

  return;
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

// the the max value of mat subject to the requirement that the corresponding
// element (same row/col) in refMat is at least refMin
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

void pushStat(statVector &sv, size_t row, size_t col, double score, double diff,
	      double var, double prevScore, double prevDiff, double prevVar)
	      
{
  /*
  size_t row,col;
  double score, diff, var;
  double prevScore, prevDiff, prevVar;
  */
  stat s;

  s.row = row;
  s.col = col;
  s.score = score;
  s.diff = diff;
  s.var = var;
  s.prevScore = prevScore;
  s.prevDiff = prevDiff;
  s.prevVar = prevVar;

  sv.push_back(s);
}

// combine "dir/file.name.ext" and "ins" to make "dir/file.name.ins.ext"
std::string makeFileName(std::string base, std::string insert) {
  stringVector spl;
  boost::split(spl, base, boost::is_any_of("."));
  if (spl.size() == 1) {
    return base+"."+insert;
  }
  std::string last = spl.back();
  spl.pop_back();
  spl.push_back(insert);
  spl.push_back(last);
  return boost::join(spl, ".");
}

void reportStat(const std::string &outFile, unsigned firstStep,
		const stringIBimap& allProteins, const statVector& sv)
{
  std::ofstream out;
  out.open(outFile);
  if (out.fail()) {
    std::cerr << "failed to write to " << outFile << std::endl;
    return;
  }

  out << "step\tprot1\tprot2\tscore\tdiff\tvar\t"
      << "prevScore\tprevDiff\tprevVar" << std::endl;

  unsigned step = firstStep;
  BOOST_FOREACH(stat s, sv) {
    out << step << "\t";
    out << allProteins.right.at(s.row) << "\t";
    out << allProteins.right.at(s.col) << "\t";
    out << s.score << "\t";
    out << s.diff << "\t";
    out << s.var << "\t";
    out << s.prevScore << "\t";
    out << s.prevDiff << "\t";
    out << s.prevVar << std::endl;
    
    step++;
  }
}


void reportEdges(const std::string &outFile, const std::string &colName,
		 const stringIBimap& allProteins, const symMat& mat, double min)
{
  std::ofstream out;
  out.open(outFile);
  if (out.fail()) {
    std::cerr << "failed to write to " << outFile << std::endl;
    return;
  }

  std::string node1, node2;
  out << "node1\tnode2\t" << colName << std::endl;
  for (symMat::const_iterator1 it1 = mat.begin1(); it1 !=mat.end1(); ++it1) {
    for (symMat::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
      if (*it2 < min) {
	continue;
      }
      node1 = allProteins.right.at(it2.index1());
      node2 = allProteins.right.at(it2.index2());

      out << node1 << "\t" << node2 << "\t" << *it2 << std::endl;
    }
  }

  return;
}

template<typename T> void dumpSize(T x) {
  std::cout << x.size() << std::endl;
}
