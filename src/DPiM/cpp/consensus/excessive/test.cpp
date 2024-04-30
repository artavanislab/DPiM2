#include "Network.hpp"
#include "rapidjson/writer.h" 
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <vector>
#include <random>
#include <map>
#include <boost/foreach.hpp>

//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

#ifndef DOUBLE_VECTOR_VECTOR
#define DOUBLE_VECTOR_VECTOR 1
typedef std::vector< doubleVector > doubleVectorVector; 
#endif

#ifndef SCORE_TYPE
#define SCORE_TYPE 1
typedef std::map< std::string, std::map< std::string, doubleVector > > scoreType; 
#endif

stringVector randomSubset(stringVector list, double fraction);
void fillRandMat(doubleVectorVector randMat, scoreType targetScores,
		 size_t targetIndex, stringVector thisProtList, 
		 double interactionProb, std::mt19937 mt );

void json(const std::string &outFile, doubleVectorVector preyByPreyHyge,
	  stringVector protList);

int main  (int argc, char* argv[]) {

  size_t nProt = 10;
  stringVector protList;
  for (size_t i=0; i < nProt; ++i) {
    char name[20];
    sprintf(name, "protein%04lu", i);
    std::string n(name);
    protList.push_back(n);
  }

  size_t nTest = 10;
  
  scoreType targetScores; // Network should reproduce this
  BOOST_FOREACH(std::string p1, protList) {
    BOOST_FOREACH(std::string p2, protList) {
      targetScores.at(p1).at(p2).resize(nTest, 0);
    }
  }
  
  std::random_device rd;
  std::mt19937 mt(rd());
  
  stringVector files;
  double interactionProb = 0.2; // fairly dense
  for (size_t i=0; i < nTest; ++i) {
    stringVector thisProtList = randomSubset(protList, protList.size() * 0.8);
    doubleVectorVector randMat;
    fillRandMat(randMat, targetScores, i, thisProtList, interactionProb, mt);
    
    std::string outFile = std::tmpnam(nullptr);
    json(outFile, randMat, thisProtList);
    files.push_back(outFile);
  }
  /*
  po::options_description desc("Usage");
  desc.add_options()
    ("in", po::value<string>()->required(),
     "")
    ("out", po::value<string>()->required(),
     "")
  */
  //Network n(jsonList);
}

// return a random subset of list
stringVector randomSubset(stringVector list, size_t sizeOfSubset) {
  std::random_shuffle(list.begin(), list.end());
  stringVector ret;
  while (ret.size() < sizeOfSubset) {
    ret.push_back(list.back());
    list.pop_back();
  }

  return ret;
}

// create a (Erdos-Renyi) random score network
void fillRandMat(doubleVectorVector randMat, scoreType targetScores,
		 size_t targetIndex, stringVector thisProtList, 
		 double interactionProb, std::mt19937 mt)
{
  size_t nProt = thisProtList.size();
  randMat.resize(nProt);
  BOOST_FOREACH(doubleVector dv, randMat) {
    dv.resize(nProt);
  }

  double maxScore = 10;
  std::bernoulli_distribution prob(interactionProb);
  std::uniform_real_distribution<double> score(1.0, maxScore);

  for (size_t i; i < nProt-1; ++i) {
    std::string p1 = thisProtList.at(i);
    for (size_t j=i+1; j < nProt; ++j) {
      if (prob(mt)) {
	std::string p2 = thisProtList.at(j);
	double s = score(mt);
	randMat.at(j).at(i) = randMat.at(i).at(j) = s;
	targetScores.at(p1).at(p2).at(targetIndex) =
	  targetScores.at(p2).at(p1).at(targetIndex) = s;
      }
      std::cout << "\trandMat.at(" << j << ").at(" << i << ") = " 
		<< randMat.at(i).at(j) << std::endl;
    }
  }

  return;
}

void json(const std::string &outFile, doubleVectorVector preyByPreyHyge,
	  stringVector protList)
{

  std::ofstream out;
  out.open(outFile.c_str());
  if (out.fail()) {
    std::cerr << "failed to write to " << outFile << std::endl;
  }
  
  rapidjson::StringBuffer sb;
  rapidjson::Writer<rapidjson::StringBuffer> writer(sb);
  //rapidjson::PrettyWriter approximately doubles the disk space
  writer.StartObject();

  // report HGScore -- only the non-zero elements
  const double minScore = 0.00001; // report only elements with score higher than
  writer.String("score");
  writer.StartObject();

  for (size_t i = 0; i < preyByPreyHyge.size()-1; i++) {
    char iString[10];
    sprintf(iString, "%lu", i);
    writer.String(iString);
    writer.StartObject();
    for (size_t j = (i + 1); j < preyByPreyHyge[0].size(); j++) {
      double s = preyByPreyHyge[i][j];
      if (s > minScore) {
	char jString[10];
	sprintf(jString, "%lu", j);
	writer.String(jString);
	writer.Double(preyByPreyHyge[i][j]);
      }
    }
    writer.EndObject();
  }
  writer.EndObject();

  // report the rownames
  writer.String("ids");
  writer.StartObject();
  //for ();
  writer.EndObject();
  writer.EndObject();
  
  out << sb.GetString() << std::endl;
  out.close();
}
