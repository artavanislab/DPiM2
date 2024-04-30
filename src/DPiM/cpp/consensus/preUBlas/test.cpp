#include "Network.hpp"
#include "rapidjson/writer.h" 
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdio>
#include <vector>
#include <random>
#include <cmath>
#include <map>

//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector< double > doubleVector; 
#endif
#ifndef DOUBLE_VECTOR_VECTOR
#define DOUBLE_VECTOR_VECTOR 1
typedef std::vector< doubleVector > doubleVectorVector; 
#endif
#ifndef SCORE_TYPE
#define SCORE_TYPE 1
typedef std::map< std::string, std::map< std::string, doubleVector > > scoreType; 
#endif

namespace rj = rapidjson;

stringVector randomSubset(stringVector list, size_t sizeOfSubset);
void fillRandMat(doubleVectorVector &randMat, scoreType &allScores,
		 size_t scoreIndex, stringVector thisProtList, 
		 double interactionProb, std::mt19937 &mt );
void writeJson(const std::string &outFile,
	       const doubleVectorVector& preyByPreyHyge,
	       const stringVector& protList);
std::string createJson(const doubleVectorVector& preyByPreyHyge,
		       const stringVector& protList);
void simpleJsonTest(const char json[]);
std::string slurp(const std::string &inFile);

int main  (int argc, char* argv[]) {
  
  size_t nProt = 100;
  stringVector protList;
  for (size_t i=0; i < nProt; ++i) {
    char name[20];
    sprintf(name, "protein%04lu", i);
    std::string n(name);
    protList.push_back(n);
  }

  size_t nTest = 100;
  
  scoreType allScores; 
  BOOST_FOREACH(std::string p1, protList) {
    BOOST_FOREACH(std::string p2, protList) {
      allScores[p1][p2].resize(nTest, 0);
    }
  }
  
  std::random_device rd;
  std::mt19937 mt(rd());
  std::srand(rd());
  
  stringVector files;
  double interactionProb = 0.2; 
  for (size_t i=0; i < nTest; ++i) {
    stringVector thisProtList = randomSubset(protList, protList.size() * 0.8);
    doubleVectorVector randMat;
    fillRandMat(randMat, allScores, i, thisProtList, interactionProb, mt);

    std::string json = createJson(randMat, thisProtList);
    simpleJsonTest(json.c_str());

    std::string outFile = std::tmpnam(nullptr);
    writeJson(outFile, randMat, thisProtList);
    files.push_back(outFile);
    
    std::string s = slurp(outFile);
    simpleJsonTest(s.c_str());
    std::cout << outFile << std::endl;
  }
  
  std::cout << "making Network...\n";
  Network n(files);
  std::cout << "\t...made it\n";

  bool passFlag = true;
  double thresh = 0.000001;
  BOOST_FOREACH(std::string p1, protList) {
    BOOST_FOREACH(std::string p2, protList) {
      if (! p1.compare(p2)) {
	continue;
      }
      doubleVector v = allScores.at(p1).at(p2);
      std::cout << p1 << "\t" << p2 << "\t" << std::endl; 
      double sum = std::accumulate(v.begin(), v.end(), 0.0);
      double mean = sum / v.size();

      double sumSq = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
      double sd = std::sqrt(sumSq / v.size() - mean * mean);
      sd *= std::sqrt(v.size() / (v.size() - 1.0));

      double nMean = n.mean(p1, p2);
      double nSd = n.sd(p1, p2);

      std::cout << std::setprecision(3) << "\tmean - " << mean << " vs "
		<< nMean << std::endl;
      std::cout << std::setprecision(3) << "\tsd - " << sd << " vs "
		<< n.sd(p1, p2) << std::endl;
      if (std::abs(1 - nMean/mean) > thresh || std::abs(1 - nSd/sd) > thresh) {
	passFlag = false;
	break;
      }
    }
    if (! passFlag) {
      break;
    }
  }

  if (passFlag) {
    std::cout << "All tests passed" << std::endl;
    return 0;
  } else {
    std::cerr << "Failed test" << std::endl;
    return 1;
  }
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
void fillRandMat(doubleVectorVector &randMat, scoreType &allScores,
		 size_t scoreIndex, stringVector thisProtList, 
		 double interactionProb, std::mt19937 &mt)
{
  size_t nProt = thisProtList.size();
  doubleVector blank(nProt, 0);
  randMat.resize(nProt, blank);

  double maxScore = 10;
  std::bernoulli_distribution prob(interactionProb);
  std::uniform_real_distribution<double> score(1.0, maxScore);

  for (size_t i=0; i < nProt-1; ++i) {
    std::string p1 = thisProtList.at(i);
    for (size_t j=i+1; j < nProt; ++j) {
      if (prob(mt)) {
	std::string p2 = thisProtList.at(j);
	double s = score(mt);
	randMat.at(j).at(i) = randMat.at(i).at(j) = s;
	allScores.at(p1).at(p2).at(scoreIndex) =
	  allScores.at(p2).at(p1).at(scoreIndex) = s;
	/*
	if (!(p1.compare("protein0000") || p2.compare("protein0001")) ||
	    !(p1.compare("protein0001") || p2.compare("protein0000")) )
	{
	  std::cout << "allScores.at(" << p1 << ").at(" << p2 << ") = " 
		    << s << std::endl;
	}
	DEBUG
	*/
      }
    }
  }

  return;
}

void writeJson(const std::string &outFile,
	       const doubleVectorVector& preyByPreyHyge,
	       const stringVector& protList)
{

  std::ofstream out;
  out.open(outFile);
  if (out.fail()) {
    std::cerr << "failed to write to " << outFile << std::endl;
  }
  
  out << createJson(preyByPreyHyge, protList) << std::endl;
  out.close();
}

std::string createJson(const doubleVectorVector& preyByPreyHyge,
		       const stringVector& protList)
{
  rj::StringBuffer sb;
  rj::Writer<rj::StringBuffer> writer(sb);
  //rj::PrettyWriter approximately doubles the disk space
  writer.StartObject();
  // report HGScore -- only the non-zero elements
  const double minScore = 0.00001; // report only scores higher than ~
  writer.String("score");
  writer.StartObject();
  for (size_t i = 0; i < preyByPreyHyge.size()-1; ++i) {
    char iString[10];
    sprintf(iString, "%lu", i);
    writer.String(iString);
    writer.StartObject();
    for (size_t j = (i + 1); j < preyByPreyHyge[0].size(); ++j) {
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
  for (size_t i = 0; i < protList.size(); ++i) {
    char gString[10]; // I'm so clever
    sprintf(gString, "%lu", i);
    writer.String(gString);
    writer.String(protList.at(i).c_str());
  }
  writer.EndObject();
  writer.EndObject();
  return sb.GetString();
}

void simpleJsonTest(const char json[]) {
  rj::Document jsonObj;
  jsonObj.Parse(json);
  //std::cout << json << std::endl;
  
  assert(jsonObj.IsObject());
  /*
  assert(jsonObj.HasMember("score"));
  assert(jsonObj["score"].IsObject());
  assert(jsonObj.HasMember("ids"));
  assert(jsonObj["ids"].IsObject());
  */
}

std::string slurp(const std::string &inFile) {
  std::ifstream in(inFile, std::ios::in | std::ios::binary);
  if (in) {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return(contents);
  }
  throw(errno);
}
