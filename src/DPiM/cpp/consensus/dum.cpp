#include "Network.hpp"
#include "preUBlas/RunningStat.hpp"
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
#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector< double > doubleVector;
#endif
#ifndef DOUBLE_VECTOR_VECTOR
#define DOUBLE_VECTORP_VECTOR 1
typedef std::vector< doubleVector > doubleVectorVector; 
#endif
#ifndef STRING_SET
#define STRING_SET 1
typedef std::set< std::string > stringSet;
#endif
#ifndef STRING_RS_MAP
#define STRING2RSMAP 1
typedef std::map< std::string, RunningStat > stringRSMap;
#endif
#ifndef PROT_PROT_MAP
#define PROT_PROT_MAP 1
typedef std::map< std::string, stringRSMap > protProtMap;
#endif

namespace rj = rapidjson;

stringVector readList(const std::string inFile);
std::string slurp(const std::string& inFile);

void addJson(const std::string &f, unsigned nRun, stringSet &allProteins,
	     protProtMap& stats);
stringIBimap fillFbgnMap(rj::Value::ConstMemberIterator begin,
			 rj::Value::ConstMemberIterator end,
			 stringSet &allProteins);
void updateStats(rapidjson::Value::ConstMemberIterator begin,
		 rapidjson::Value::ConstMemberIterator end,
		 const stringIBimap& indexMap, unsigned nRun,
		 protProtMap& stats );

  
int main (int argc, char* argv[]) {
  //std::string jsFile = "/tmp/file149VQ0";
  std::string jsList = "/home/glocke/DPiM/cpp/consensus/dummyData/json.list";
  stringVector sv = readList(jsList);
  while(sv.size() > 3) {
    sv.pop_back();
  }
  
  for (int i=0; i<2; ++i) {
    std::cout << i << std::endl;
    unsigned nRun=0;
    stringSet allProteins;
    protProtMap stats;
    BOOST_FOREACH(const std::string f, sv) {
      ++nRun;

      std::cout << "\treading " << f << "..." << std::endl;
      std::cout << "\t\tslurping" << std::endl;
      std::string json = slurp(f);
      std::cout << "\t\t\tjson.size() = " << json.size() << std::endl;

      // parse
      std::cout << "\t\tparsing" << std::endl;
      rj::Document jsonObj;
      jsonObj.Parse(json.c_str());

      std::cout << "\t\tverify..." << std::endl;
      // verify that the json looks correct
      assert(jsonObj.IsObject());
      assert(jsonObj.HasMember("score"));
      assert(jsonObj["score"].IsObject());
      assert(jsonObj.HasMember("ids"));
      assert(jsonObj["ids"].IsObject());

      std::cout << "\t\tfillFbgnMap...\n";
      stringIBimap fbgnMap = fillFbgnMap(jsonObj["ids"].MemberBegin(),
					 jsonObj["ids"].MemberEnd(),
					 allProteins);

      std::cout << "\t\tupdateStats..\n";
      updateStats(jsonObj["score"].MemberBegin(), jsonObj["score"].MemberEnd(),
		  fbgnMap, nRun, stats);
      
    }
  }

  for (int i=0; i<1; ++i) {
    std::cout << i << std::endl;
    Network n(sv, 0);
  }
}


std::string slurp(const std::string& inFile) {
  std::ifstream in(inFile);
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

void addJson(const std::string &f, unsigned nRun, stringSet &allProteins,
	     protProtMap& stats)
{
  std::cout << "\t\tslurping" << std::endl;
  std::string json = slurp(f);
  std::cout << "\t\t\tjson.size() = " << json.size() << std::endl;

  // parse
  std::cout << "\t\tparsing" << std::endl;
  rj::Document jsonObj;
  jsonObj.Parse(json.c_str());

  std::cout << "\t\tverify..." << std::endl;
  // verify that the json looks correct
  assert(jsonObj.IsObject());
  assert(jsonObj.HasMember("score"));
  assert(jsonObj["score"].IsObject());
  assert(jsonObj.HasMember("ids"));
  assert(jsonObj["ids"].IsObject());

  std::cout << "\t\tfillFbgnMap...\n";
  stringIBimap fbgnMap = fillFbgnMap(jsonObj["ids"].MemberBegin(),
				     jsonObj["ids"].MemberEnd(),
				     allProteins);

  std::cout << "\t\tupdateStats..\n";
  updateStats(jsonObj["score"].MemberBegin(), jsonObj["score"].MemberEnd(),
	      fbgnMap, nRun, stats);
  /*
  */
}

stringIBimap fillFbgnMap(rj::Value::ConstMemberIterator begin,
			 rj::Value::ConstMemberIterator end,
			 stringSet &allProteins)
{
  stringIBimap ret;
  typedef stringIBimap::value_type keyValue;
  std::string fbgn;
  size_t id;
  rj::Value::ConstMemberIterator p1;
  size_t maxId = 0;
  for (p1 = begin; p1 != end; ++p1) {
    id = boost::lexical_cast<size_t>(p1->name.GetString());
    fbgn = p1->value.GetString();
    //std::cout << "(fbgn, id) = ( "<< fbgn << ", " << id << ")\n"; // DEBUG
    ret.insert(keyValue(fbgn, id));
    allProteins.insert(fbgn);
    if (maxId < id) {
      maxId = id;
    }
  }

  // add in any protein we've seen before 
  // this way, it's easy to keep track of proteins that have no edges 
  stringSet::iterator all;
  for (all = allProteins.begin(); all != allProteins.end(); ++all) {
    if (ret.left.count(*all) == 0) {
      ret.insert(keyValue(*all, ++maxId));
    }
  }
  return ret;
}

void updateStats(rj::Value::ConstMemberIterator begin,
		 rj::Value::ConstMemberIterator end,
		 const stringIBimap& fbgnMap, unsigned nRun,
		 protProtMap& stats )
{
  // the json is sparse, but here we have to keep track of the zeroes
  // so we start by filling a matrix with zeroes
  doubleVector blanks(fbgnMap.size(), 0);
  doubleVectorVector scores(fbgnMap.size(), blanks);

  rj::Value::ConstMemberIterator p1, p2;
  size_t i1, i2;
  for (p1 = begin; p1 != end; ++p1) {
    i1 = boost::lexical_cast<size_t>(p1->name.GetString());
    //std::cout << "\t\t\tseek fbgnMap.right.at(" << i1 << ")\n"; // DEBUG
    for (p2 = p1->value.MemberBegin(); p2 != p1->value.MemberEnd(); ++p2) {
      i2 = boost::lexical_cast<size_t>(p2->name.GetString());

      if (i1 > i2) {
	std::cerr << "i1 > i2!!!" << std::endl;
      }
      //std::cout << "\t\t\tseeking scores.at(" << i1 << ").at(" << i2 << ")\n";
      scores.at(i1).at(i2) = p2->value.GetDouble();
    }
  }
  
  std::string prot1, prot2;
  for (size_t i1 = 0; i1 < fbgnMap.size()-1; ++i1) {
    for (size_t i2 = i1+1; i2 < fbgnMap.size(); ++i2) {
      //std::cout << "seek fbgnMap.right.at(" << i2 << ")\n"; // DEBUG
      prot1 = fbgnMap.right.at(i1);
      prot2 = fbgnMap.right.at(i2);
      // disambiguate stats[prot1][prot2] vs stats[prot2][prot1]
      if (prot1.compare(prot2) > 0) {
	std::swap(prot1, prot2);
      }
      //std::cout << "\t\t\tstats[" << prot1 << "][" << prot2 << "]\n"; // DEBUG
      if (stats.count(prot1) == 0 || stats.at(prot1).count(prot2) == 0) {
	// brackets automatically instantiate 
	stats[prot1][prot2] = RunningStat(nRun, scores.at(i1).at(i2));
      } else {
	stats.at(prot1).at(prot2).Push(scores.at(i1).at(i2));
      }
    }
  }

  return;
}
