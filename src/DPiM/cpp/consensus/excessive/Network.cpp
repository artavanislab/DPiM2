#include "Network.hpp"
#include "rapidjson/document.h"
#include <fstream>
#include <iostream>
#include <boost/lexical_cast.hpp>
//#include <boost/algorithm/string.hpp>

#ifndef INT_STRING_MAP
#define INT_STRING_MAP 1
typedef std::map< int, std::string > intStringMap;
#endif

namespace rj = rapidjson;

Network::Network() { };

Network::Network(const stringVector &jsonFiles) {
  fbgnMap = new stringIBimap;
  scoreHist = new doubleVectorPVectorPVector;
  
  stringVector::const_iterator f = jsonFiles.begin();
  firstJson(*f); // initilize various things

  ++f;
  for (; f != jsonFiles.end(); ++f) {
    addJson(*f);
  }
}
  

// return all scores for this pair of proteins
// returning a pointer is potentially unsafe, but probably that doesn't matter
doubleVector*
Network::pairScores(const std::string prot1, const std::string prot2) {
  size_t i1 = fbgnMap->left.at(prot1);
  size_t i2 = fbgnMap->left.at(prot2);
  return scoreHist->at(i1)->at(i2);
  // scoreHist->at(prot1).at(prot2)->at(run) = score of this interaction
}

void Network::firstJson(const std::string jsonFile) {
  char* json;
  {
    std::ifstream t;
    int length;
    t.open(jsonFile.c_str());        
    t.seekg(0, std::ios::end); // go to the end
    length = t.tellg();        // report location (this is the length)
    t.seekg(0, std::ios::beg); // go back to the beginning
    json = new char[length];   // 
    t.read(json, length);      // read the whole file into the buffer
    t.close();                 
  }
  rj::Document jsonObj;
  jsonObj.Parse(json);
  delete json;
  
  assert(jsonObj.IsObject());
  assert(jsonObj.HasMember("score"));
  assert(jsonObj["score"].IsObject());
  assert(jsonObj.HasMember("ids"));
  assert(jsonObj["ids"].IsObject());

  fillFbgnMap(fbgnMap, jsonObj["ids"].MemberBegin(),
	      jsonObj["ids"].MemberEnd());

  nProt = fbgnMap->size();

  /*
   * initialize the score matrix
   */
  doubleVectorPVector * scoreMat = new doubleVectorPVector;
  iIMap identity;
  for (size_t i = 0; i < nProt; ++i) {
    scoreMat->push_back( new doubleVector(nProt, 0) );
    identity.at(i) = i;
  }

  // fill it
  fillScoreMat(scoreMat, jsonObj["score"].MemberBegin(),
	       jsonObj["score"].MemberEnd(), identity);

  // add to history
  scoreHist->push_back(scoreMat);
  
  std::cout << "I love peanuts!\n";
}
void Network::addJson(const std::string jsonFile) {
}

void Network::fillFbgnMap(stringIBimap *idMap,
		 rj::Value::ConstMemberIterator begin,
		 rj::Value::ConstMemberIterator end)
{
  typedef stringIBimap::value_type keyValue;
  std::string fbgn;
  size_t id;
  rj::Value::ConstMemberIterator p1;
  for (p1 = begin; p1 != end; ++p1) {
    id = boost::lexical_cast<size_t>(p1->name.GetString());
    fbgn = p1->value.GetString();
    idMap->insert(keyValue(fbgn, id));
  }
  /*
    std::cout << "idMap[FBgn0000032] = " << idMap->left.at("FBgn0000032")
    << std::endl;
    std::cout << "idMap[FBgn0000039] = " << idMap->left.at("FBgn0000039")
    << std::endl;
  */
  return;
}

void Network::fillScoreMat(doubleVectorPVector *scoreMat,
		  rj::Value::ConstMemberIterator begin,
		  rj::Value::ConstMemberIterator end,
		  const iIMap &indexMap )
{
  rj::Value::ConstMemberIterator p1, p2;
  size_t i1, i2, tmp1, tmp2;
  double s;
  for (p1 = begin; p1 != end; ++p1) {
    tmp1 = boost::lexical_cast<size_t>(p1->name.GetString());
    i1 = indexMap.at(tmp1);
    for (p2 = p1->value.MemberBegin(); p2 != p1->value.MemberEnd(); ++p2) {
      tmp2 = boost::lexical_cast<size_t>(p2->name.GetString());
      i2 = indexMap.at(tmp2);
      s = p2->value.GetDouble();
      scoreMat->at(i1)->at(i2) = s;
    }
  }
}
