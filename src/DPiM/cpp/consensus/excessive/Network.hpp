/* 
 * File:   Network.hpp
 * Author: glocke
 * Description: read in the json
 * 
 * Created on July 20, 2010, 5:50 PM01-12-2016
 * Last Modified July 20, 2010, 5:50 PM
 */

#ifndef _NETWORK_HPP
#define	_NETWORK_HPP 1

#include "rapidjson/document.h"
#include <vector>
#include <string>
#include <map>
#include <boost/bimap.hpp>

#ifndef STRING_I_BIMAP
#define STRING_I_BIMAP 1
typedef boost::bimap< std::string, size_t > stringIBimap;
#endif
#ifndef I_I_MAP
#define I_I_MAP 1
typedef std::map< size_t, size_t > iIMap;
#endif
#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector< double > doubleVector;
#endif
#ifndef DOUBLE_VECTORP_VECTOR
#define DOUBLE_VECTORP_VECTOR 1
typedef std::vector< doubleVector* > doubleVectorPVector; 
#endif
#ifndef DOUBLE_VECTORP_VECTORP_VECTOR
#define DOUBLE_VECTORP_VECTORP_VECTOR 1
typedef std::vector< doubleVectorPVector* > doubleVectorPVectorPVector;
#endif
#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector< std::string > stringVector;
#endif

class Network {
public:
  Network();
  Network(const stringVector &jsonFiles);
  

  // return all scores for this pair of proteins
  doubleVector* pairScores(const std::string prot1, const std::string prot2);
  
private:
  stringIBimap* fbgnMap;
  size_t nProt;
  
  doubleVectorPVectorPVector * scoreHist; // Hist as in History
  // scoreHist->at(prot1).at(prot2)->at(run) = score of this interaction
  // put pointers at the end so that scoreHist[i][j] equals [j][i] without
  // storing the same data twice.

  //doubleVector blankVector; // for constructing scoreHist
  
  void firstJson(const std::string jsonFile);
  void addJson(const std::string jsonFile);

  void fillFbgnMap(stringIBimap* idMap,
		   rapidjson::Value::ConstMemberIterator begin,
		   rapidjson::Value::ConstMemberIterator end);
  void fillScoreMat(doubleVectorPVector* scoreMat,
		    rapidjson::Value::ConstMemberIterator begin,
		    rapidjson::Value::ConstMemberIterator end,
		    const iIMap &indexMap );

  void resizeScoreHist(doubleVectorPVector* scoreMat);
  
  int readList(stringVector *ret, const std::string file);
  
};

#endif	/* _JSONNET_HPP */

/*
Network::
fillFbgnMap(boost::bimaps::bimap<std::string, unsigned long, mpl_::na, mpl_::na, mpl_::na>&,
	    rapidjson::GenericMemberIterator<true, rapidjson::UTF8<char>, rapidjson::MemoryPoolAllocator<rapidjson::CrtAllocator> >,
	    rapidjson::GenericMemberIterator<true, rapidjson::UTF8<char>, rapidjson::MemoryPoolAllocator<rapidjson::CrtAllocator> >*/
