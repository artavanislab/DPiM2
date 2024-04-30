#include "rapidjson/prettywriter.h" 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <map>
#include <limits>
#include <math.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>

using namespace std;
namespace po = boost::program_options;
namespace rj = rapidjson;

//using namespace rapidjson;


/*
 * notes on compiling
 * g++ dum.cxx -lboost_program_options  -I/home/glocke/DPiM/cpp/rapidjson/rapidjson-master/include
 * don't forget to add boost module
 */
int main(int argc, char** argv) {

  string inFile, outFile, rownamesFile;
  
  string rownamesDefault = "out.rownames.json";
  po::options_description desc("Usage");
  desc.add_options()
    ("in", po::value<string>()->required(),
     "REQUIRED: APMS input data")
    ("out", po::value<string>()->required(),
     "REQUIRED: write HGScore here (JSON)")
    ("rownames", po::value<string>()->default_value(rownamesDefault),
     "records the FBgn's corresponding to rows/cols of output")
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

  inFile = opts["in"].as<string>();
  outFile = opts["out"].as<string>();
  rownamesFile = opts["rownames"].as<string>();

  if (rownamesFile.compare(rownamesDefault)==0) {
    rownamesFile = outFile;
    boost::replace_last(rownamesFile, "json", "rownames.json");
    if (rownamesFile.compare(outFile)==0) {
      // if the above had no effect
      rownamesFile += ".rownames.json";
    }
  }

  vector< string > filenames;
  filenames.push_back(inFile);
  filenames.push_back(outFile);
  filenames.push_back(rownamesFile);
  
  //StringBuffer sb;
  //PrettyWriter<StringBuffer> writer(sb);
  rj::StringBuffer sb;
  rj::PrettyWriter<rj::StringBuffer> writer(sb);
  writer.StartObject();
  writer.String("a");
  writer.StartArray();
  for (vector<string>::iterator f=filenames.begin(); f!= filenames.end(); ++f) {
#if RAPIDJSON_HAS_STDSTRING
    writer.String(*f);
#else
    writer.String(f->c_str(), static_cast<rj::SizeType>(f->length())); // Supplying length of string is faster.
#endif
  }
  writer.EndArray();
  writer.String("b");
  writer.Uint(-1);
  writer.EndObject();

  ofstream out;
  out.open("test.out");
  out << sb.GetString() << std::endl;
  out.close();
  return 0;
  

  if (0) {
    double x = -100.5;
    int y = ceil(x - 0.5);
    cout << x << "\t" << y << "\n";
    return 0;
  }  

  {
    double m =std::numeric_limits<double>::min();
    double mLog = -log(m);
    double dm =std::numeric_limits<double>::denorm_min();
    double dmLog = -log(dm);
  
    std::cout << std::boolalpha;
    std::cout << "Min double: " << m << '\n';
    std::cout << "-log(Min double): " << mLog << '\n';

    std::cout << "denorm Min double: " << dm << '\n';
    std::cout << "-log(denorm Min double): " << dmLog << '\n';

    return 0;
  }
  
  map<string, int> myMap;
  cout << "size of myMap = " << myMap.size() << endl;

  for (int i=0; i < 10; ++i) {
    string s = "axk"; // not working? to_string(i);
    //string s = to_string(i);
    cout << "myMap.count(" << s << ") = " << myMap.count(s) << endl;
  }

  cout << "size of myMap = " << myMap.size() << endl;

  cout << "myMap.insert(\"x\", myMap.size())" << endl;

  typedef pair<string, int> idMapPair;
  myMap.insert(idMapPair("x", myMap.size()));
  cout << "myMap[\"x\"] = " << myMap["x"]  << endl;

  cout << "setting myMap[\"y\"] = myMap.size()" << endl;
  myMap["y"]=myMap.size();
  cout << "myMap[\"y\"] = " << myMap["y"]  << endl;
  cout << "size of myMap = " << myMap.size() << endl;
}
