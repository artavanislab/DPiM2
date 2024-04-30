#include "Network.hpp"
#include <iostream>
#include <string>

#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector< std::string > stringVector;
#endif

readList(stringVector * ret, const std::string inFile);

int main  (int argc, char* argv[]) {
  std::string jsonList = "/home/glocke/DPiM/dpim4/consTest0_01-11-2016/json.list";
  stringVector  * jsonFiles = new stringVector;
  readList(jsonFiles, jsonList);
  Network n(*jsonFiles);
}

int readList(stringVector * ret, const std::string inFile) {

  std::string line;
  std::ifstream myFile;
  myFile.open(inFile.c_str());
  
  if (!myFile.is_open()) {
    std::cerr << "failed to open list file " << inFile << std::endl;
    return -1;
  }
  
  while (! myFile.eof() ) {
    getline (myFile,line);
    if (line[0]=='#') { continue; } // comments
    if (line.length()==0) { continue; } // eof

    //boost::algorithm::trim_right(line); chomp apparently unnecessary
    ret->push_back(line);
  }

  return 0;
}

