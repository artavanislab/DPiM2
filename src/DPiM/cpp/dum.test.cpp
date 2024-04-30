  {   // DEBUG 2016-12-19
    cerr << "search_id " << endl;
    //for(auto const& it : searchId2Idx) {
    typedef std::map<std::string, int>::iterator sid_iter;

    for (sid_iter it = searchId2Idx.begin(); it!= searchId2Idx.end();
	 ++it) {
      string sid = it->first;
      int sidIdx = it->second;
      int tscSum = std::accumulate(searchByPrey.at(sidIdx).begin(),
				   searchByPrey.at(sidIdx).end(), 0);
      cerr << sid << "\t" << searchByPrey.at(sidIdx).size() << "\t"
	   << tscSum << endl;
    /*
    */  
    }
    exit(0);
  }
