/* 
 * File:   BaitByPrey.cpp
 * Author: jmintser
 * 
 * Created on July 20, 2010, 5:50 PM
 */

#include "BaitByPrey.hpp"

BaitByPrey::BaitByPrey(string filename) {
  loadData(filename);
  verbose=false;
}


BaitByPrey::BaitByPrey(string filename, string len_corr_file_name) {
  loadData(filename);
  loadPreyLength(len_corr_file_name);
  verbose=false;
}

void BaitByPrey::print () {
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      cout << (int) searchByPrey[i][j] << " ";
    }
    cout << endl;
  }
}

/**
 * Load length correction numbers from a tab-delimited file
 * with the format "#fbgn\tavg_exp_tryptic\tavg_length"
 * If the value happens to be 0 (as when using predicted # of tryptic peptides)
 * record a value of 1 instead
 */
void BaitByPrey::loadPreyLength(string input_file) {
  string prot_ref;
  double avg_exp_tryp_peps;
  double length;

  ifstream lenFile(input_file.c_str());
  if (!lenFile.is_open()) {
    cout << "Error opening length correction file " << input_file << endl;
    exit(0);
    return;
  }
  string line;
  while (getline(lenFile, line, '\n')) {
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    vector<stringstream*> field;
    string token;
    while (getline(ss, token, '\t')) {
      field.push_back(new stringstream(token));
    }

    *field[0] >> prot_ref;
    *field[1] >> avg_exp_tryp_peps;
    *field[2] >> length;
    if ((*field[1]).fail() || length == 0) {
      // glocke replaced replaced | with || in the above code
      length = 1;
      avg_exp_tryp_peps = 1;
    }

    if (preyLengthMap.count(prot_ref) == 0) {
      preyLengthMap[prot_ref] = length;
    } else {
      cout << "ERROR: Found multiple length correction values for " << prot_ref << endl << endl;
      exit(1);
    }
        
    // Clean-up memory to avoid nasty memory lea
    for(vector<stringstream *>::iterator inx = field.begin(); inx != field.end(); ++inx)
      delete *inx;
  }
  lenFile.close();
  cout << preyLengthMap.size() << " length correction values recorded" << endl; ////
  return;
}

/*
 * read the input file
 * 
 */
void BaitByPrey::loadData(string input_file, bool countGhostBait) {
  vector<string> searchIdVec;
  vector<int> totalPeptideVec;
  vector<string> preyReferenceVec;
  set<string> preyReferenceSet;
  set<string>::iterator setIt;

  int total_peptides;
  string tap_id, search_id, sample_date, bait_ref, prey_ref;

  
  ifstream ppiFile(input_file.c_str());
  if (!ppiFile.is_open()) {
    cout << "Error opening PPI file " << input_file << endl;
    exit(0);
    return;
  }
  string line;
  while (getline(ppiFile, line, '\n')) {
    if (line[0] == '#')
      continue;
        
    stringstream ss(line);
    vector<stringstream*> field;
    string token;
    while (getline(ss, token, '\t')) {
      field.push_back(new stringstream(token));
    }

    *field[1] >> search_id;
    *field[3] >> total_peptides;
    if ((*field[3]).fail()) {
      total_peptides = 1;
    }
    *field[4] >> bait_ref;
    *field[5] >> prey_ref;

    //cout << search_id << "\t" << total_peptides << "\t" << bait_ref << "\t" << prey_ref << endl; ////
    searchIdVec.push_back(search_id);
    totalPeptideVec.push_back(total_peptides);
    preyReferenceVec.push_back(prey_ref);

    preyReferenceSet.insert(prey_ref);
    if (countGhostBait) {
      preyReferenceSet.insert(bait_ref);
    }

    if (searchId2Idx.count(search_id) == 0) {
      //cout << search_id << "\t" << searchId2Idx.size() << endl;
      searchId2Idx.insert(idMapPair(search_id, searchId2Idx.size()));
      // \\=> glocke changed 07-24-2015, potential off-by-one error  
      searchId2BaitReference[search_id] = bait_ref;
      searchIdReferenceVec.push_back(search_id);
    } 
    for(vector<stringstream *>::iterator inx = field.begin(); inx != field.end(); ++inx)
      delete *inx;
  }
  ppiFile.close();
  
  // create the index map
  cout << preyReferenceSet.size() << " prey refs" << endl; ////
  cout << searchId2BaitReference.size() << " BaitRef search_ids" << endl; ////
  cout << searchId2Idx.size() << " search_ids" << endl; ////
  int cnt = 0;
  //uniqPreyReferenceVec.resize(preyReferenceSet.size(), 0);

  /**
   * This implementation of preyRef2Idx from preyReferenceSet ensures that
   * the index corresponds to essentially a numerically sorted FBgn ids.  
   * Combined with the fact that I am only using the upper triangle of the
   * preyByPrey matrix this then futher results in FBgn id pairs in the final 
   * output being sorted such that the FBgn_id1 is always smaller than FBgn_id2 
   * in every pair
   */
  for (setIt = preyReferenceSet.begin(); setIt != preyReferenceSet.end(); ++setIt) {
    //cout << *setIt << endl;
    preyRef2Idx[*setIt] = cnt;
    uniqPreyReferenceVec.push_back(*setIt);
    cnt++;
  }

  // fill in the search(bait) X prey matrix
  cout << "... creating a " << searchId2Idx.size() << " X " << preyReferenceSet.size() << " search(bait)-by-prey data matrix" << endl;
  searchByPrey.resize(searchId2BaitReference.size(), vector<double>( preyReferenceSet.size(), 0 ));
  cout << "searchByPrey complete\n";
  //bool first = true;
  for (unsigned int i = 0; i < searchIdVec.size(); i++) {
    
    // Debug start
    //cout << searchId2Idx[searchIdVec[i]] << "\t" << preyRef2Idx[preyReferenceVec[i]] << endl; ////
    //cout << searchIdVec[i] << endl;
    //if(searchId2Idx.find(searchIdVec[i]) == searchId2Idx.end()){
    //    cout << "Its not there" << endl;
    //}
    //else{
    //    cout << "Here it is " << searchId2Idx.find(searchIdVec[i])->second << endl;;
    //}
    // Debug end
        
    int sbp_i = searchId2Idx.at(searchIdVec.at(i));
    // sbp_i is the sbp row index for this search_id
    int sbp_j = preyRef2Idx.at(preyReferenceVec.at(i));
    // sbp_j is the sbp col index for this prey
    
    searchByPrey.at(sbp_i).at(sbp_j) += totalPeptideVec.at(i);
    
    /*
    // RECORD DIRECTLY OBSERVED INTERACTIONS
    string baitRef = searchId2BaitReference.at(searchIdVec.at(i));
    string preyRef = preyReferenceVec.at(i);
    if (searchByPrey.at(sbp_i).at(sbp_j) > 0) {
      if (preyRef > baitRef) {
	directBait_PreyPair.insert(baitRef + "\t" + preyRef);
      } else if (baitRef > preyRef) {
	directBait_PreyPair.insert(preyRef + "\t" + baitRef);
      } // ELSE BAIT == PREY
    }
    */
  }
  
  // add a fake prey count for each ghost bait
  if (countGhostBait) {
    for (map<string,int>::iterator mapIt = searchId2Idx.begin(); mapIt != searchId2Idx.end(); mapIt++) {
      if (searchByPrey[ (*mapIt).second ][ preyRef2Idx[ searchId2BaitReference[(*mapIt).first] ] ] == 0) {
	searchByPrey[ (*mapIt).second ][ preyRef2Idx[ searchId2BaitReference[(*mapIt).first] ] ] = 1;
	////cout << "...adding fake prey count for ghost bait " << searchId2BaitReference[(*mapIt).first] << endl; ////
      }
    }
  }

  shuffled = 0;
  return;
}

/**
 * Create a simulated randomized dataset that is a better null model than the output of shuffle()
 * Draws from a set of preys according the the real spectral count distribution.  Therefore, preys
 * that are abundant in the real dataset will also be abundant in the simulated one.
 * The number of uniq preys and total spectral counts for each run should equal the real dataset
 */
void BaitByPrey::simulateDataset() {
  cout << "... simulating " << shuffled << endl;
  // fill in prey spectral count vector
  if (totalSpecCountPreyIdxVec.size() == 0) {
    totalSpecCountPreyIdxVec.clear();
    totalSpecCountPreyIdxVec.reserve(searchByPrey.size()); // try to save some push_back allocation time
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      int col_sum = colSum(j);
      totalSpecCountPreyIdxVec.push_back(col_sum);
    }
  }
  Histogram specCountPreyHist(totalSpecCountPreyIdxVec);

  // fill in total search spectral count vector
  if (totalSpecCountSearchIdxVec.size() == 0) {
    for (unsigned int i = 0; i < searchByPrey.size(); i++) {
      totalSpecCountSearchIdxVec.push_back( rowSum(i) );
    }
  }

  // fill in uniq search spectral count vector
  if (uniqSpecCountSearchIdxVec.size() == 0) {
    for (unsigned int i = 0; i < searchByPrey.size(); i++) {
      int uniqPreys = 0;
      for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
	if (searchByPrey[i][j] > 0)
	  uniqPreys++;
      }
      uniqSpecCountSearchIdxVec.push_back(uniqPreys);
    }
  }

  // clear and re-fill with data drawn from specCountPreyIdxVec
  searchByPrey.clear();
  searchByPrey.resize(searchId2Idx.size(), vector<double>(uniqPreyReferenceVec.size(), 0));
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    int total_spec_count = totalSpecCountSearchIdxVec[i];
    int uniq_spec_count = uniqSpecCountSearchIdxVec[i];

    // draw random prey until have target number of uniques
    int cur_total_spec_counts = 0;
    map<int,int> cur_uniq_preyMap; // preyIdx => numTotSimulatedSpecCounts
    while ((int)cur_uniq_preyMap.size() < uniq_spec_count) {
      int cur_prey_idx = specCountPreyHist.sample();
      ////cout << cur_prey_idx << endl; ////
      cur_uniq_preyMap[cur_prey_idx] ++;
      cur_total_spec_counts ++;
    }
        
    // fill in a local specCountPreyIdxVec containing only relevant preys
    Histogram cur_spec_count_prey_hist(cur_uniq_preyMap); // will draw from this to finish filling in cur_uniq_preyMap

    // draw from current prey until have target number of total spec counts
    while (cur_total_spec_counts < total_spec_count) {
      int cur_prey_idx = cur_spec_count_prey_hist.sample();
      ////cout << cur_prey_idx << endl; ////
      cur_uniq_preyMap[cur_prey_idx] ++;
      cur_total_spec_counts ++;
    }

    // write back to searchByPrey[i]
    for (map<int,int>::const_iterator mapIt = cur_uniq_preyMap.begin(); mapIt != cur_uniq_preyMap.end(); mapIt++) {
      searchByPrey[i][(*mapIt).first] = (*mapIt).second;
      ////cout << i << "\t" << (*mapIt).first << "\t" << (*mapIt).second << endl; ////
    }
    ////exit(0);
  }
  shuffled++;
  return;
}


// Useful for outputting shuffled datasets
void BaitByPrey::outputRawScores() {
  cout << "#tap_id\tsearch_id\tsample_date\ttotal_peptides\tbait_ref\tprey_ref" << endl;
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      if (searchByPrey[i][j] == 0)
	continue;
      cout << "1\t"
	   << searchIdReferenceVec[i]
	   << "\t0000-00-00\t"
	   << (int) searchByPrey[i][j] << "\t"
	   << searchId2BaitReference[ searchIdReferenceVec[i] ] << "\t"
	   << uniqPreyReferenceVec[j] << endl;
    }
  }
  return;
}

void BaitByPrey::binary() {
  cout << "... re-scaling: Binary" << endl;
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {

      searchByPrey[i][j] = searchByPrey[i][j] > 0 ? 1 : 0;
    }
  }
  return;
}

/**
 * Convert the matrix spectral counts into Normalized Spectral Abundance Factors (NSAF)
 * Paoletti et al., PNAS (2006)
 * Essentially, this accounts for protein length
 */
void BaitByPrey::nsaf() {
  cout << "... re-scaling: NSAF" << endl;
  if (searchByPrey[0].size() > preyLengthMap.size()) {
    cout << "ERROR: Some or all length correction values appear to be missing!" << endl << endl;
    exit(1);
  }
  double row_sum;
  for (unsigned int i = 0; i < searchByPrey.size(); i++) { // loop over search_id
    // Normalize by length factor
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      if (preyLengthMap.count( uniqPreyReferenceVec[j] ) == 0) {
	cout << "ERROR: Missing length correction value for " << uniqPreyReferenceVec[j] << endl << endl;
	exit(1);
      }
      searchByPrey[i][j] /= preyLengthMap[ uniqPreyReferenceVec[j] ];
    }

    // Normalize by total run counts
    row_sum = rowSum(i);
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      searchByPrey[i][j] /= row_sum;
    }
  }

  // re-scale such that the smallest NSAF value is 1
  double mat_min = minMin();
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      //searchByPrey[i][j] = ceil(searchByPrey[i][j]/mat_min - 0.5);
      // GL removed rounding 1-7-2016
      searchByPrey[i][j] = searchByPrey[i][j]/mat_min;
    }
  }

  return;
}

void BaitByPrey::ln() {
  cout << "... re-scaling: Log" << endl;
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      searchByPrey[i][j] = (searchByPrey[i][j] <= 1) ? searchByPrey[i][j] : ceil(log(searchByPrey[i][j])); // ignore 0's and 1's
    }
  }
  return;
}

void BaitByPrey::logBase10() {
  cout << "... re-scaling: Log10" << endl;
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      searchByPrey[i][j] = (searchByPrey[i][j] <= 1) ? searchByPrey[i][j] : ceil(log10(searchByPrey[i][j])); // ignore 0's and 1's
    }
  }
  return;
}

void BaitByPrey::sqrt() {
  cout << "... re-scaling: Sqrt" << endl;
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      //searchByPrey[i][j] = ceil(pow(searchByPrey[i][j], 0.5) - 0.5);
      // GL removed rounding 1-7-2016
      searchByPrey[i][j] = pow(searchByPrey[i][j], 0.5);
    }
  }
  return;
}

// round results to the nearest integer
// this code relies on the fact that searchByPrey is non-negative
void BaitByPrey::round() {
  if (verbose) {
    cout << "... re-scaling: round to nearest integer" << endl;
  }
  for (unsigned int i = 0; i < searchByPrey.size(); i++) {
    for (unsigned int j = 0; j < searchByPrey[0].size(); j++) {
      searchByPrey[i][j] = ceil(searchByPrey[i][j] - 0.5);
    }
  }
  return;
}

void  BaitByPrey::reScale(RescaleType sc) {
  switch (sc) {
  case Binary:
    binary();
    break;
  case NSAF:
    nsaf();
    break;
  case Log:
    ln();
    break;
  case Log10:
    logBase10();
    break;
  case Sqrt:
    sqrt();
    break;
  case Round:
    round();
    break;
  default:
    cerr << "BaitByPrey::reScale doesn't recognize RescaleType\n";
  }
  return;
}

double BaitByPrey::rowSum(int sid) {
  return accumulate(searchByPrey[sid].begin(), searchByPrey[sid].end(), 0.0);
}

double BaitByPrey::sumSum() {
  double sum_sum = 0.0;
  for (vector< vector<double> >::const_iterator vecIt = searchByPrey.begin(); vecIt != searchByPrey.end(); vecIt++) {
    sum_sum += accumulate((*vecIt).begin(), (*vecIt).end(), 0.0);
  }
  return sum_sum;
}

double BaitByPrey::colSum(int prey_idx) {
  double col_sum = 0.0;
  for (vector< vector<double> >::const_iterator vecIt = searchByPrey.begin();
       vecIt != searchByPrey.end(); vecIt++)
  {
    col_sum += (*vecIt)[prey_idx];
  }
  return col_sum;
}

// Find smallest non-zero count in the matrix
double BaitByPrey::minMin(bool excludeZero) {
  double min_min = numeric_limits<double>::max();
  double row_min;
  for (vector< vector<double> >::const_iterator vecIt = searchByPrey.begin(); vecIt != searchByPrey.end(); vecIt++) {
    for (vector<double>::const_iterator vecIt2 = (*vecIt).begin(); vecIt2 != (*vecIt).end(); vecIt2++) {
      row_min = *vecIt2;
      if (row_min < min_min && (excludeZero && row_min > 0)) {
	min_min = row_min;
      }
    }
  }
  return min_min;
}

double BaitByPrey::maxMax() {
  double max_max = 0;
  // formerly numeric_limits<double>::min(), which is equivalent to zero for
  // these purposes
  double row_max;
  for (vector< vector<double> >::const_iterator vecIt = searchByPrey.begin(); vecIt != searchByPrey.end(); vecIt++) {
    row_max = *max_element((*vecIt).begin(), (*vecIt).end());
    if (row_max > max_max) {
      max_max = row_max;
    }
  }
  return max_max;
}


double BaitByPrey::vectorMean(vector<double> &vec) {
  double sum = accumulate(vec.begin(), vec.end(), 0.0);
  return sum/vec.size();
}


double BaitByPrey::vectorStd(vector<double> &vec) {
  double mean = vectorMean(vec);
  vector<double> mean_diff_squared = vec;
  for (unsigned int i = 0; i < mean_diff_squared.size(); i++) {
    mean_diff_squared[i] = pow((mean_diff_squared[i] - mean), 2);
  }

  double SSE = accumulate(mean_diff_squared.begin(), mean_diff_squared.end(), 0.0);
  return pow(SSE/(mean_diff_squared.size() - 1), 0.5);
}


BaitByPrey::~BaitByPrey() {
}
