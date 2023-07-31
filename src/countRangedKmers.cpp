#include <Rcpp.h>
#include "countRangedKmers.h"

//' Count k-mers in given ranges of a sequence.
//'
//' Slide and update the cummulated table count.
//'
//' @param sequence A sequence to count.
//' @param starts Start positions.
//' @param ends End positions.
//' @param k K-mer size.
//' @return A k-mer-named vector of count.
//'
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string,int> countRangedKmers(std::string sequence,
    std::vector<int> starts, std::vector<int> ends, int k) {
  
  // The first is kmer and the second is the count.
  std::unordered_map<std::string,int> counts;
  
  // Create pair of start (adjust to 0-based index) and end to sort
  std::deque<std::pair<int,int>> coor(starts.size());
  for (int i = 0; i < int(starts.size()); ++i) {
    coor[i] = std::make_pair(starts[i] - 1, ends[i]);
  }
  std::sort(coor.begin(), coor.end());
  
  // Range to loop in sequence i.e. min start and max end
  int min_start = *std::min_element(starts.begin(), starts.end()) - 1;
  int max_end = *std::max_element(ends.begin(), ends.end()) - k + 1;
  
  int i = min_start;
  while (i < max_end) {
    
    // naive approach - loop every coor records
    // for (int j = 0; j < int(coor.size()); ++j) {
    //   if (i >= coor[j].first & i < coor[j].second) {
    //     std::string kmer = sequence.substr(i, k);
    //     counts[kmer]++;
    //   }
    // }
    
    // Loop every record
    // Break if sequence index is below coordinate start
    // Erase coordinate that has past sequence index
    auto it = coor.begin();
    while (it != coor.end()) {
      if (i < it->first) {
        break;
      } else if ((i >= it->first) & (i < it->second - k + 1)) {
        std::string kmer = sequence.substr(i, k);
        counts[kmer]++;
        ++it;
      } else if (it == coor.begin() && i >= it->second - k + 1) {
        // If std::vector cannot handle erase, change to list container.
        // 'it' will be set to the next element in coor
        //it = coor.erase(it);
        coor.pop_front();
        it = coor.begin();
      } else if (i >= it->second - k + 1) {
        ++it;
      }
    }
    
    // Jump if sequence index is lower than the next coordinate.
    if (!coor.empty() && i < coor[0].first) {
      i = coor[0].first;
    } else {
      i++;
    }
    
  }

  return counts;
}
