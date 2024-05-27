#include <Rcpp.h>

//' Count k-mers in given ranges of a sequence.
//'
//' Slide and update the cummulated table count.
//'
//' @param sequence A sequence to count.
//' @param starts Start positions.
//' @param ends End positions.
//' @param k K-mer size.
//' @param mid_pattern A middle pattern to search for.
//' @return A k-mer-named vector of count.
std::unordered_map<std::string,int> countMidPatternRangedKmers(
    std::string sequence, std::vector<int> starts, std::vector<int> ends, int k,
    std::string mid_pattern) {

  // The first is kmer and the second is the count.
  std::unordered_map<std::string,int> counts;

  int len_flank = (k - mid_pattern.size()) / 2;

  // Create pair of start (adjust to 0-based index) and end to sort
  // Additional adjustment based on middle pattern position
  std::deque<std::pair<int,int>> coor(starts.size());
  for (int i = 0; i < int(starts.size()); ++i) {
    coor[i] = std::make_pair(starts[i] - 1 + len_flank,
                             ends[i] - len_flank - mid_pattern.size() + 1);
  }
  std::sort(coor.begin(), coor.end());

  // Range to loop in sequence i.e. min start and max end
  int min_start = *std::min_element(starts.begin(), starts.end()) - 1 +
    len_flank - 1;
  int max_end = *std::max_element(ends.begin(), ends.end()) - len_flank -
    mid_pattern.size() + 1;

  // Find mid_idx to jump
  int mid_idx = sequence.find(mid_pattern, min_start);

  while ((mid_idx < max_end) & (mid_idx != int(std::string::npos))) {

    // Loop every record
    // Break if sequence index is below coordinate start
    // Erase coordinate that has past sequence index
    auto it = coor.begin();
    while (it != coor.end()) {
      if (mid_idx < it->first) {
        break;
      } else if ((mid_idx >= it->first) && (mid_idx < it->second)) {
        std::string kmer = sequence.substr(mid_idx - len_flank, k);
        counts[kmer]++;
        ++it;
      } else if (it == coor.begin() && mid_idx >= it->second) {
        // If std::vector cannot handle erase, change to list container.
        // 'it' will be set to the next element in coor
        //it = coor.erase(it);
        coor.pop_front();
        it = coor.begin();
      } else if (mid_idx >= it->second) {
        ++it;
      }
    }

    // Jump if sequence index is lower than the next coordinate.
    // string_view would be better to prevent overflow.
    if (!coor.empty() && mid_idx < coor[0].first) {
      mid_idx = sequence.find(mid_pattern, coor[0].first);
    } else {
      mid_idx = sequence.find(mid_pattern, mid_idx + 1);
    }

  }

  return counts;
}