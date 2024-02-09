#include <Rcpp.h>
using namespace Rcpp;

#include "count_substring_fixed.h"

// [[Rcpp::plugins(openmp)]]
//' Count sequence content in a given sequence.
//'
//' stringi has no function that search within substring without memory copy it.
//' This function has two versions. One without the need to memory copy denoted
//' as ***. The only downside to this is std::string::find cannot stop searching
//' past end of substring. I manage to at least stop it as soon as possible. If
//' the pattern is long and rare, it won't stop until it find post-substring
//' pattern. The other version is memory copy substring but as this operation is
//' in the loop, the memory is still within comfortable range. c++17 has
//' std::string_view that solve this but still new and not widely available. Use
//' count_substring_regex to avoid memory copy.
//'
//' @param sequence A sequence to map.
//' @param start Start positions.
//' @param end End positions.
//' @param pattern A pattern to search for.
//' @return A numeric vector of count.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> count_substring_fixed(std::string sequence,
    std::vector<int> start, std::vector<int> end, std::string pattern) {

  int vec_sz = start.size();
  std::vector<int> cnt(vec_sz);

  //                                    ONLY for ***
  // Convert to zero-based position and adjust end to not go over pattern size.
  for (int i = 0; i < vec_sz; ++i) {
    start[i]--;
    end[i]--;
    //***end[i] = end[i] - pattern.length() + 1 - 1;
    if (end[i] >= int(sequence.length())) {
      end[i] = sequence.length() - 1;
    }
    if (start[i] < 0) {
      start[i] = 0;
    }
  }

  for (int i = 0; i < vec_sz; ++i) {

    //***while ((start[i] = sequence.find(pattern, start[i])) !=
    int idx = 0;
    std::string substr = sequence.substr(start[i], end[i]-start[i]+1);
    while ((idx = substr.find(pattern, idx)) !=
        int(std::string::npos)) {

     //***if (start[i] > end[i]) {
     //***  break;
     //***}

      cnt[i]++;

     //***if (start[i] == end[i]) {
     //***  break;
     //***}

      //***start[i]++;
      idx++;
    }
  }
  return cnt;
}
