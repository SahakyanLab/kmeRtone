#include <Rcpp.h>
using namespace Rcpp;

#include <omp.h>
#include <regex>
#include "count_substring_regex.h"

// [[Rcpp::plugins(cpp11)]]
//' Count sequence content in a given sequence.
//'
//' stringi has no function that search within substring without memory creating
//' it. This function solve that. Unlike count_substring_fixed, this function
//' does not need to memory copy substring.
//'
//' @param sequence A sequence to map.
//' @param start Start positions.
//' @param end End positions.
//' @param pattern A regex pattern to search for within start and end positions.
//' @return A numeric vector of count.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> count_substring_regex(std::string sequence,
    std::vector<int> start, std::vector<int> end, std::string pattern) {

  int vec_sz = start.size();
  std::vector<int> cnt(vec_sz);
  std::regex rgx(pattern);

  if (pattern == "") {
    stop("Empty search patterns are not supported");
  }

  // Convert to zero-based position
  for (int i = 0; i < vec_sz; ++i) {
    start[i]--; end[i]--;
    if (end[i] >= int(sequence.length())) {
      end[i] = sequence.length() - 1;
    }
    if (start[i] < 0) {
      start[i] = 0;
    }
  }

  #pragma omp parallel for
  for (int i = 0; i < vec_sz; ++i) {

    auto idx_start = sequence.begin() + start[i];
    auto idx_end = sequence.begin() + end[i] + 1;
    std::ptrdiff_t match_num = std::distance(
        std::sregex_iterator(idx_start, idx_end, rgx),
        std::sregex_iterator());
    cnt[i] = match_num;

  }
  return cnt;
}
