#include <Rcpp.h>
using namespace Rcpp;

//' Count sequence content in a sliding window of a sequence.
//'
//' @param sequence A sequence to slide.
//' @param window Size of a window.
//' @param pattern A pattern to search for.
//' @return A numeric vector of count.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> countSlidingWindow(std::string sequence, int window,
    std::string pattern) {

  std::vector<int> counts(sequence.length() - window + 1);
  std::vector<int> cnt_pos(sequence.length());

  int counter = 0;
  int idx = 0;

  for (int i = 0; i < int(sequence.length() - pattern.length()) + 1; ++i) {

    if (sequence.substr(i, pattern.length()) == pattern) {
      counter++;
      cnt_pos[i]++;
    }

    // Start assigning
    if (i - idx + 1 == window - int(pattern.length()) + 1) {
      counts[idx] = counter;
      counter -= cnt_pos[idx];
      idx++;
    }
  }

  return counts;
}
