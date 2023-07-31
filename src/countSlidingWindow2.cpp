#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' Count sequence content in a sliding window of a sequence.
//'
//' @param sequence A sequence to slide.
//' @param window Size of a window.
//' @param patterns Patterns of the same size to search for.
//' @return Named vector of frequency of count.
//'
//' @export
// [[Rcpp::export]]
std::unordered_map<int,int> countSlidingWindow2(std::string sequence,
    int window, std::vector<std::string> patterns) {

  std::unordered_map<int,int> counts;
  std::vector<int> cnt_pos(sequence.length());
  int len = patterns[0].length();

  // Error checking
  for (int i = 0; i < int(patterns.size()); ++i) {
    if (int(patterns[i].length()) != len) {
      stop("Patterns must be similiar in length.");
    }
  }

  int counter = 0;

  for (int i = 0, j = 0; i < int(sequence.length()) - len + 1; ++i) {

    if (std::find(patterns.begin(), patterns.end(), sequence.substr(i, len))
          != patterns.end()) {
      counter++;
      cnt_pos[i]++;
    }

    // Start assigning
    if (i - j + 1 == window - len + 1) {
      counts[counter]++;
      counter -= cnt_pos[j];
      j++;
    }
  }

  return counts;
}
