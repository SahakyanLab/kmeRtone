#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' Ccount sequence context of given point positions.
//'
//' @param sequence A sequence to slide.
//' @param points Middle point positions.
//' @param len Length of the middle point.
//' @param window Size of a surrounding window.
//' @param context_patterns Context patterns to search for.
//' @return A named vector of frequency of counts.
//'
//' @export
// [[Rcpp::export]]
std::unordered_map<int,int> countPointContext2(std::string sequence,
    std::vector<int> points, int len, int window,
    std::vector<std::string> context_patterns) {

  // Error checking
  for (int i = 0; i < int(points.size()); ++i) {
    if ((points[i] < 1) | (points[i] > int(sequence.length()))) {
      stop("Input points are out of range.");
    }
  }

  std::unordered_map<int,int> counts;
  std::sort(points.begin(), points.end());

  int context_len = context_patterns[0].length();

  // Impose restriction to patterns' length to be the same'
  for (int i = 0; i < int(context_patterns.size()); ++i) {
    if (int(context_patterns[i].length()) != context_len) {
      stop("Patterns must be similiar in length.");
    }
  }

  std::vector<int> cnt_pos(sequence.length());
  int counter = 0;

  int idx = 0;
  int min_start = std::max(0, points[0] - window/2 + 1 - 1);
  int max_end = std::min(int(sequence.length()) - 1,
      points[points.size() - 1] + len + window/2 - 1 - 1);

  for (int i = min_start, j = min_start; i <= max_end - context_len + 1; ++i) {

    // Slide and count every position
    if (std::find(context_patterns.begin(), context_patterns.end(),
          sequence.substr(i, context_len)) != context_patterns.end()) {
      counter++;
      cnt_pos[i]++;
    }


    // Add count when reaching window end
    while (i == std::min(int(sequence.length()) - 1, points[idx] + len +
            window/2 - 1 - 1) - context_len + 1) {

      counts[counter]++;

      if (idx == int(points.size()) - 1) {
        break;
      }

      idx++;
      // Jump - adjust i and j if the next idx is not overlap i.e. far away
      if (i < points[idx] - window/2 - 1) {
        i = points[idx] - window/2 - 2; // main loop will add one.
        j = i + 1;
        counter = 0;
      }
    }

    // Sliding count adjustment
    if (i - j + 1 == window + len - context_len + 1) {
      counter -= cnt_pos[j];
      j++;
    }
  }

  return counts;
}
