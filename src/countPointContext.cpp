#include <Rcpp.h>
using namespace Rcpp;

//' Ccount sequence context of given point positions.
//'
//' @param sequence A sequence to slide.
//' @param points Middle point positions.
//' @param len Length of the middle point. Default to 1.
//' @param window Size of a surrounding window.
//' @param context_pattern A context pattern to search for.
//' @return A numeric vector of count.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> countPointContext(std::string sequence,
    std::vector<int> points, int len, int window,
    std::string context_pattern) {

  // Sort and save original index
  // Initially sorted_counts contain window ends in zero-based
  int idx = 0;
  std::vector<std::pair<int,int>> sorted_counts;
  for (int i = 0; i < int(points.size()); ++i) {
    sorted_counts.push_back(std::make_pair(std::min(int(sequence.length()),
            points[i] + len + window/2 - 1 - 1), i));
  }
  std::sort(sorted_counts.begin(), sorted_counts.end());

  std::vector<int> cnt_pos(sequence.length());
  int counter = 0;

  int min_start = std::max(0, sorted_counts[0].first - window - len + 1);
  int max_end = sorted_counts[sorted_counts.size() - 1].first;

  for (int i = min_start, j = min_start; i <= max_end; ++i) {

    // Slide and count every position
    if (i <= int(sequence.length() - context_pattern.length()) + 1) {
      if (sequence.substr(i, context_pattern.length()) == context_pattern) {
        counter++;
        cnt_pos[i]++;
      }
    }

    // Add count when reaching window end
    if (i == sorted_counts[idx].first - int(context_pattern.length()) + 1) {
      sorted_counts[idx].first = counter;
      if (idx == int(sorted_counts.size()) - 1) {
        break;
      }

      idx++;
      // Jump - adjust i and j if the next idx is not overlap i.e. far away
      if (i < sorted_counts[idx].first - window - len + 1) {
        i = sorted_counts[idx].first - window - len; // main loop will add one.
        j = i + 1;
        counter = 0;
      }
    }

    // Sliding count adjustment
    if (i - j + 1 == window + len - int(context_pattern.length()) + 1) {
      counter -= cnt_pos[j];
      j++;
    }
  }

  std::vector<int> counts(points.size());
  for (int i = 0; i < int(points.size()); ++i) {
    counts[sorted_counts[i].second] = sorted_counts[i].first;
  }

  return counts;
}
