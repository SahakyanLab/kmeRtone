#include <Rcpp.h>
using namespace Rcpp;

//' Locate a middle sequence pattern and count its sequence context.
//'
//' @param sequence A sequence to slide.
//' @param mid_pattern A middle pattern to search for.
//' @param window Size of a surrounding window.
//' @param context_pattern A context pattern to search for.
//' @return A numeric vector of count.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> countMidPatternContext(std::string sequence,
    std::string mid_pattern, int window, std::string context_pattern) {

  std::vector<int> counts;
  std::vector<int> cnt_pos(sequence.length());

  int counter1 = 0;
  int counter2 = 0;
  int win1 = window / 2;
  int win2 = (window + mid_pattern.length()) - win1;

  // The first is count and the second is window end position
  std::deque<std::pair<int,int>> mid_counter;

  // Find mid_idx to jump
  int mid_idx = sequence.find(mid_pattern);

  for (int i = std::max(0, mid_idx - win1), j = i, k = i;
      i < int(sequence.length()); ++i) {

    // Slide and count every position
    if (i <= int(sequence.length() - context_pattern.length()) + 1) {
      if (sequence.substr(i, context_pattern.length()) == context_pattern) {
        counter1++;
        counter2++;
        cnt_pos[i]++;
      }
    }

    // Add right window
    if (!mid_counter.empty() && i == mid_counter[0].second -
        int(context_pattern.length()) + 1) {
      counts.push_back(mid_counter[0].first + counter2);
      mid_counter.pop_front();

      // If there is no other mid_counter, try to jump.
      if (mid_counter.empty()) {
        mid_idx = sequence.find(mid_pattern, i);
        if (mid_idx == int(std::string::npos)) {
          break;
        } else if (mid_idx - win1 > i) {
          i = mid_idx - win1 - 1; // -1 for main loop ++i
          j = k = i + 1;
          counter1 = counter2 = 0;
        }
      }
    }

    // Check for mid_pattern
    if (i + 1 <= int(sequence.length() - mid_pattern.length())) {
      if (sequence.substr(i + 1, mid_pattern.length()) == mid_pattern) {
        mid_counter.push_back(std::make_pair(counter1, i + win2));
      }
    }

    // Sliding count adjustment for win1 - allow overflow
    if (i - j + 1 == win1) {
      counter1 -= cnt_pos[j];
      j++;
    }
    // Sliding count adjustment for win2
    if (i - k + 1 == win2 - int(context_pattern.length()) + 1) {
      counter2 -= cnt_pos[k];
      k++;
    }
  }

  return counts;
}
