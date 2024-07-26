#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' Locate a middle sequence pattern and count its sequence context.
//'
//' This function searches for a specified middle pattern within a given sequence. 
//' It then counts the occurrences of specific context patterns within a defined window
//' size around the middle pattern. The function returns a map where keys are the 
//' counts of context patterns found and values are the frequencies of these counts.
//'
//' @param sequence A string representing the sequence to be analyzed.
//' @param mid_pattern A string representing the middle pattern to search for within the sequence.
//' @param window An integer specifying the size of the surrounding window around the middle pattern.
//' @param context_patterns A vector of strings representing the context patterns to search for within the window.
//' @return A std::unordered_map<int,int> where keys are the counts of context patterns found 
//'         and values are the frequencies of these counts.
//'
//'
//' @export
// [[Rcpp::export]]
std::unordered_map<int,int> countMidPatternContext2(std::string sequence,
    std::string mid_pattern, int window,
    std::vector<std::string> context_patterns) {

  std::unordered_map<int,int> counts;
  std::vector<int> cnt_pos(sequence.length());
  int context_len = context_patterns[0].length();

  // Error checking
  for (int i = 0; i < int(context_patterns.size()); ++i) {
    if (int(context_patterns[i].length()) != context_len) {
      stop("Patterns must be similiar in length.");
    }
  }

  int counter1 = 0;
  int counter2 = 0;
  int win1 = window / 2;
  int win2 = (window + mid_pattern.length()) - win1;

  // |--win1---|----win2---|
  // ---------XX------------

  // The first is count and the second is window end position
  std::deque<std::pair<int,int>> mid_counter;

  // Find mid_idx to jump
  int mid_idx = sequence.find(mid_pattern);

  for (int i = std::max(0, mid_idx - win1), j = i, k = i;
      i < int(sequence.length()); ++i) {

    // Slide and count every position
    if (std::find(context_patterns.begin(), context_patterns.end(),
          sequence.substr(i, context_len)) != context_patterns.end()) {
      counter1++;
      counter2++;
      cnt_pos[i]++;
    }

    // Add right window
    if (!mid_counter.empty() && i == mid_counter[0].second -
        context_len + 1) {
      counts[mid_counter[0].first + counter2]++;
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
        mid_counter.push_back(std::make_pair(counter1,
              std::min(int(sequence.length()) - 1, i + win2)));
      }
    }

    // Sliding count adjustment for win1 - allow overflow
    if (i - j + 1 == win1) {
      counter1 -= cnt_pos[j];
      j++;
    }
    // Sliding count adjustment for win2
    if (i - k + 1 == win2 - context_len + 1) {
      counter2 -= cnt_pos[k];
      k++;
    }
  }

  return counts;
}
