#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
//' Count relavent k-mers with specified middle pattern from sequence string(s)
//' using a simple hash table.
//'
//' Locate a middle sequence pattern and register and count k-mer.
//'
//' @param sequence A sequence to slide.
//' @param k K-mer size.
//' @param mid_pattern A middle pattern to search for.
//' @return A k-mer-named vector of count.
//'
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string,int> countMidPatternKmers(
    std::vector<std::string> sequences, int k, std::string mid_pattern) {

  // The first is kmer and the second is the count.
  std::unordered_map<std::string,int> counts;

  int len_flank = (k - mid_pattern.size()) / 2;

  for (int i = 0; i < int(sequences.size()); ++i) {

    // Find mid_idx to jump
    int mid_idx = sequences[i].find(mid_pattern, len_flank);

    while((mid_idx <= int(sequences[i].length() - len_flank -
      mid_pattern.size())) & (mid_idx != int(std::string::npos))) {

      counts[sequences[i].substr(mid_idx - len_flank, k)]++;

      mid_idx = sequences[i].find(mid_pattern, mid_idx + 1);

    }

  }

  return counts;
}