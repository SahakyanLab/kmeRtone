#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
//' Count Relevant K-mers with Specified Middle Pattern from Sequence String(s)
//'
//' This function scans through each sequence in the provided vector, locating a specified middle pattern.
//' For each occurrence of the middle pattern, the function extracts and counts the surrounding k-mers. 
//' The k-mers are identified based on the given k-mer size and centered around the middle pattern.
//'
//' @param sequences A vector of strings, each representing a sequence to be analyzed.
//' @param k An integer specifying the size of the k-mers to be extracted and counted.
//' @param mid_pattern A string representing the middle pattern to search for within each sequence.
//' @return A std::unordered_map with k-mers as keys and their counts as values.
//'
//' @examples
//' sequences <- c("ATCGATCGA", "GCGCATGCA")
//' k <- 5
//' mid_pattern <- "CG"
//' countMidPatternKmers(sequences, k, mid_pattern)
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