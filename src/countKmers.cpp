#include <Rcpp.h>
using namespace Rcpp;

//' Count k-mers from string(s) using a simple hash table.
//' 
//' Count only observed k-mers. Biostrings::oligonucleotideFrequency reports all
//' possible k-mers. For k > 12, the memory for creating empty k-mer counts
//' spiked and crashed the R session.
//'
//' @param sequences Sequence strings.
//' @param k Size of k-mer.
//' @return A vector of k-mer counts. The counts of multiple sequences are
//'    combined, similar to Biostrings::oligonucleotideFrequency simplify.as
//'    "collapsed".
//'
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string,int> countKmers(
    std::vector<std::string> sequences, int k) {
  
  // The first is kmer and the second is the count.
  std::unordered_map<std::string,int> counts;
  
  for (int i = 0; i < int(sequences.size()); ++i) {
    for (int j = 0; j < int(sequences[i].length()) - k; ++j) {
      counts[sequences[i].substr(j, k)]++;
    }
  }
  
  return counts;
}