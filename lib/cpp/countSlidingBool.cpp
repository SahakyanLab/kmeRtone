#include <Rcpp.h>
using namespace Rcpp;

// This function calculate the "TRUE" percentage of sliding windows of a boolean vector.
//    - The boolean vector is result of single base pattern matching (e.g. G|C, G, etc.) 
//      of a DNA sequence.
//    - width is the size of the sliding windows
//    - Boolean vector is smaller in size compare to DNA sequence. e.g. chr1 size is 2 GB but
//      it's G|C boolean vector is only 900 MB

// [[Rcpp::export]]
NumericVector countSlidingBool(LogicalVector bool_vector, int width) {
  
  // total DNA sliding windows
  int n = bool_vector.length() - width + 1 ;
  
  // initialise bool_percent vector
  NumericVector bool_count( n ) ;
  
  // calculate the percentage
  for( int i=0; i<n; i++) {
    bool_count[i] = sum(bool_vector[ Range(i, i+width-1) ]) ;
  }
  
  return bool_count ;
}