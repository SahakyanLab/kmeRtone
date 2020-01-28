#include <Rcpp.h>
using namespace Rcpp;

// This function calculate the "TRUE" percentage at a given coordinate start and end.
//    - The bool vector is result of single base pattern matching (e.g. G|C, G, etc.) 
//      of a DNA sequence.
//    - Boolean vector is smaller in size compare to DNA sequence. e.g. chr1 size is 2 GB but
//      it's G|C boolean vector is only 900 MB

// [[Rcpp::export]]
NumericVector countWindowBool(LogicalVector bool_vector, NumericVector start, NumericVector end) {
  
  // total DNA window
  int n = start.length() ;
  
  // initialise bool_percent vector
  NumericVector bool_percent( n ) ;
  
  // minus one to follow c++ indexing style
  start = start - 1 ;
  end = end - 1 ;
  
  // calculate the percentage
  for( int i=0; i<n; i++) {
    bool_percent[i] = sum(bool_vector[ Range(start[i], end[i]) ]) / (end[i]-start[i]+1) * 100 ;
  }
  
  return bool_percent ;
}