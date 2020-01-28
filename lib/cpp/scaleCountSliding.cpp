#include <Rcpp.h>
using namespace Rcpp;

// This function calculate pattern matching count from scaled-up
// sliding window. It save times by reducing the length subsetting
// if we were to calculate from the beginning. Even c++ is slow
// for subsetting large vector e.g. DNA.sequence[ 1:100000 ]
// Parallelisation is possible but not memory efficient because a
// large count vector 1 need to be copy to each cluster and a large
// count vector 2 need to be created in each of them.

// [[Rcpp::export]]
NumericVector scaleCountSliding(NumericVector count_vector_1, int width_1, int width_2) {
  
  int scale = width_2 / width_1 ;
  
  // total count_vector_2
  int n = count_vector_1.length() + width_1 - 1 - width_2 + 1 ;
  
  // initialise count_vector_2
  NumericVector count_vector_2( n ) ;
  
  // count
  for( int i=0; i<width_2; i++) {
    
    // find non-overlap bin as if binning method
    int cnt = 0 ;
    NumericVector bin_1( (count_vector_1.length() + width_1 - 1 - i) / width_1 ) ;
    for (int idx = i; idx < count_vector_1.length(); idx += width_1) {
      bin_1[cnt] = idx ;
      cnt++ ;
    }

    // last bin for vector 2
    int bin_2_idx_end = (bin_1.length() / scale * scale) - 1 ;

    // absorb suitable bin1 to bin2
    NumericVector bin_2_idx = bin_1[ Range(0, bin_2_idx_end) ] ;
    
    // absorb the count
    NumericVector bin_2 = count_vector_1[ bin_2_idx ] ;
    
    // merge bin1
    cnt = 0 ;
    NumericVector bin_2_merge( bin_2.length() / scale  ) ;
    for (int j = 0; j < bin_2.length(); j += scale) {
      
      bin_2_merge[cnt] = sum( bin_2[ Range(j, j+scale-1)] ) ;
      cnt++ ;
    }
    
    // position the bin appropriately to reflect sliding (s) window
    cnt = 0 ;
    for (int s_idx = i; s_idx < count_vector_2.length(); s_idx += width_2) {
      count_vector_2[s_idx] = bin_2_merge[cnt] ;
      cnt++ ;
    }
    
  }
  
  return count_vector_2 ;
}