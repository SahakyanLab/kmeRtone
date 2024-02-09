#include <Rcpp.h>
#include <stringi.h>
#include <random>
#include "countRangedKmers.h"

//' Simulate a population given ranges of chromosome sequence to mutate.
//'
//' This is internal function - not exported.
//'
//' @param chrom_seq A chromosome sequence.
//' @param start Start positions.
//' @param end End positions.
//' @param strand Strand type: "+" or "-".
//' @param snv_df A table of SNV frequency. Columns: position, base, count.
//' @param pop_size Size of population.
//' @param top_kmers Extreme k-mers i.e. highly susceptible k-mers.
//' @param central_pattern K-mer central pattern.
//' @param k K-mer size.
//' @return A count matrix with 4 rows for total top k-mers and susceptible
//'    k-mers in sense and antisense. Columns correspond to population
//'    individuals.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix simulatePopulation(Rcpp::CharacterVector chrom_seq,
    std::vector<int> starts, std::vector<int> ends,
    std::string strand, Rcpp::DataFrame snv_df, int pop_size,
    Rcpp::CharacterVector top_kmers, Rcpp::CharacterVector central_pattern,
    Rcpp::NumericVector k) {
  
  // 4: top_kmer_count sense, top_kmer_count antisense,
  //    susc_kmer_count antisense, susc_kmer_count antisense
  Rcpp::NumericMatrix counts(4, pop_size);
  
  Rcpp::NumericVector snv_df_position = snv_df["position"];
  Rcpp::CharacterVector snv_df_base = snv_df["base"];
  Rcpp::NumericVector snv_df_count = snv_df["count"];
  Rcpp::NumericVector snv_pos = Rcpp::unique(snv_df_position);
  std::sort(snv_pos.begin(), snv_pos.end());
  int snv_tot = snv_pos.length();
  Rcpp::NumericVector cp_len = stri_length(central_pattern);
  Rcpp::NumericVector flank = (k - cp_len) / 2;
  
  // Reverse complement top_kmers and central_pattern
  Rcpp::CharacterVector top_kmers_rc = stri_trans_char(top_kmers,
                                      Rcpp::CharacterVector::create("ACGT"),
                                      Rcpp::CharacterVector::create("TGCA"));

  top_kmers_rc = stri_reverse(top_kmers_rc);
  
  Rcpp::CharacterVector central_pattern_rc = stri_trans_char(central_pattern,
                                            Rcpp::CharacterVector::create("ACGT"),
                                            Rcpp::CharacterVector::create("TGCA"));

  central_pattern_rc = stri_reverse(central_pattern_rc);

  // Start simulation
  for (int i = 0; i < pop_size; ++i) {
    
    // Sample SNV based on frequency
    // The loop here is slow. Cannot use openmp as the Rcpp-wrapped R object is
    // not thread safe. data.table groupby is faster than this!
    Rcpp::CharacterVector bases(snv_tot);
    for (int j = 0; j < snv_tot; ++j) {
      
      Rcpp::LogicalVector idx = snv_df_position == snv_pos[j];
      Rcpp::NumericVector base_counts = snv_df_count[idx];
      // Rcpp::NumericVector base_probs = base_counts / std::accumulate(base_counts.begin(), base_counts.end(), 0);
      Rcpp::CharacterVector base_names = snv_df_base[idx];

      //   sample(const Vector<RTYPE>& x, int size, bool replace = false, sugar::probs_t probs = R_NilValue)
      // Rcpp::CharacterVector b = Rcpp::sample(base_names, 1, false, base_probs);
      std::default_random_engine generator;
      std::discrete_distribution<int> dd{base_counts.begin(), base_counts.end()};
      bases[j] = base_names[dd(generator)];
    }

    // Mutate genome
    Rcpp::List list_snv_pos = Rcpp::List::create(snv_pos);
    Rcpp::CharacterVector mut_chrom = stri_sub_replacement_all(
      chrom_seq, list_snv_pos, R_NilValue, Rcpp::List::create(1),
      Rcpp::LogicalVector(1, false), Rcpp::List::create(bases),
      Rcpp::LogicalVector(1, true));
    //SEXP from, SEXP to, SEXP length, SEXP omit_na, SEXP value, SEXP use_matrix
    
    // Count k-mers
    Rcpp::NumericVector kmers = Rcpp::wrap(countRangedKmers(Rcpp::as<std::string>(mut_chrom),
                                             starts, ends, Rcpp::as<int>(k)));

    // Take relevant kmers in reference strand
    Rcpp::NumericVector top_cnts_ref = kmers[top_kmers];
    double top_cnt_ref = std::accumulate(top_cnts_ref.begin(), top_cnts_ref.end(), 0);
    // stri_sub(SEXP str, SEXP from, SEXP to, SEXP length, SEXP use_matrix,
    //          SEXP ignore_negative_length)
    Rcpp::CharacterVector cp_kmers = stri_sub(kmers.names(), flank, R_NilValue,
                             stri_length(central_pattern),
                             Rcpp::LogicalVector(1, true),
                             Rcpp::LogicalVector(1, false));
    Rcpp::LogicalVector idx_susc_ref = cp_kmers == central_pattern;
    Rcpp::NumericVector susc_cnts_ref = kmers[idx_susc_ref];
    double susc_cnt_ref = std::accumulate(susc_cnts_ref.begin(), susc_cnts_ref.end(), 0);

    // Take relevant kmers in opposite strand
    Rcpp::NumericVector top_cnts_opp = kmers[top_kmers_rc];
    double top_cnt_opp = std::accumulate(top_cnts_opp.begin(), top_cnts_opp.end(), 0);
    Rcpp::LogicalVector idx_susc_opp = cp_kmers == central_pattern_rc;
    Rcpp::NumericVector susc_cnts_opp = kmers[idx_susc_opp];
    double susc_cnt_opp = std::accumulate(susc_cnts_opp.begin(), susc_cnts_opp.end(), 0);

    // Update the count
    if (strand == "+") {
      
      counts(0,i) = counts(0,i) + top_cnt_ref;
      counts(1,i) =  counts(1,i) + top_cnt_opp;
      counts(2,i) = counts(2,i) + susc_cnt_ref;
      counts(3,i) = counts(3,i) + susc_cnt_opp;

    } else if (strand == "-") {
      
      counts(0,i) = counts(0,i) + top_cnt_opp;
      counts(1,i) =  counts(1,i) + top_cnt_ref;
      counts(2,i) = counts(2,i) + susc_cnt_opp;
      counts(3,i) = counts(3,i) + susc_cnt_ref;

    }
    std::cout << i << std::endl;
  }
  
  return counts;
}
