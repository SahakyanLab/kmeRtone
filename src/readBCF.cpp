#include <Rcpp.h>
#include <htslib/synced_bcf_reader.h>

//' Read BCF files and extract chromosome information
//'
//' This function reads BCF (Binary Variant Call Format) files and extracts chromosome information.
//' @param fname The file name or URL of the BCF file.
//' @param regions A string specifying the regions to be read from the BCF file.
//' @param is_file Boolean indicating if the regions parameter is a file.
//' @return A list containing chromosome information.
Rcpp::List readBCF(const std::string& fname, const std::string& regions, int is_file) {
  
  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
  bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

  if (bcf_sr_set_regions(sr, regions.c_str(), is_file) == -1) {
    Rcpp::stop("Failed to set the regions.");
  }

  if (!bcf_sr_add_reader(sr, fname.c_str())) {
    if (sr->errnum)
      Rcpp::stop("Error adding reader: %s\n", bcf_sr_strerror(sr->errnum));
    else
      Rcpp::stop("Unknown error adding reader.");
  }

  std::vector<int32_t> chrom;

  while (bcf_sr_next_line(sr)) {
    bcf1_t *line = bcf_sr_get_line(sr, 0);
    chrom.push_back(line->rid);
  }
  
  if (sr->errnum) {
    Rcpp::stop("Error: %s\n", bcf_sr_strerror(sr->errnum));
  }

  bcf_sr_destroy(sr);
  
  return Rcpp::List::create(Rcpp::Named("CHROM") = chrom);
}