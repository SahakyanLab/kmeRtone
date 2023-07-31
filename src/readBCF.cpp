#include <Rcpp.h>
#include <htslib/synced_bcf_reader.h>

//' @export
// [[Rcpp::export]]
Rcpp::List readBCF(std::string fname, std::string regions, int is_file) {
  
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
      Rcpp::stop("Unkown error adding reader.");
  }
  
  // To be continued on weekend...
  // Log running this in R.
  // a <- readBCF(fname = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr10.vcf.bgz", regions = "chr10:31340490-31340499", is_file = 0L)
  // 1) Failed to read index file. I need to download manually the index file, then
  //    it runs.
  // 2) Failed to get chrom number. It did not even go inside the loop. So far this
  //    function only return chrom to test.

  std::vector<int32_t> chrom;
  std::cout << 7 << std::endl;
  std::cout << bcf_sr_has_line(sr, 0) << std::endl;
  while ( bcf_sr_next_line(sr) )
  {
    
    std::cout << 8 << std::endl;
    bcf1_t *line = bcf_sr_get_line(sr, 0);
    std::cout << 9 << std::endl;
    std::cout << line->rid << std::endl;
    std::cout << 10 << std::endl;
    chrom.push_back(line->rid);
    std::cout << 11 << std::endl;
    
  }
  
  std::cout << 1 << std::endl;
  
  if ( sr->errnum ) Rcpp::stop("Error: %s\n", bcf_sr_strerror(sr->errnum));
  bcf_sr_destroy(sr);
  
  return Rcpp::List::create(Rcpp::Named("CHROM") = chrom);
}