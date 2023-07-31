#include <Rcpp.h>
using namespace Rcpp;

#include <zlib.h>
#include <regex>

// [[Rcpp::plugins(cpp11)]]
//' Read fast single-header FASTA file.
//'
//' @param file_path A path to FASTA file.
//' @param unmask Capitalize all base letters?
//' @return A single sequence string without header.
//'
//' @export
// [[Rcpp::export]]
std::string readSingleFASTA(std::string file_path, std::string mask = "none") {

  gzFile infile;
  infile = gzopen(file_path.c_str(), "rb");

  std::string seq;

  if (!infile) {
    stop("Failed to open FASTA file " + file_path);
  }

  char outbuffer[1024*16];

  while(!gzeof(infile)) {
    int len = gzread(infile, outbuffer, sizeof(outbuffer));
    seq.insert(seq.end(), outbuffer, outbuffer+len);
  }

  // Remove header on the first line
  seq.erase(0, seq.find('\n') + 1);

  // Remove newline character
  seq.erase(std::remove(seq.begin(), seq.end(), '\n'),
      seq.end());

  // Unmask the soft mask i.e. capitalize base letters
  if (mask == "none") {
    std::transform(seq.begin(), seq.end(), seq.begin(),
        (int(*)(int))( std::toupper ));
  } else if (mask == "hard") {
    std::regex rgx ("[a-z]");
    seq = std::regex_replace(seq, rgx, "N");
  }

  gzclose(infile);

  return seq;
}

