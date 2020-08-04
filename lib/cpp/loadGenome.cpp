#include <fstream>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

// Fast retrieving DNA sequence from a fasta file with one single-line header.
// User can select a sequence as speficied by start and end coordinate and output
// it in a string or vector format. Besides that, the user can get the size of the
// entire fasta sequence.


// [[Rcpp::export]]
CharacterVector loadGenomeCPP(std::string chromosome, std::string genome_path, std::string genome_prefix="",
                           std::string genome_suffix=".fa.gz", std::string form="vector",
                           std::string letter_case="upper", std::string FULL_PATH="") {
  
  std::string fname ;
  // get chromosome path
  if (FULL_PATH == "") {
    fname = genome_path + genome_prefix + chromosome + genome_suffix ;
  } else {
    fname = FULL_PATH ;
  }
  
  // open file
  std::ifstream fasta ;
  fasta.open(fname, std::ios::binary) ;
  if(!fasta) Rcout << ("Can't open file for reading.");
  
  // get 1st line to get number of character to skip using seekg
  std::string line_1 ;
  std::getline(fasta, line_1) ;
  
  // get length of file to calculate total length of chromosome
  fasta.seekg (0, std::ios::end) ;
  int total_size = fasta.tellg() ;
  
  fasta.clear() ;
  
  // move pointer to the start position (add more +1 to include \n)
  fasta.seekg(line_1.length() + 2) ; 
  
  // resize to the sequence length
  std::string contents ;
  contents.resize( total_size - line_1.length() + 2 + 1 ) ; 
  
  // get the sequence
  fasta.read((&contents[0]), contents.size()) ;
  
  // close the file
  fasta.close() ;
  
  contents.erase(std::remove(contents.begin(), contents.end(), '\n'), contents.end());
  
  if (form == "vector") {
    // convert std::string to std::vector<char> to CharacterVector
    std::vector<char> contents_vector(contents.begin(), contents.end()) ;
    CharacterVector contents_vector_ = Rcpp::wrap(contents_vector);
    return contents_vector_ ;
  }
  
  return contents ;
}