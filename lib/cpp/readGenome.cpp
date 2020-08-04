#include <fstream>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

// Fast retrieving DNA sequence from a fasta file with one single-line header.
// User can select a sequence as speficied by start and end coordinate and output
// it in a string or vector format. Besides that, the user can get the size of the
// entire fasta sequence.


// [[Rcpp::export]]
CharacterVector readGenome(std::string chromosome, std::string genome_path, std::string genome_prefix="",
                           std::string suffix=".fa.gz", std::string form="vector",
                           std::string letter_case="upper", std::string FULL_PATH="") {
  
  std::string fname ;
  // get chromosome path
  if (FULL_PATH == "") {
    fname = genome_path + prefix + chromosome + suffix ;
  } else {
    fname = FULL_PATH ;
  }
  
  // open file
  std::ifstream fasta ;
  fasta.open(fname, std::ios::binary) ;
  if(!fasta) Rcout << ("Can't open file for reading.");
  
  // ignore first line
  // in.ignore(100000,'\n') ;
  
  // get 1st line to get number of character to skip using seekg
  std::string line_1 ;
  std::getline(fasta, line_1) ;
  
  // get 2nd line to know the number of bases per line
  std::string line_2 ;
  std::getline(fasta, line_2) ;
  
  // get length of file to calculate total length of chromosome
  fasta.seekg (0, std::ios::end) ;
  int total_size = fasta.tellg() ;
  
  // check last line
  std::string line_fin ;
  std::getline(fasta, line_fin) ;
  
  int cnt = -1 ;
  while (line_fin.length() == 0) {
    fasta.clear() ;
    fasta.seekg (cnt, std::ios::end) ;
    total_size = fasta.tellg() ;
    getline(fasta, line_fin) ;
    cnt-- ;
  }
  
  // calculate total character
  total_size = total_size - line_1.length() + 1 ;
  
  // total "\n"
  int total_newline = total_size / (line_2.length() + 1) + ((total_size % line_2.length()+1) > 0) ;
  
  // offset the newline character
  total_size = total_size - total_newline ;
    
  if (size == true) {
    return std::to_string(total_size) ;
  }
  
  // if size is out of range stop
  if (end > total_size | start < 1) {
    stop("Coordinate is out of range.") ;
  } else if (end < start) {
    stop("End coordinate cannot be less than stop coordinate.") ;
  }
  
  fasta.clear() ;
  
  // new start and end index; -1 because c++ zero-index; include the newline character
  int start_ = start - 1 + (start / line_2.length() - 1) + (start % line_2.length() > 0) ;
  int end_ = end -1 + (end / line_2.length() - 1) + (end % line_2.length() > 0) ;

  // calculate start position of the requested DNA (include the newline character)
  int start_position = line_1.length() + 1 + start_;
  
  // move pointer to the start position
  fasta.seekg(start_position) ; 
  
  // resize to the sequence length
  std::string contents ;
  contents.resize( end_ - start_ + 1 ) ; 
  
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