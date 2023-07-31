#include <fstream>
#include <Rcpp.h>
using namespace Rcpp;

#include <zlib.h>

//' Read fast single-header FASTA file.
//'
//' Slower than data.table::fread
//' 
//' @param file_path A path to FASTA file.
//' @param unmask Capitalize all base letters?
//' @return A single sequence string without header.
//'
//' @export
// [[Rcpp::export]]
void splitFASTA2(std::string file_path, std::string output_dir) {
  
  std::ofstream outfile;
  gzFile infile;
  infile = gzopen(file_path.c_str(), "rb");
  
  std::string chunk;
  std::string line;
  std::string header;
  std::string fa_name;
  
  if (!infile) {
    stop("Failed to open FASTA file " + file_path);
  }
  
  char outbuffer[1024*16];
  
  while(!gzeof(infile)) {

    // Read by chunk
    int len = gzread(infile, outbuffer, sizeof(outbuffer));
    chunk.insert(chunk.end(), outbuffer, outbuffer+len);
    
    // Count number of newline
    int cnt = std::count(chunk.begin(), chunk.end(), '\n');

    // Read chunk as stream file
    std::stringstream ss(chunk);
    for (int i = 1; i <= cnt; ++i) {
      std::getline(ss, line, '\n');
      
      // Check for header
      if (line.at(0) == '>') {
        header = line;
        fa_name = header.substr(1, header.find(" ")-1) + ".fa";
        
        // Open connection
        // gzclose(outfile);
        // outfile = gzopen((output_dir + "/" + fa_name).c_str(), "ab");
        outfile.close();
        outfile = std::ofstream((output_dir + "/" + fa_name).c_str(),
                                std::ios::app);
      }
      
      // Save file
      outfile << line << std::endl;
      // gzwrite(outfile, line.c_str(), strlen(line.c_str()));
      
    }
    
    // Save the last line without \n character for the next loop
    if (std::getline(ss, line)) {
      chunk = line;
    } else {
      chunk.clear();
    }
  }
  
  gzclose(infile);
  outfile.close();

  return;
}