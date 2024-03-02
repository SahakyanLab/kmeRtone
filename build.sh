rm ../kmeRtone_0.1.tar.gz
R --slave -e 'Sys.setenv("PKG_CXXFLAGS"="-fopenmp -std=c++11");library(devtools);document();build();install()'
 #-I/project/sahakyanlab/adibabdu/Programs/conda/lib/R/library/stringi/include
# devtools::install_git("git@github.com:SahakyanLab/kmeRtone.git", ref = "rpackage")
