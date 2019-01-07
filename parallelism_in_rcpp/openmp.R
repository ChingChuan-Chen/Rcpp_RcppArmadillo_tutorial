# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)

# cpp file
Rcpp::sourceCpp("example_openmp.cpp")

system.time(long_computation(1e4, 1L))
# Number of threads=1
#   user  system elapsed 
#   11.04    0.02   11.22 
system.time(long_computation(1e4, 4L))
# Number of threads=4
#   user  system elapsed 
#   14.22    0.03    3.75
