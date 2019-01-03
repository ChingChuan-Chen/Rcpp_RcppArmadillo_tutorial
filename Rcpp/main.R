# hello
Rcpp::sourceCpp("hello.cpp")

# rnorm
Rcpp::sourceCpp("rnorm.cpp")

set.seed(100)
cpp_rnorm(5)
set.seed(100)
cpp_rnorm(5)

# multiplication
Rcpp::sourceCpp("multiplication.cpp")

# sapply
Rcpp::sourceCpp("cpp_function.cpp")

# rnorm with options
Rcpp::sourceCpp("rnorm_with_opts.cpp")

# modify_df
Rcpp::sourceCpp("modify_df.cpp")

