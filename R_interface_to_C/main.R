# 1. hello
## build
system("build.bat 1-hello.c")

## load dll
dyn.load("1-hello.dll")

## check whether function is loaded
is.loaded("hello") # TRUE

## run
.C("hello", 5L)

## unload dll
## build
dyn.unload("1-hello.dll")

# 2. get normal random numbers
system("build.bat 2-rnorm.c")

## load dll
dyn.load("2-rnorm.dll")

## check whether function is loaded
is.loaded("C_rnorm") # TRUE

## run
### output X
x <- vector("double", 5L)
### run C
result <- .C("C_rnorm", length(x), 1, 2, x)
### print result
print(result[[4]])

## test set.seed
set.seed(123)
.C("C_rnorm", length(x), 1, 2, x)[[4]]
set.seed(123)
.C("C_rnorm", length(x), 1, 2, x)[[4]]

## unload dll
dyn.unload("2-rnorm.dll")
