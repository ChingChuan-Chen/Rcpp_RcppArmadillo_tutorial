##### 1. hello #####
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

##### 2. get normal random numbers #####
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


##### 3. hello #####
## build
system("build.bat 3-hello-again.c")

## load dll
dyn.load("3-hello-again.dll")

## check whether function is loaded
is.loaded("hello_again") # TRUE

## run
.Call("hello_again", 5L)

## unload dll
dyn.unload("3-hello-again.dll")


##### 4. get normal random numbers #####
## build
system("build.bat 4-rnorm2.c")

## load dll
dyn.load("4-rnorm2.dll")

## check whether function is loaded
is.loaded("rnorm2") # TRUE

## run
.Call("rnorm2", 5L, 1, 2)

## set seed
set.seed(100)
.Call("rnorm2", 5L, 1, 2)
set.seed(100)
.Call("rnorm2", 5L, 1, 2)

## unload dll
dyn.unload("4-rnorm2.dll")
