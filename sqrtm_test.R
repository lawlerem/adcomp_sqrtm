library(TMB)
compile("sqrtm_test.cpp")
dyn.load(dynlib("sqrtm_test"))

n<- 5
mat<- rWishart(1, n, diag(n))[, , 1]

obj<- MakeADFun(
  data = list(
    idx = c(0, 0)
  ),
  para = list(
    mat = mat
  ),
  DLL = "sqrtm_test"
)

cpp_sqrt<- obj$report()$mat_sqrt
cpp_grad<- function(i, j) {
  obj$env$data$idx[c(1, 2)]<- c(as.integer(i - 1), as.integer(j - 1))
  obj$retape()
  return( c(obj$gr()) )
}

f<- function(M) return( expm::sqrtm(M) )
r_sqrt<- f(mat)
r_grad<- function(i, j) {
  func<- function(x, i, j) {
    xmat<- matrix(x, nrow = nrow(mat), ncol = ncol(mat))
    return( f(xmat)[i, j])
  }
  return( numDeriv::grad(func, x = c(mat), i = i, j = j) )
}

max_error<- mat
for( i in seq(n) ) {
  for( j in seq(n) ) {
    max_error[i, j]<- max( r_grad(i, j) - cpp_grad(i, j) )
  }
}

# Looks good, maximum error for derivatives is ~1e-10
max_error



# Check standardization using sqrtm
library(MASS)
X<- mvrnorm(
  n = 1000,
  mu = numeric(n),
  Sigma = mat
)
Zsqrt<- X %*% solve(t(cpp_sqrt))

# Should be ~identity matrix
cov(Zsqrt)