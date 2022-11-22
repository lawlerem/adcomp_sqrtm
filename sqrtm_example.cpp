#include <TMB.hpp>
#include <adcomp_sqrtm-main/sqrtm.hpp>
using contrib::lawlerem::sqrtm;


template<class Type>
Type objective_function<Type>::operator() () {
  static const Type pi = 3.14159;
  DATA_MATRIX(X); // n x k matrix

  PARAMETER(log_sd);
  PARAMETER(logit_ar1);

  int n = X.rows();
  int k = X.cols();
  matrix<Type> Z(X);
  matrix<Type> Q(X);

  Type sd = exp(log_sd); // (-Inf, +Inf) --> (0, +Inf)
  Type ar1 = 2 * invlogit(logit_ar1) - 1; // (-Inf, +Inf) --> (-1, +1)

  matrix<Type> Sigma(X.cols(), X.cols());
  for(int i = 0; i < k; i++) {
    for(int j = 0; j < k; j++) {
      Sigma(i, j) = pow(sd, 2) * pow(ar1, abs(i - j));
    }
  }
  matrix<Type> sqrtSigma = sqrtm(Sigma);
  matrix<Type> sqrtPrecision = atomic::matinv(sqrtSigma);

  Type nll = 0.0;
  for(int i = 0; i < n; i++) {
    Z.row(i) = sqrtPrecision * vector<Type>(X.row(i));
    for(int j = 0; j < k; j++) {
      Q(i, j) = pnorm(Z(i, j), Type(0.0), Type(1.0));
    }
  }

  // ll = sum_Z (1 / 2) * ( -k * 2 * pi - logdet(Sigma) - Z^T * Z)
  //    = (n / 2) * ( -2 * k * pi - logdet(Sigma)) - (1 / 2) * sum_Z Z^T * Z

  nll -= Type(n) * 0.5 * ( -2.0 * Type(k) * pi - atomic::logdet(Sigma));
  for(int i = 0; i < n; i++) {
    nll -= -0.5 * Type(Z.row(i).squaredNorm());
  }

  SIMULATE{
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < k; j++) {
        Q(i, j) = runif(0.0, 1.0);
        Z(i, j) = qnorm(Q(i, j), Type(0.0), Type(1.0));
      }
      X.row(i) = sqrtSigma * vector<Type>(Z.row(i));
    }
    REPORT(Q);
    REPORT(Z);
    REPORT(X);
  }

  REPORT(X);
  REPORT(Z);
  REPORT(Q);
  ADREPORT(sd);
  ADREPORT(ar1);

  return nll;
}
