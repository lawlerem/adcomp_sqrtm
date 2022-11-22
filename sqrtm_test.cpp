#include <TMB.hpp>
#include <adcomp_sqrtm-main/sqrtm.hpp>
using contrib::lawlerem::sqrtm;

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_IVECTOR(idx); // (i, j)
  PARAMETER_MATRIX(mat);

  matrix<Type> mat_sqrt = sqrtm(mat);

  REPORT(mat);
  REPORT(mat_sqrt);

  return mat_sqrt(idx(0), idx(1));
}
