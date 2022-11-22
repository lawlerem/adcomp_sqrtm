A matrix square root function for use in TMB C++ template functions. See https://github.com/kaskr/example_cpp, https://github.com/kaskr/adcomp/pull/370, and https://github.com/kaskr/adcomp/.

This is an atomic version of the function that is compatible with automatic differentiation.
This function is useful when standardizing a general multivariate normal vector via the transformation **Z** = **A**<sup>-1/2</sup> * (**X** - **mu**), where **A** is the covariance matrix of the original vector **X**.

To use this function, do the following:

This installs the contributed code:
```R
TMB:::install.contrib("https://github.com/lawlerem/adcomp_sqrtm/archive/master.zip")
```

After that you can use it in a TMB model by adding the following two lines to your C++ file:
```C++
#include <TMB.hpp>
#include <adcomp_sqrtm-main/sqrtm.hpp>
using contrib::lawlerem::sqrtm;
```

A test of the implementation is included in sqrtm_test.R and sqrtm_test.cpp, and a toy example using the function is given in sqrtm_example.R and sqrtm_example.cpp.
The example is a simple longitudinal model, where the sqrtm function is used to transform the dependent time series of subject into a sequence of standard normal variables.
Note that the example could be estimated much more compactly using TMB's MVNORM_t objects, but the sqrtm version allows you to calculate the quantiles of each observation which could then, for example, be fed into a copula model to model the dependence between subjects.

As currently implemented, sqrtm should only be used on small matrices since it is computationally very inefficient (see https://github.com/kaskr/adcomp/pull/370#issuecomment-1317603340).