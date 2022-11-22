library(TMB)
compile("sqrtm_example.cpp")
dyn.load(dynlib("sqrtm_example"))

n<- 20
k<- 5

sd<- 2
ar1<- 0.6

simobj<- MakeADFun(
  data = list(
    X = matrix(0, nrow = n, ncol = k)
  ),
  para = list(
    log_sd = log(sd),
    logit_ar1 = qlogis(0.5 * (ar1 + 1))
  ),
  DLL = "sqrtm_example"
)
sim<- simobj$simulate()


obj<- MakeADFun(
  data = list(
    X = sim$X
  ),
  para = list(
    log_sd = log(1),
    logit_ar1 = qlogis(0.5 * (0 + 1))
  ),
  DLL = "sqrtm_example"
)
opt<- nlminb(obj$par, obj$fn, obj$gr)
sdr<- sdreport(obj, opt$par)

# Parameter estimates are good!
summary(sdr)