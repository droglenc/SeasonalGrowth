################################################################################
## Main Function
##   Linf, t0 as usual
##   Kpr = K-prime as defined in Pauly et al. (units are diff than usual K)
##   ts = start of sinusoidal growth (maximum growth rate)
##   NGT = "No Growth Time" = "fraction of a year where no growth occurs"
##   tpr = "t-prime" = actual age (t) minus cumulative NGT prior to t
################################################################################

VBSCGF <- function(t,Linf,Kpr=NULL,t0=NULL,ts=NULL,NGT=NULL) {
  if (length(Linf)==5) { Kpr <- Linf[[2]]; t0 <- Linf[[3]]
  ts <- Linf[[4]]; NGT <- Linf[[5]]
  Linf <- Linf[[1]] }
  tpr <- iCalc_tpr(t,ts,NGT)
  q <- Kpr*(tpr-t0) +
    (Kpr*(1-NGT)/(2*pi))*sin((2*pi)/(1-NGT)*(tpr-ts)) -
    (Kpr*(1-NGT)/(2*pi))*sin((2*pi)/(1-NGT)*(t0-ts))
  Linf*(1-exp(-q))
}
