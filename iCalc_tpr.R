################################################################################
## internal function to compute t-prime
################################################################################
iCalc_tpr <- function(t,WP,NGT) {
  ## Step 1
  SNG <- WP-NGT/2
  tmp.t <- t-SNG
  ## Step 2 (in parentheses) and Step 3
  tmp.t2 <- (tmp.t-floor(tmp.t)) - NGT
  ## Step 4
  tmp.t2[tmp.t2<0] <- 0
  ## Step 5 (in parentheses) and Step 6 (also returns value)
  (floor(tmp.t)*(1-NGT)+tmp.t2) + SNG
}

