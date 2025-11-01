##==========================================================================
## Guidelines:

# Required quantities in general:
#   1. gamma.11: Response rate in the first stage
#   2. nu.d1: Skewness of the data
#   3. alpha: Size of the test. Taken 0.05 here.

# Required for precision-based formulas:
#   1. psi.d1: Precision factor

# Required for power-based formulas:
#   1. eff.size: Standardized effect sizes
#   2. pow: Desired power of the test. Takes 0.8 here.

# Required for formulas with Setting 2:
#   1. gamma.12: Response rate in the second stage
#   2. r1 and r2: Scaling factors of sigma and kappa in DTR3 
#                 with respect to DTR1. Both taken as 1 here.

##==========================================================================
## Precision-based (Setting 1)

samp.size.prec.s1 <- function(gam.11, psi.d1, nu.d1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- (4 - 2*gam.11) * qnorm(1-alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ceiling(ret))
}


##==========================================================================
## Precision-based (Setting 2)

samp.size.prec.s2 <- function(gam.11, gam.12, psi.d1, nu.d1, r1 = 1, r2 = 1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  
  ret <- ( (4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12) ) * qnorm(1-alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ceiling(ret))
}

##==========================================================================
## Power-based (Setting 2)

samp.size.pow.s2 <- function(gam.11, gam.12, nu.d1, eff.size,
                             pow = 0.8, alpha = 0.05, r1 = 1, r2 = 1)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- 2 * (qnorm(pow) + qnorm(1- alpha/2))^2 * ((4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12)) * kap.d1 / (eff.size^2 * (1+r1))
  return(ret)
}
