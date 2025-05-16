#==========================================================================
#Overall Sample sizes

rm(list = ls())

#------------------------------
#the function needed

samp.size.pow.s2 <- function(gam.11 = 0.35, gam.12 = 0.35, nu.d1 = 1, eff.size = 0.2,
                              pow = 0.8, alpha = 0.05, r1 = 1, r2 = 1)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- 2 * (qnorm(pow) + qnorm(1- alpha/2))^2 * ((4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12)) * kap.d1 / (eff.size^2 * (1+r1))
  return(ret)
}


#------------------------------
#tabulating values

eff.vec <- c(0.1, 0.15, 0.2, 0.25, 0.3)
gam.11.vec <- c(0.35, 0.5, 0.65)
gam.12.vec <- c(0.35, 0.5, 0.65)
nu.d1.vec <- c(-100, -50, -25, -10, -8, -5, -4, -3, -2.5, -2, -1.5, -1, -0.75, -0.5, -0.25,  0,
               0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 25, 50, 100)

store.n <- array(dim = c(length(eff.vec),length(gam.11.vec),length(nu.d1.vec)))

for (i.eff in 1:length(eff.vec))
{
  for (i.gam in 1:length(gam.11.vec))
  {
    for (i.nu.d1 in 1:length(nu.d1.vec))
    {
      gam.11 <- gam.11.vec[i.gam]
      gam.12 <- gam.12.vec[i.gam]
      nu.d1 <- nu.d1.vec[i.nu.d1]
      eff <- eff.vec[i.eff]
      
      #-- - -- --
      
      foo <- samp.size.pow.s2(gam.11, gam.12, nu.d1, eff)
      
      store.n[i.eff, i.gam, i.nu.d1] <- foo
      
    }
  }
}



#==========================================================================
#Operation Characteristics and the plots

library(ggplot2)
library(latex2exp)

rm(list = ls())

#------------------------------
#the function needed (with non-integer output)

samp.size.pow.s2 <- function(gam.11 = 0.35, gam.12 = 0.35, nu.d1 = 1, eff.size = 0.2,
                              pow = 0.8, alpha = 0.05, r1 = 1, r2 = 1)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  ret <- 2 * (qnorm(pow) + qnorm(1- alpha/2))^2 * ((4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12)) * kap.d1 / (eff.size^2 * (1+r1))
  return(ret)
}

#------------------------------
#Dependencies over psi

eff.size <- seq(0.10,0.4,0.01)

vals.eff.size <- samp.size.pow.s2(eff.size = eff.size)

ggplot() +
  geom_line(mapping = aes(x = eff.size, y= vals.eff.size), lwd = 1) +
  labs(x = TeX(r"($\Delta$)"), y = TeX(r"($N$ )")) +
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 21))

#------------------------------
#Dependencies over response

gam.1 <- seq(0.05,0.95,0.05)

vals.gam.1 <- samp.size.pow.s2(gam.11 = gam.1, gam.12 = gam.1)

ggplot() +
  geom_line(mapping = aes(x = gam.1, y= vals.gam.1), lwd = 1) +
  labs(x = TeX(r"($\gamma_{11}$)"), y = TeX(r"($N$ )")) +
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 21))

#------------------------------
#Dependencies over skewness

nu.d1 <- c(-10, -8, -5, -4, -3, -2.5, -2, -1.5, -1, -0.75, -0.5, -0.25,  0,
           0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10)

vals.nu.d1 <- samp.size.pow.s2(nu.d1 = nu.d1)

ggplot() +
  geom_line(mapping = aes(x = nu.d1, y= vals.nu.d1), lwd = 1) +
  labs(x = TeX(r"($\nu_{d_1}$)"), y = TeX(r"($N$ )")) +
  theme(axis.text = element_text(size = 17),
        axis.title = element_text(size = 19))

