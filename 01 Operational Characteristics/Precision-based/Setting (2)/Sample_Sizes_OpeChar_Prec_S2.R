#==========================================================================
#Overall Sample sizes

rm(list = ls())

#------------------------------
#the function needed

samp.size.prec.s2 <- function(gam.11, gam.12, psi.d1, nu.d1, r1 = 1, r2 = 1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  
  ret <- ( (4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12) ) * qnorm(1-alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ceiling(ret))
}


#------------------------------
#tabulating values

psi.vec <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
gam.11.vec <- c(0.35, 0.5, 0.65)
gam.12.vec <- c(0.35, 0.5, 0.65)
nu.d1.vec <- c(-100, -50, -25, -10, -8, -5, -4, -3, -2.5, -2, -1.5, -1, -0.75, -0.5, -0.25,  0,
               0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 25, 50, 100)

store.n <- array(dim = c(length(psi.vec),length(gam.11.vec),length(nu.d1.vec)))

for (i.psi in 1:length(psi.vec))
{
  for (i.gam in 1:length(gam.11.vec))
  {
    for (i.nu.d1 in 1:length(nu.d1.vec))
    {
      gam.11 <- gam.11.vec[i.gam]
      gam.12 <- gam.12.vec[i.gam]
      nu.d1 <- nu.d1.vec[i.nu.d1]
      psi <- psi.vec[i.psi]
      
      #-- - -- --
      
      foo <- samp.size.prec.s2(gam.11, gam.12, psi, nu.d1)
      
      store.n[i.psi, i.gam, i.nu.d1] <- foo
      
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

samp.size.prec.s2 <- function(gam.11 = 0.35, gam.12 = 0.35, psi.d1 = 0.4, nu.d1 = 5, r1 = 1, r2 = 1, alpha = 0.05)
{
  xi.d1 <- nu.d1 / sqrt(1+nu.d1^2)
  kap.d1 <- (1 - 2*xi.d1^2 / pi) 
  
  ret <- ( (4 - 2*gam.11) + r1 * r2 * (4 - 2*gam.12) ) * qnorm(alpha/2)^2 / (psi.d1^2 * kap.d1) 
  
  return(ret)
}

#------------------------------
#Dependencies over psi

psi <- seq(0.20,0.9,0.01)

vals.psi <- samp.size.prec.s2(psi = psi)

ggplot() +
  geom_line(mapping = aes(x = psi, y= vals.psi), lwd = 1) +
  labs(x = TeX(r"($\psi$)"), y = TeX(r"($N$ )")) +
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 21)) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), limits = c(0,1))

#------------------------------
#Dependencies over response

gam.11 <- seq(0.05,0.95,0.05)

vals.gam.11 <- samp.size.prec.s2(gam.11 = gam.11)

ggplot() +
  geom_line(mapping = aes(x = gam.11, y= vals.gam.11), lwd = 1) +
  labs(x = TeX(r"($\gamma_{11}$)"), y = TeX(r"($N$ )")) +
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 21)) +
  scale_y_continuous(breaks = c(350,390,430)) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), limits = c(0,1))

#------------------------------
#Dependencies over skewness

nu.d1 <- c(-15, -10, -8, -5, -4, -3, -2.5, -2, -1.5, -1, -0.75, -0.5, -0.25,  0,
           0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 15)

vals.nu.d1 <- samp.size.prec.s2(nu.d1 = nu.d1)

ggplot() +
  geom_line(mapping = aes(x = nu.d1, y= vals.nu.d1), lwd = 1) +
  labs(x = TeX(r"($\nu_{d_1}$)"), y = TeX(r"($N$ )")) +
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 21)) + 
  scale_y_continuous(breaks = c(150,240,330,420))
