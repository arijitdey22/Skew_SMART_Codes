library(ggplot2)
library(latex2exp)

SN_vis <- function(mu, sigma, nu){
  X <- seq(mu - 3 * sigma, mu + 3 * sigma, 0.001) 
  Y.axis <- 2 / sigma * dnorm((X - mu) / sigma) * pnorm(nu * (X - mu) / sigma)
  
  ggplot() +
    geom_line(mapping = aes(x = X, y = Y.axis), col = "red", lwd = 0.5) +
    labs(x = TeX(r"($X$)"),
         y = TeX(r"($f(x)$ )")) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))
}

SN_vis(0, 0.3, 1)

