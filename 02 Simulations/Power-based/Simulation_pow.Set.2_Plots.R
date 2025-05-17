library(ggplot2)
library(latex2exp)
library(tidyverse)
library(patchwork)

rm(list = ls())

data <- load("Simulation_pow.Set.2.Rdata")

nu.d1.vec <- c(-25, -10, -5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 5, 10, 25)

mat.all <- matrix(as.vector(store.pow), ncol = 9, byrow = T)
mat.all <- cbind(mat.all, nu.d1.vec)

df.3 <- as.data.frame(mat.all[,c(1,4,7,10)])
df.4 <- as.data.frame(mat.all[,c(2,5,8,10)])
df.5 <- as.data.frame(mat.all[,c(3,6,9,10)])

colnames(df.5)  <- colnames(df.4)  <- colnames(df.3) <- c("0.35", "0.5", "0.65", "nu.d1")

p1 <- df.3 |> 
  pivot_longer(cols = !nu.d1, names_to = "gamma", values_to = "pow.hat") |>
  ggplot() +
  geom_line(aes(x = nu.d1, y = pow.hat, lty = gamma)) +
  labs(x = TeX(r"($\nu_{d_1}$)"),
       y = TeX(r"($\hat{\psi}$)"),
       lty = TeX(r"($\gamma_{11} = \gamma_{12}$       )")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "top") +
  scale_y_continuous(limits = c(0.7, 0.9),
                     breaks = c(0.74, 0.8, 0.86))

p2 <- df.4 |> 
  pivot_longer(cols = !nu.d1, names_to = "gamma", values_to = "pow.hat") |>
  ggplot() +
  geom_line(aes(x = nu.d1, y = pow.hat, lty = gamma)) +
  labs(x = TeX(r"($\nu_{d_1}$)"),
       y = TeX(r"($\hat{\psi}$)"),
       lty = TeX(r"($\gamma_{11} = \gamma_{12}$       )")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "none") +
  scale_y_continuous(limits = c(0.7, 0.9),
                     breaks = c(0.74, 0.8, 0.86))

p3 <- df.5 |> 
  pivot_longer(cols = !nu.d1, names_to = "gamma", values_to = "pow.hat") |>
  ggplot() +
  geom_line(aes(x = nu.d1, y = pow.hat, lty = gamma)) +
  labs(x = TeX(r"($\nu_{d_1}$)"),
       y = TeX(r"($\hat{\psi}$)"),
       lty = TeX(r"($\gamma_{11} = \gamma_{12}$       )")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "none") +
  scale_y_continuous(limits = c(0.7, 0.9),
                     breaks = c(0.74, 0.8, 0.86))

p1 + p2 + p3 +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "top"))
