library(ggplot2)
library(latex2exp)
library(tidyverse)
library(patchwork)

rm(list = ls())

data <- load("Simulation_prec.Set.2.Rdata")

nu.d1.vec <- c(-25, -10, -5, -2.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2.5, 5, 10, 25)

mat.all <- matrix(as.vector(store.psi_hat), ncol = 9, byrow = T)
mat.all <- cbind(mat.all, nu.d1.vec)

df.3 <- as.data.frame(mat.all[,c(1,4,7,10)])
df.4 <- as.data.frame(mat.all[,c(2,5,8,10)])
df.5 <- as.data.frame(mat.all[,c(3,6,9,10)])

colnames(df.5)  <- colnames(df.4)  <- colnames(df.3) <- c("0.35", "0.5", "0.65", "nu.d1")

p1 <- df.3 |> 
  pivot_longer(cols = !nu.d1, names_to = "gamma", values_to = "psi.hat") |>
  ggplot() +
  geom_ribbon(aes(x = nu.d1, ymin = 0.285, ymax = 0.315), fill = "gray80") +
  geom_line(aes(x = nu.d1, y = psi.hat, lty = gamma)) +
  labs(x = TeX(r"($\nu_{d_1}$)"),
       y = TeX(r"($\hat{\psi}$)"),
       lty = TeX(r"($\gamma_{11} = \gamma_{12}$       )")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "top",
        plot.background = element_rect(fill = "gray100")) +
  scale_y_continuous(breaks = c(0.27, 0.3, 0.33), limits = c(0.27, 0.33))

p2 <- df.4 |> 
  pivot_longer(cols = !nu.d1, names_to = "gamma", values_to = "psi.hat") |>
  ggplot() +
  geom_ribbon(aes(x = nu.d1, ymin = 0.38, ymax = 0.42), fill = "gray80") +
  geom_line(aes(x = nu.d1, y = psi.hat, lty = gamma)) +
  labs(x = TeX(r"($\nu_{d_1}$)"),
       y = TeX(r"($\hat{\psi}$)"),
       lty = TeX(r"($\gamma_{11} = \gamma_{12}$       )")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "none",
        plot.background = element_rect(fill = "gray100")) +
  scale_y_continuous(breaks = c(0.36, 0.4, 0.44), limits = c(0.36, 0.44))

p3 <- df.5 |> 
  pivot_longer(cols = !nu.d1, names_to = "gamma", values_to = "psi.hat") |>
  ggplot() +
  geom_ribbon(aes(x = nu.d1, ymin = 0.475, ymax = 0.525), fill = "gray80") +
  geom_line(aes(x = nu.d1, y = psi.hat, lty = gamma)) +
  labs(x = TeX(r"($\nu_{d_1}$)"),
       y = TeX(r"($\hat{\psi}$)"),
       lty = TeX(r"($\gamma_{11} = \gamma_{12}$       )")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "none",
        plot.background = element_rect(fill = "gray100")) +
  scale_y_continuous(breaks = c(0.45, 0.5, 0.55), limits = c(0.45, 0.55))

p1 + p2 + p3 +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "top"))

