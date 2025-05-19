rm(list = ls())
library(ggplot2)

data <- read.csv("obj2_weight_difference.csv")
wgt <- na.omit(data$weight_difference_kg)

breaks <- seq(floor(min(wgt))-1, ceiling(max(wgt)), by = 2)
hist.data <- hist(wgt, breaks = breaks, right = FALSE, plot = FALSE)

df <- data.frame(vals = hist.data$mids, freq = hist.data$counts)

ggplot(df, aes(x = vals, y = freq)) +
  geom_line(color = "dodgerblue4", linewidth = 0.5) +
  geom_point(color = "white", size = 4, shape = 21, fill = "white") + 
  geom_point(color = "dodgerblue4", size = 2, shape = 1) + 
  labs(x = "Weight difference (kgs)",
    y = "Frequency") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 13, hjust = 0.5))
