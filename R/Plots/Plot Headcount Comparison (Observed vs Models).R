library(ggplot2)
library(dplyr)
library(tidyr)

data(HC_limited_data)

HC_long <- HC_limited_data %>%
  pivot_longer(cols = c(Observed, FISK_H, LN_H, NP_H),
               names_to = "Model",
               values_to = "Headcount")

## if only for certain poverty line threshold
HC_long2 = subset(HC_long, threshold == 0.2)

ggplot(HC_long2, aes(x = Country, y = Headcount*100, color = Model, group = Model)) +
  geom_point(size = 1.5, alpha=0.5) +
  geom_line(alpha=0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "Poverty Headcount Comparison: Observed vs. Fitted Models",
       x = "Country",
       y = "Poverty Headcount (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
