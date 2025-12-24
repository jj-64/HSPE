# library(ggplot2)
# library(dplyr)
# library(tidyr)

# Limited_data <- read_excel("data/Limited data.xlsx", sheet = "Parameters")
# Param_limited_data = as.data.frame(Limited_data)
# save(Param_limited_data, file = "data/Param_limited_data.rda")

data("HC_limited_data")
HC_long <- HC_limited_data %>%
  tidyr::pivot_longer(cols = c(observed_HC, HC),
               names_to = "Model",
               values_to = "Headcount")


## if only for certain poverty line threshold
HC_long = subset(HC_limited_data, threshold == 0.5)

ggplot(HC_long, aes(x = Country, y = HC*100, color = model, group = model)) +
  geom_point(size = 1.5, alpha=0.5) +
  geom_point(aes(x=Country, y = observed_HC*100), color = "gray", size=1.5, alpha = 0.9) +
  geom_line(alpha=0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "Poverty Headcount Comparison: Observed vs. Fitted Models",
       x = "Country",
       y = "Poverty Headcount (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("lightblue", "lightgreen", "plum")) +
  theme_classic()
