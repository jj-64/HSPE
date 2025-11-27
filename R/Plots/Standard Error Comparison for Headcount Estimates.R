HC_se_long <- HC_limited_data %>%
  pivot_longer(
    cols = c(FISK_H_SE, LN_H_SE, NP_H_SE),
    names_to = "Model",
    values_to = "SE"
  )

ggplot(HC_se_long, aes(x = Country, y = SE, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 14) +
  labs(title = "Standard Errors of Model-based Headcount Ratios",
       x = "Country",
       y = "Standard Error") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
