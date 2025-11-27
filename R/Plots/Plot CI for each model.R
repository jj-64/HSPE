data(CI_limited_data)

CI_long <- CI_limited_data %>%
  pivot_longer(
    cols = ends_with("_lower") | ends_with("_upper"),
    names_to = c("Model", ".value"),
    names_pattern = "(.+)_(lower|upper)"
  )

ggplot(CI_long, aes(x = threshold, color = Model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model),
              alpha = 0.25, color = NA) +
  geom_line(aes(y = (lower + upper) / 2), size = 1.2) +
  facet_wrap(~ Country) +
  theme_minimal(base_size = 14) +
  labs(title = "95% Confidence Intervals for Poverty Headcounts",
       x = "Poverty Threshold (fraction of mean income)",
       y = "Headcount") +
  guides(fill = guide_legend(title = "Model"),
         color = guide_legend(title = "Model"))
