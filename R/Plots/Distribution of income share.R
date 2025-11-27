data(SumData)

dec_cols <- paste0("d", 1:10)

dec_long <- SumData %>%
  select(Country, all_of(dec_cols)) %>%
  pivot_longer(cols = starts_with("d"),
               names_to = "Decile",
               values_to = "IncomeShare")

ggplot(dec_long, aes(x = Decile, y = IncomeShare, group = Country)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  theme_minimal(base_size = 14) +
  labs(title = "Income Distribution: Decile Shares",
       x = "Decile",
       y = "Income Share (%)")
