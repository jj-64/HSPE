ggplot(SumData, aes(x = Mean, y = Gini)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "darkblue") +
  theme_minimal(base_size = 14) +
  labs(title = "Gini Index vs Mean Income",
       x = "Mean Disposable Income",
       y = "Gini Index")
