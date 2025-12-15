######### general theme
theme_set(theme_void(base_family = "Roboto"))

theme_update(
  axis.text.x = element_text(color = "black", face = "bold", size = 12, margin = margin(t = 6)),
  axis.text.y = element_text(color = "black", size = 12, hjust = 1, margin = margin(r = 6), family = "Roboto Mono"),
  axis.line.x = element_line(color = "black", size = 1),panel.grid.major.y = element_line(color = "grey90", size = .6),
  plot.background = element_rect(fill = "white", color = "white"),
  plot.margin = margin(rep(20, 4)))

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")
##"#7F3C8D" "#11A579" "#3969AC" "#F2B701" "#E73F74" "#80BA5A" "#E68310" "#A5AA99"
## 1:purple, 2:green,3:blue,4:yellow, 5:red, 6:limegreen, 7:orange, 8:gray

############
load("Diag_grouped_data")

## AIC
aggregate(AIC ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= AIC, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1e6,1e6)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")

## MSE
aggregate(MSE ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= MSE, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")

## RMSE
aggregate(RMSE ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= RMSE, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")

## MSE
aggregate(MSE ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= MSE, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")

## BIC
aggregate(BIC ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= BIC, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1e5,1e6)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")


## KS
aggregate(KS ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= KS, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")

## LogLik
aggregate(logLik ~ model, data = Diag_grouped_data, FUN = mean)
ggplot(Diag_grouped_data, aes(x=model, y= logLik, color = model, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1e5, 1e5)) +
  scale_color_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none")
