## readxl || data wrangling + ggplot2 || adjust colors || Carto palettes ||sina plots || halfeye plots
## ridgeline plots || beeswarm plots || off-set jitter || custom fonts ||reshape

# List of required packages
required_packages <- c("tidyverse", "colorspace", "rcartocolor", "ggforce",
                       "ggdist","ggridges","ggbeeswarm","gghalves","systemfonts","tidyr")

# Install missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages)
}

# Load all required libraries
lapply(required_packages, library, character.only = TRUE)

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

pov_line = 0.8

## Limited Data -------------------
data("HC_limited_data")

data_limited <- HC_limited_data %>%
  na.omit() %>%
  filter(threshold == pov_line) %>%
  mutate(
    abs = observed_HC - HC,
    rel = abs / observed_HC
  )

## boxplot absolute
ggplot(data_limited, aes(x = model, y = abs, color = model, fill = model)) +
  scale_y_continuous() + #n.breaks = c(-10,0,10,20,30,40)
  scale_color_manual(values = my_pal[c(2,4,5)], guide = "none") +
  scale_fill_manual(values = my_pal[c(2,4,5)], guide = "none") +
  geom_violin(
    #aes(fill = model,fill = after_scale(colorspace::lighten(fill, .5))),
    #size = 1, bw=0.9
  )

## boxplot Relative
ggplot(data_limited, aes(x = model, y = rel*100, color = model, fill = model)) +
  scale_y_continuous(limits = c(-100,100)) + #n.breaks = c(-10,0,10,20,30,40)
  scale_color_manual(values = my_pal[c(2,4,5)], guide = "none") +
  scale_fill_manual(values = my_pal[c(2,4,5)], guide = "none") +
  geom_violin(
    #aes(fill = model,fill = after_scale(colorspace::lighten(fill, .5))),
    #size = 1, bw=0.9
  )

## Grouped Data -------------------
data("HC_grouped_data")

data_grouped <- HC_grouped_data %>%
  na.omit() %>%
  filter(threshold == pov_line) %>%
  mutate(
    abs = observed_HC - HC,
    rel = abs / observed_HC
  )

## boxplot absolute
ggplot(data_grouped, aes(x = model, y = abs, color = model, fill = model)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) + #n.breaks = c(-10,0,10,20,30,40)
  scale_color_manual(values = my_pal[c(8,1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8, 1:7)], guide = "none") +
  geom_violin(
    #aes(fill = model,fill = after_scale(colorspace::lighten(fill, .5))),
    #size = 1, bw=0.9
  )

## boxplot relative
ggplot(data_grouped, aes(x = model, y = rel*100, color = model, fill = model)) +
  scale_y_continuous(limits = c(-100,100)) + #n.breaks = c(-10,0,10,20,30,40)
  scale_color_manual(values = my_pal[c(8,1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8,1:7)], guide = "none") +
  geom_violin(
    #aes(fill = model,fill = after_scale(colorspace::lighten(fill, .5))),
    #size = 1, bw=0.9
  )

## Micro Data -------------------
data("HC_micro_data")

data_micro <- HC_micro_data %>%
  na.omit() %>%
  filter(threshold == pov_line, model != "B2") %>%
  mutate(
    abs = observed_HC - HC,
    rel = abs / observed_HC
  )

data_micro[data_micro$model == "BetaPrime", "model"] = "B2"

## boxplot absolute
ggplot(data_micro, aes(x = model, y = abs, color = model, fill = model)) +
  scale_y_continuous(limits = c(-0.5,0.5)) + #n.breaks = c(-10,0,10,20,30,40)
  scale_color_manual(values = my_pal[c(8,1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8,1:7)], guide = "none") +
  geom_violin(
    #aes(fill = model,fill = after_scale(colorspace::lighten(fill, .5))),
    #size = 1, bw=0.9
    )

## boxplot relative
ggplot(data_micro, aes(x = model, y = rel*100, color = model, fill = model)) +
  scale_y_continuous(limits = c(-100,100)) + #n.breaks = c(-10,0,10,20,30,40)
  scale_color_manual(values = my_pal[c(8,1:7)], guide = "none") +
  scale_fill_manual(values = my_pal[c(8,1:7)], guide = "none") +
  geom_violin(
    #aes(fill = model,fill = after_scale(colorspace::lighten(fill, .5))),
    #size = 1, bw=0.9
  )
# g+ ggdist::stat_halfeye(
#   aes(fill = Group, fill = after_scale(colorspace::lighten(fill, .7))),
#   adjust = .2, position = position_nudge(x = -.3))
