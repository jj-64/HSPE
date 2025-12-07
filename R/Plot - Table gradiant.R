## readxl || data wrangling + ggplot2 || adjust colors || Carto palettes ||sina plots || halfeye plots
## ridgeline plots || beeswarm plots || off-set jitter || custom fonts ||reshape

# List of required packages
required_packages <- c("reshape2","ggplot2","readxl","dplyr")

# Install missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages)
}

# Load all required libraries
lapply(required_packages, library, character.only = TRUE)

############## Load Data ########################
# HC_micro_data <- as.data.frame(read_excel("data/Micro data.xlsx"))
# save(HC_micro_data, file = "data/HC_micro_data.rda")

load("data/HC_micro_data.rda")

data_H = HC_micro_data
data_H = na.omit(data_H)

data_H$abs_delta = abs(data_H$observed_HC - data_H$HC)
data_H$rel_delta = abs(data_H$observed_HC - data_H$HC)/data_H$observed_HC

data_H$abs_group = cut(data_H$abs_delta, breaks= c(0, 0.005, 0.01,0.03,0.05, Inf),
                       labels = c("∆ < 0.5","0.5 ≤ ∆ < 1","1 ≤ ∆ < 3","3 ≤ ∆ < 5","∆ ≥ 5"),
                       include.lowest = TRUE,
                       right = FALSE)

data_H$rel_group = cut(data_H$rel_delta, breaks= c(0, 0.05, 0.1,0.3,0.5, Inf),
                       labels = c("ε < 5%","5% ≤ ε < 10%","10% ≤ ε < 30%","30% ≤ ε < 50%","ε ≥ 50%"),
                       include.lowest = TRUE,
                       right = FALSE)

## create a summarize table for absolute difference ------------------

abs_delta_Sum <- data_H %>%
  group_by(model, abs_group) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

abs_delta_Sum <- as.data.frame(abs_delta_Sum )

# Plot using ggplot2
ggplot(abs_delta_Sum, aes(x = model, y = abs_group, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(fill = "Intensity") +
  theme_minimal() +
  coord_fixed() +
  theme_classic()+
  #scale_y_reverse() + # Flips y-axis for standard matrix view
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())  +
  geom_text(aes(label = round(prop*100, 0)), color = "white", size = 3) + # Adds cell values
  #scale_x_continuous(breaks = 1:nrow(Data1), labels = x_labels) + # Custom x-axis labels
  #scale_y_continuous(breaks = 1:(ncol(Data1)-1), labels = y_labels) + # Custom y-axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(angle = 360, vjust = 0)) +
  ylab("") +
  xlab("")

## create a summarize table for relative difference ------------------

rel_delta_Sum <- data_H %>%
  group_by(model, rel_group) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

rel_delta_Sum <- as.data.frame(rel_delta_Sum )

# Plot using ggplot2
ggplot(rel_delta_Sum, aes(x = model, y = rel_group, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(fill = "Intensity") +
  theme_minimal() +
  coord_fixed() +
  theme_classic()+
  #scale_y_reverse() + # Flips y-axis for standard matrix view
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())  +
  geom_text(aes(label = round(prop*100, 0)), color = "white", size = 3) + # Adds cell values
  #scale_x_continuous(breaks = 1:nrow(Data1), labels = x_labels) + # Custom x-axis labels
  #scale_y_continuous(breaks = 1:(ncol(Data1)-1), labels = y_labels) + # Custom y-axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(angle = 360, vjust = 0)) +
  ylab("") +
  xlab("")

