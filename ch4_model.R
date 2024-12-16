
# install.packages("GGally")
# install.packages("tidyverse")
# install.packages("readr")
# install.packages("corrplot")
# install.packages("mice")
# install.packages("psych")

# Load required libraries
library(GGally)
library(tidyverse)
library(readr)
library(corrplot)

# Readme: Outputs are saved in the Output Folder and follow the naming system; [checking, analysis]_[plot no]_[plot or table type].
# Example: check_p1_scatter.png for data checking, plot 1, and scatter plot

# Set the working directory
setwd("/Users/phoopyezone/Documents/IRRI_MSI/data")

# Create output folder if it doesn't exist
if (!dir.exists("outputs")) {
  dir.create("outputs")
}

# Load the data
data <- read_csv("first_run.csv") %>%
  filter(ch4 >= 0)

# Check the data structure
str(data)

############ Data Exploration and Checking ################


# Count the levels of categorical variables
cat_vars <- data %>% select_if(is.character) %>% summarise_all(~n_distinct(.))
print(cat_vars)

# No of obs by country
sample_size <- data %>% group_by(country) %>% summarise(count = n())
print(sample_size)

# Missing values check
missing <- data %>% select_if(is.numeric) %>% summarise_all(~sum(is.na(.)))
print(missing)

# Scatter plot matrix
check_p1_scatter <- ggpairs(sub_data,
                            mapping = aes(color = treatment),
                            title = "Scatter Plot Matrix of Numerical Variables Colored by Treatment",
                            lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
                            diag = list(continuous = wrap("densityDiag")),
                            upper = list(continuous = wrap("cor", size = 3))
)
print(check_p1_scatter)

# CH4 vs SOC colored by treatment
check_p2_scatter <- ggplot(sub_data, aes(x = SOC, y = ch4, color = treatment)) +
  geom_point(size = 2, alpha = 0.6) +
  theme_classic() +
  labs(title = "Scatter Plot of CH4 vs SOC Colored by Treatment",
       x = "Soil Organic Carbon (SOC)",
       y = "CH4 Emission (ch4)")
print(check_p2_scatter)

# CH4 vs SOC colored by country
check_p3_scatter <- ggplot(sub_data, aes(x = SOC, y = ch4, color = country)) +
  geom_point(size = 2, alpha = 0.6) +
  theme_classic() +
  labs(title = "Scatter Plot of CH4 vs SOC Colored by Country",
       x = "Soil Organic Carbon (SOC)",
       y = "CH4 Emission (ch4)")
print(check_p3_scatter)

# CH4 boxplot, colored by treatment
check_p4_boxplot <- ggplot(sub_data, aes(x = treatment, y = ch4, color = treatment)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Box Plot of CH4 Emissions by Treatment",
       x = "Treatment",
       y = "CH4 Emission (ch4)")
print(check_p4_boxplot)

# CH4 box plot, colored by treatment with data points colored by country
check_p5_boxplot <- ggplot(sub_data, aes(x = treatment, y = ch4, color = treatment)) +
  geom_boxplot() +
  geom_jitter(aes(color = country), width = 0.2, size = 1.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "Box Plot of CH4 Emissions by Treatment with Data Points Colored by Country",
       x = "Treatment",
       y = "CH4 Emission (ch4)")
print(check_p5_boxplot)

# CH4 histogram
check_p6_histogram <- ggplot(sub_data, aes(x = ch4)) +
  geom_histogram(bins = 30) +
  theme_minimal() +
  ggtitle("Distribution of CH4")
print(check_p6_histogram)

#N rate histogram
check_p7_histogram <- ggplot(sub_data, aes(x = N_rate)) +
  geom_histogram(bins = 30) +
  theme_minimal() +
  ggtitle("Distribution of N_rate")
print(check_p7_histogram)



# Boxplot of CH4 by country, colored by country
check_p8_boxplot <- ggplot(sub_data, aes(x = country, y = ch4, color = country)) +
  geom_boxplot() +
  geom_jitter(aes(color = treatment), width = 0.2, size = 1, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Box Plot of CH4 Emissions by Country with Points Colored by Treatment",
       x = "Country",
       y = "CH4 Emission (ch4)")
print(check_p8_boxplot)




# Correlation matrix for numeric variables
cor_matrix <- cor(sub_data[c("N_rate", "SOC", "clay", "ch4")], use = "complete.obs")
check_p9_corrplot <- corrplot(cor_matrix, method = "color")

# Save all plots
# Scatter plot matrix
ggsave(filename = "outputs/check_p1_scatter.png", plot = check_p1_scatter, width = 10, height = 10)

# CH4 vs SOC colored by treatment
ggsave(filename = "outputs/check_p2_scatter.png", plot = check_p2_scatter, width = 10, height = 7)

# CH4 vs SOC colored by country
ggsave(filename = "outputs/check_p3_scatter.png", plot = check_p3_scatter, width = 10, height = 7)

# Box plot for CH4 by treatment
ggsave(filename = "outputs/check_p4_boxplot.png", plot = check_p4_boxplot, width = 10, height = 7)

# Box plot for CH4 by treatment with data points colored by country
ggsave(filename = "outputs/check_p5_boxplot.png", plot = check_p5_boxplot, width = 15, height = 10)

# Boxplot of CH4 by country, colored by country with axis labels at 45 degrees
ggsave(filename = "outputs/check_p8_boxplot.png", plot = check_p8_boxplot, width = 10, height = 7)

# Histogram for CH4
ggsave(filename = "outputs/check_p6_histogram.png", plot = check_p6_histogram, width = 10, height = 7)

# Histogram for N_rate
ggsave(filename = "outputs/check_p7_histogram.png", plot = check_p7_histogram, width = 10, height = 7)



########## Assumptions check ############

# Assumption 1: Linearity
# Plot residuals vs. predicted values
model <- lm(ch4 ~ SOC + N_rate, data = sub_data)
residuals_vs_fitted <- ggplot(model, aes(.fitted, .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals")
print(residuals_vs_fitted)

# Assumption 2: Independence of Observations
# Durbin-Watson test for autocorrelation
# install.packages("lmtest")
library(lmtest)
durbin_watson_result <- dwtest(model)
print(durbin_watson_result)

# Assumption 3: Homoscedasticity
# Breusch-Pagan test
# install.packages("lmtest")
library(lmtest)
bptest_result <- bptest(model)
print(bptest_result)

# Assumption 4: Normality of Residuals
# Q-Q plot of residuals
qq_plot <- ggplot(model, aes(sample = .resid)) +
  stat_qq() +
  stat_qq_line() +
  theme_classic() +
  labs(title = "Q-Q Plot of Residuals")
print(qq_plot)

# Additional check: Multicollinearity
# Calculate VIF (Variance Inflation Factor)
# install.packages("car")
library(car)
vif_values <- vif(model)
print(vif_values)

# Save assumption plots
# Residuals vs Fitted Values
ggsave(filename = "outputs/assumption_residuals_vs_fitted.png", plot = residuals_vs_fitted, width = 10, height = 7)

# Q-Q Plot of Residuals
ggsave(filename = "outputs/assumption_qq_plot.png", plot = qq_plot, width = 10, height = 7)



########## Analysis ##########

# Summary statistics for numerical variables
# install.packages("psych")
library(psych)
describe_stats <- describe(sub_data)[, c("n", "mean", "sd", "median", "trimmed", "min", "max", "range", "skew", "se")]
print(describe_stats)