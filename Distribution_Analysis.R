############################################
# Setup and Data Input
############################################

# Load required packages
library(ggplot2)
library(fitdistrplus)
library(dplyr)
library(tidyr)
library(scales)
library(VGAM)

# Simulated example: heavy head and long tail
# set.seed(123)
# data_vector <- c(rep(1, 250), rep(2, 50), rep(3, 30), sample(4:500, 71, replace = TRUE))


####

############################################
# 1. Fit Candidate Distributions using fitdistrplus and poweRlaw
############################################

# (A) Lognormal
fit_lognorm <- fitdist(data_vector, "lnorm")

# (B) Exponential
fit_exp <- fitdist(data_vector, "exp")

# (C) Half-Normal: Define a custom density function and wrapper
dhalfnorm <- function(x, sigma, log = FALSE) {
  out <- rep(0, length(x))
  valid <- x >= min(data_vector)  
  if (!log) {
    out[valid] <- (sqrt(2) / sqrt(sigma^2 * pi)) * exp(- (x[valid] - min(data_vector))^2 / (2 * sigma^2))
  } else {
    out[valid] <- log(sqrt(2) / sqrt(sigma^2 * pi)) - (x[valid] - min(data_vector))^2 / (2 * sigma^2)
  }
  return(out)
}
dhalfnorm_wrapper <- function(x, sigma, log = FALSE) {
  dhalfnorm(x, sigma, log)
}
start_sigma <- sd(data_vector) / 2
fit_halfnorm <- fitdist(data_vector, "halfnorm_wrapper", start = list(sigma = start_sigma))

# (D) Power Law with poweRlaw

# Create power law object and set x_min = min(data_vector)
m_pl_cont <- conpl$new(data_vector)
m_pl_cont$setXmin(min(data_vector))
est <- estimate_pars(m_pl_cont)
m_pl_cont$setPars(est$pars)
alpha_cont <- est$pars[1]
cat("Estimated alpha for power law:", alpha_cont, "\n")
ll_cont <- dist_ll(m_pl_cont)
aic_pl_cont <- -2 * ll_cont + 2 * length(est$pars)

############################################
# 2. AIC Summary Table
############################################
aic_table <- data.frame(
  Model = c("Lognormal", "Exponential", "Half Normal", "Power Law"),
  AIC = c(AIC(fit_lognorm), AIC(fit_exp), AIC(fit_halfnorm), aic_pl_cont)
)
print(aic_table)

############################################
# 3. Plotting: CDF and CCDF for All Models
############################################

# Define colors and later relabel them in the legend
model_colors <- c("Empirical" = "black", 
                  "Lognormal" = "firebrick", 
                  "Exponential" = "seagreen", 
                  "HalfNormal" = "steelblue", 
                  "PowerLaw" = "darkorchid")

# (A) Prepare Empirical CDF Data
sorted_data <- sort(data_vector)
empirical_cdf <- ecdf(data_vector)(sorted_data)
df_cdf <- data.frame(Size = sorted_data, Empirical = empirical_cdf)

# (B) Define Fitted CDF Functions
cdf_lognorm <- function(x) { plnorm(x, meanlog = fit_lognorm$estimate["meanlog"],
                                    sdlog = fit_lognorm$estimate["sdlog"]) }
cdf_exp <- function(x) { pexp(x, rate = fit_exp$estimate["rate"]) }
cdf_halfnorm <- function(x) { 2 * pnorm((x - min(data_vector)) / fit_halfnorm$estimate["sigma"]) - 1 }
# For power law, the CDF is:
#   F(x) = 1 - (x/xmin)^(1-alpha)
cdf_pl_cont <- function(x) { 
  ifelse(x < min(data_vector), 0, 1 - (x / min(data_vector))^(1 - alpha_cont))
}

# (C) Create a CDF Plot.
p_cdf <- ggplot() +
  geom_line(data = df_cdf, aes(x = Size, y = Empirical, color = "Empirical"), size = 1) +
  stat_function(fun = cdf_lognorm, aes(color = "Lognormal"), size = 1) +
  stat_function(fun = cdf_exp, aes(color = "Exponential"), size = 1) +
  stat_function(fun = cdf_halfnorm, aes(color = "HalfNormal"), size = 1) +
  stat_function(fun = cdf_pl_cont, aes(color = "PowerLaw"), size = 1) +
  scale_color_manual(
    values = model_colors,
    breaks = c("Empirical", "Lognormal", "Exponential", "HalfNormal", "PowerLaw"),
    labels = c("Empirical", "Lognormal", "Exponential", "Half Normal", "Power Law"),
    drop = FALSE
  ) +
  labs(title = "CDF: Empirical vs. Fitted Distributions",
       x = "Size", y = "Cumulative Probability", color = "Model") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_cdf)

# (D) CCDF Plot: Plot CCDF against log size.
empirical_ccdf <- 1 - (seq_along(sorted_data) - 1) / length(sorted_data)
df_ccdf <- data.frame(Size = sorted_data, CCDF = empirical_ccdf)

ccdf_lognorm <- function(x) { 1 - cdf_lognorm(x) }
ccdf_exp <- function(x) { 1 - cdf_exp(x) }
ccdf_halfnorm <- function(x) { 1 - cdf_halfnorm(x) }
ccdf_pl_cont <- function(x) { 1 - cdf_pl_cont(x) }

x_vals <- seq(min(data_vector), max(data_vector), length.out = 1000)
df_ccdf_fit <- data.frame(
  x = x_vals,
  Lognormal = ccdf_lognorm(x_vals),
  Exponential = ccdf_exp(x_vals),
  HalfNormal = ccdf_halfnorm(x_vals),
  PowerLaw = ccdf_pl_cont(x_vals)
)
df_ccdf_fit_long <- pivot_longer(df_ccdf_fit, cols = -x, names_to = "Model", values_to = "CCDF")
df_ccdf_fit_long$Model <- factor(df_ccdf_fit_long$Model, 
                                 levels = c("Empirical", "Lognormal", "Exponential", "HalfNormal", "PowerLaw"))

# Ensure the data are sorted by x within each model
df_ccdf_fit_long <- df_ccdf_fit_long %>% arrange(Model, x)

p_ccdf <- ggplot() +
  geom_point(data = df_ccdf, aes(x = Size, y = CCDF, color = "Empirical"), size = 1.5, alpha = 0.7) +
  geom_line(data = df_ccdf_fit_long, aes(x = x, y = CCDF, color = Model), size = 1) +
  scale_color_manual(
    values = model_colors,
    limits = c("Empirical", "Lognormal", "Exponential", "HalfNormal", "PowerLaw"),
    labels = c("Empirical", "Lognormal", "Exponential", "Half Normal", "Power Law"),
    drop = FALSE
  ) +
  scale_x_continuous(trans = "log10", labels = scales::comma_format()) +
  scale_y_log10(limits = c(1e-4, NA), labels = scales::comma_format()) +
  labs(title = "CCDF: Empirical vs. Fitted Distributions (Log Scale)",
       x = "Size (log scale)", y = "CCDF (log scale)", color = "Model") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_ccdf)


