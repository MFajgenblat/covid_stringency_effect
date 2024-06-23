#-------------------------------------------------------------------------------
# Packages
#-------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(reshape2)
library(zoo)
library(MASS)
library(abind)
library(cmdstanr)
library(splines)
library(posterior)
library(tidybayes)
library(ggtext)
library(patchwork)

#-------------------------------------------------------------------------------
# Data selection and modelling choices
#-------------------------------------------------------------------------------

# Country codes to be used
selected_country_codes <- c("ALB", "AUT", "BEL", "BGR", "BIH", "BLR", "CHE", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", "GRC", "HRV", "HUN", "IRL", "ISL", "ITA", "LTU", "LUX", "LVA", "MDA", "MLT", "NLD", "NOR", "POL", "PRT", "ROU", "RUS", "SRB", "SVK", "SVN", "SWE", "UKR")

# Starting date
startdate <- as.Date("2020-01-01")

# Ending date
enddate <- as.Date("2023-02-01")

# Maximal lag for the explanatory variables (days prior to the confirmed case date)
maxlag <- 45

# Maximal lag for the explanatory variables (days after the confirmed case date)
maxlead <- 10

# Number of basis functions along each dimension for the response surface
n_bf <- 4

#-------------------------------------------------------------------------------
# Data retrieval and preparation
#-------------------------------------------------------------------------------

data <- read.csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_nat_latest.csv") %>%
  # Select the relevant columns
  dplyr::select(CountryName, CountryCode, Date, ConfirmedCases, StringencyIndex_Average) %>%
  # Filter the selected countries
  filter(CountryCode %in% selected_country_codes) %>%
  # Convert the dates to the date format
  mutate(Date = as.Date(as.character(Date), "%Y%m%d")) %>%
  group_by(CountryCode) %>%
  # Decumulate the number of confirmed cases (i.e. compute the daily number)
  mutate(ConfirmedCases = ConfirmedCases - lag(ConfirmedCases, 1),
         # Compute the centered 7-day rolling mean
         ConfirmedCases = rollmeanr(ConfirmedCases, 7, fill = NA, align = "center")) %>%
  # Join with vaccination data
  left_join(mutate(read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv"),
                   date = as.Date(date)),
            by = c("CountryCode" = "iso_code", "Date" = "date")) %>%
  # Impute the fraction of people fully vaccinated for Luxemburg and Switzerland based the a 28-day lagged one-dose data
  mutate(FullyVaccinated = case_when(CountryCode %in% c("LUX", "CHE") ~ lag(people_vaccinated_per_hundred, 28),
                                     TRUE ~ people_fully_vaccinated_per_hundred)) %>%
  # Select the relevant columns
  dplyr::select(CountryName, CountryCode, Date, ConfirmedCases, StringencyIndex_Average, FullyVaccinated) %>%
  group_by(CountryCode) %>%
  # Compute the centered daily changes in confirmed cases
  mutate(Y = ConfirmedCases/lag(ConfirmedCases, n = 1) - 1,
         FullyVaccinated = case_when(FullyVaccinated == 0 ~ as.numeric(NA),
                                     TRUE ~ FullyVaccinated),
         # Set the fraction of people fully vaccinated to 0 at the first day
         FullyVaccinated = replace(FullyVaccinated, row_number() < min((FullyVaccinated > 0) * row_number(), na.rm = T), 0),
         # Fill gaps of missing vaccination data through linear interpolation, up to 50 consecutive missing data points
         FullyVaccinated = na.approx(FullyVaccinated, maxgap = 50, rule = 2)) %>%
  ungroup() %>%
  # Rename columns 
  rename(X1 = StringencyIndex_Average,
         X2 = FullyVaccinated) %>%
  # Select the relevant columns
  dplyr::select(CountryName, CountryCode, Date, Y, X1, X2) %>%
  # Add the lagged predictor columns 
  cbind(., do.call(cbind, lapply(maxlead:1, function(i) dplyr::select(mutate(., across(matches("X[1-2]$"), ~lead(., i), .names = paste0("{.col}_lead", i))), contains("lead"))))) %>%
  # Add the lagged predictor columns 
  cbind(., do.call(cbind, lapply(0:maxlag, function(i) dplyr::select(mutate(., across(matches("X[1-2]$"), ~lag(., i), .names = paste0("{.col}_lag", i))), contains("lag"))))) %>%
  # Remove the original predictor columns
  dplyr::select(-X1, -X2) %>%
  # Date filter
  filter(Date >= startdate,
         Date <= enddate,
         # Filter extreme changes
         !is.infinite(Y)) %>%
  # Omit rows with missing values
  drop_na()

#-------------------------------------------------------------------------------
# Constructing data arrays with the predictor x basis function x lag information
#-------------------------------------------------------------------------------

# Lagged basis functions for the observed data
X1 <- dplyr::select(data, starts_with("X1"))/100 - 0.5
X2 <- dplyr::select(data, starts_with("X2"))/100 - 0.5
bf_1 <- abind(lapply(1:ncol(X1), function(i) bs(X1[,i], knots = seq(-0.5, 0.5, length.out = n_bf-2), degree = 3, intercept = F)[,-(n_bf+1)]), along = 3)
bf_2 <- abind(lapply(1:ncol(X2), function(i) bs(X2[,i], knots = seq(-0.5, 0.5, length.out = n_bf-2), degree = 3, intercept = F)[,-(n_bf+1)]), along = 3)
bf_2D <- abind(lapply(1:ncol(X1), function(i) t(sapply(1:nrow(data), function(j) kronecker(bf_1[j,,i], bf_2[j,,i])))), along = 3)

# Grid-projected basis functions for prediction purposes (high resolution) 
bf_pred <- bs(seq(-0.5, 0.5, length.out = 101), knots = seq(-0.5, 0.5, length.out = n_bf-2), degree = 3, intercept = F)[,-(n_bf+1)]
grid_pred <- expand.grid(x = 1:nrow(bf_pred), y = 1:nrow(bf_pred))
bf_pred_2D <- t(sapply(1:nrow(grid_pred), function(i) kronecker(bf_pred[grid_pred[i,1],], bf_pred[grid_pred[i,2],])))

# Grid-projected basis functions for prediction purposes (low resolution) 
bf_pred_lowres <- bs(seq(-0.5, 0.5, length.out = 31), knots = seq(-0.5, 0.5, length.out = n_bf-2), degree = 3, intercept = F)[,-(n_bf+1)]
grid_pred_lowres <- expand.grid(x = 1:nrow(bf_pred_lowres), y = 1:nrow(bf_pred_lowres))
bf_pred_2D_lowres <- t(sapply(1:nrow(grid_pred_lowres), function(i) kronecker(bf_pred_lowres[grid_pred_lowres[i,1],], bf_pred_lowres[grid_pred_lowres[i,2],])))

#-------------------------------------------------------------------------------
# Constructing temporal basis functions
#-------------------------------------------------------------------------------

n_bf_days <- 50
N_days <- as.numeric(max(data$Date) - min(data$Date) + 1)
day_range <- as.numeric((min(data$Date):max(data$Date))) - as.numeric(min(data$Date))
bf_days <- bs(day_range, knots = seq(min(day_range), max(day_range), length.out = n_bf_days-2), degree = 3, intercept = F)[,-(n_bf_days+1)]
plot(bf_days %*% rnorm(ncol(bf_days)))

#-------------------------------------------------------------------------------
# Constructing the input data list
#-------------------------------------------------------------------------------

datalist <- list(N = nrow(data),
                 N_countries = length(unique(data$CountryName)),
                 N_cov = dim(aperm(abind(list(X1, X2), along = 3), c(3,1,2)))[1],
                 N_bf = ncol(bf_2D),
                 N_lags = dim(aperm(abind(list(X1, X2), along = 3), c(3,1,2)))[3],
                 country = as.numeric(factor(data$CountryName)),
                 X = aperm(abind(list(X1, X2), along = 3), c(3,1,2)),
                 bf = aperm(bf_2D, c(2,1,3)),
                 lag_range = seq(-1, 1, length.out = dim(aperm(abind(list(X1, X2), along = 3), c(3,1,2)))[3]),
                 y = data$Y,
                 N_pred = nrow(bf_pred_2D),
                 bf_pred = bf_pred_2D,
                 countrystart = pull(arrange(mutate(summarise(group_by(rowid_to_column(data), CountryName), first = first(rowid)), CountryName = factor(CountryName)), CountryName), first),
                 countryend = pull(arrange(mutate(summarise(group_by(rowid_to_column(data), CountryName), last = last(rowid)), CountryName = factor(CountryName)), CountryName), last),
                 N_pred_country = nrow(grid_pred_lowres),
                 X_pred_country = (grid_pred_lowres-1)/30 - 0.5,
                 bf_pred_lowres = bf_pred_2D_lowres,
                 # Normal likelihood (t = 0) or Student's t likelihood (t = 1)?
                 t = 1,
                 # Degrees of freedom of the Student's t likelihood
                 nu_t = 3,
                 # Should the model include a temporal effect (temporal = 1) or not (temporal = 0)?
                 temporal = 0,
                 # Prior specification
                 prior_intercept = c(0,1),
                 prior_tau = c(0,0.2),
                 prior_Omega = c(2),
                 prior_nu_mvt = c(2,0.1),
                 prior_alpha_smooth = c(0,1),
                 prior_rho_lags = c(5,5),
                 prior_alpha_lags = c(0, 5),
                 prior_sigma = c(0,1),
                 N_days = N_days,
                 N_bf_days = n_bf_days,
                 bf_days = bf_days,
                 day_range = seq(-1, 1, length.out = n_bf_days),
                 day = as.numeric(data$Date - min(data$Date) + 1),
                 prior_rho_temporal = c(5,5),
                 prior_alpha_temporal = c(0,1))

#-------------------------------------------------------------------------------
# Model fitting and saving results
#-------------------------------------------------------------------------------

# Compiling the model
model <- cmdstan_model("covid_stringency_effect.stan")
# Estimating the model
fit <- model$sample(data = datalist, iter_warmup = 500, iter_sampling = 500, refresh = 10, chains = 4, parallel_chains = 4)
# Converting the MCMC output
stanfit <- read_cmdstan_csv(fit$output_files(), format = "draws_matrix")
stanfit <- stanfit$post_warmup_draws
# Saving the MCMC output for later use
saveRDS(stanfit, "covid_stringency_effect_fit.rds")
# Reading existing MCMC output
stanfit <- readRDS("covid_stringency_effect_fit.rds")

#-------------------------------------------------------------------------------
# Computing additionally generated quantities and saving results
#-------------------------------------------------------------------------------

# Compiling the model
mod_gqs <- cmdstan_model("covid_stringency_effect_gqs.stan")
# Computing
fit_gqs <- mod_gqs$generate_quantities(stanfit, data = datalist, parallel_chains = 4)
# Converting the MCMC output
stanfit_gqs <- as_draws_df(fit_gqs$draws())
# Saving the MCMC output for later use
saveRDS(stanfit_gqs, "covid_stringency_effect_preds.rds")
# Reading existing MCMC output
stanfit_gqs <- readRDS("covid_stringency_effect_preds.rds")

#-------------------------------------------------------------------------------
# Metadata preparation
#-------------------------------------------------------------------------------

metadata_covariates <- data.frame(cov = 1:4,
                                  name = c("Intercept", "SI", "Vax", "Interaction"),
                                  Cov = factor(c("Intercept", "Stringency Index", "Fraction of population fully vaccinated", "Stringency Index x Fraction of population fully vaccinated"),
                                               levels = rev(c("Intercept", "Stringency Index", "Fraction of population fully vaccinated", "Stringency Index x Fraction of population fully vaccinated"))),
                                  Covbreak = factor(c("Intercept", "Stringency\nIndex", "Fraction of population\nfully vaccinated", "Stringency Index x Fraction\nof population fully vaccinated"),
                                                    levels = c("Intercept", "Stringency\nIndex", "Fraction of population\nfully vaccinated", "Stringency Index x Fraction\nof population fully vaccinated")),
                                  Covbreak_rev = factor(c("Intercept", "Stringency\nIndex", "Fraction of population\nfully vaccinated", "Stringency Index x Fraction\nof population fully vaccinated"),
                                                        levels = rev(c("Intercept", "Stringency\nIndex", "Fraction of population\nfully vaccinated", "Stringency Index x Fraction\nof population fully vaccinated"))))

metadata_countries <- data.frame(countrycode = data$CountryCode,
                                 Country = data$CountryName) %>%
  mutate(Country = case_when(Country == "Bosnia and Herzegovina" ~ "Bosn. and Herz.",
                             TRUE ~ Country)) %>%
  distinct() %>%
  arrange(countrycode) %>%
  mutate(Country = factor(Country),
         country = as.numeric(Country))
metadata_countries$Country_rev <- factor(metadata_countries$Country,
                                         levels = rev(levels(metadata_countries$Country)))

#-------------------------------------------------------------------------------
# Figure 1
#-------------------------------------------------------------------------------

mask <- kde2d(data$X1_lag14, data$X2_lag14, h = 20, lims = c(0, 100, 0, 100), n = c(1010,1010))

Fig1A <- stanfit_gqs %>%
  spread_draws(preds[id]) %>%
  right_join(data.frame(id = 1:(101^2),
                        expand.grid(SI = 0:100,
                                    Vax = 0:100))) %>%
  group_by(SI, Vax) %>%
  summarise(pred = mean(preds)) %>%
  mutate(pred = case_when(pred > 0.051 ~ 0.051,
                          pred < -0.051 ~ -0.051,
                          T ~ pred))  %>%
  ggplot() +
  geom_tile(aes(x = SI, y = Vax, fill = pred), color = NA) +
  geom_vline(xintercept = seq(20, 90, by = 20), color = "white", size = 0.4, alpha = 0.35) +
  geom_hline(yintercept = seq(20, 90, by = 20), color = "white", size = 0.4, alpha = 0.35) +
  geom_vline(xintercept = seq(10, 90, by = 20), color = "white", size = 0.2, alpha = 0.35) +
  geom_hline(yintercept = seq(10, 90, by = 20), color = "white", size = 0.2, alpha = 0.35) +
  geom_tile(data = subset(melt(mask$z), value < quantile(mask$z, 0.33)), aes(x = ((Var1-1)/10.09)*1.01-0.5, y = ((Var2-1)/10.09)*1.01-0.5), fill = "white", color = "white") +
  geom_point(data = data, aes(x = X1_lag14, y = X2_lag14), size = 0.35, alpha = 0.05, stroke = 0, shape = 16) +
  scale_fill_gradient2("Predicted change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                       limits = c(-0.051, 0.051), breaks = seq(-0.05, 0.05, by = 0.025), labels = c("<-0.05", 0.025, 0, 0.025, ">0.05"),
                       guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"), ticks.colour = "black")) +
  scale_x_continuous("Stringency Index", expand = c(0,0), breaks = c(-0.5, 20, 40, 60, 80, 100.5), labels = seq(0, 100, by = 20)) +
  scale_y_continuous("Fraction of population fully vaccinated (%)", expand = c(0,0), breaks = c(-0.5, 20, 40, 60, 80, 100.5), labels = seq(0, 100, by = 20)) +
  facet_wrap(~ "Posterior\nmean") +
  coord_equal() +
  theme(plot.margin = margin(10, 10, 10, 10),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.ticks = element_line(size = 0.2),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key = element_blank(),
        legend.position = "bottom")

Fig1B <- stanfit_gqs %>%
  spread_draws(preds[id]) %>%
  right_join(data.frame(id = 1:(101^2),
                        expand.grid(SI = 0:100,
                                    Vax = 0:100))) %>%
  group_by(SI, Vax) %>%
  summarise(pred_lower = quantile(preds, 0.025),
            pred_upper = (quantile(preds, 0.975))) %>%
  pivot_longer(c(pred_lower,pred_upper), names_to = "type", values_to = "pred") %>%
  mutate(pred = case_when(pred > 0.051 ~ 0.051,
                          pred < -0.051 ~ -0.051,
                          T ~ pred))  %>%
  mutate(type = case_when(type == "pred_lower" ~ "Lower 95% credible\ninterval bound",
                          type == "pred_upper" ~ "Upper 95% credible\ninterval bound")) %>%
  ggplot() +
  geom_tile(aes(x = SI, y = Vax, fill = pred), color = NA) +
  geom_vline(xintercept = seq(20, 90, by = 20), color = "white", size = 0.4/2, alpha = 0.35) +
  geom_hline(yintercept = seq(20, 90, by = 20), color = "white", size = 0.4/2, alpha = 0.35) +
  geom_vline(xintercept = seq(10, 90, by = 20), color = "white", size = 0.2/2, alpha = 0.35) +
  geom_hline(yintercept = seq(10, 90, by = 20), color = "white", size = 0.2/2, alpha = 0.35) +
  geom_tile(data = subset(melt(mask$z), value < quantile(mask$z, 0.33)), aes(x = ((Var1-1)/10.09)*1.01-0.5, y = ((Var2-1)/10.09)*1.01-0.5), fill = "white", color = "white") +
  geom_point(data = data, aes(x = X1_lag14, y = X2_lag14), size = 0.35/2, alpha = 0.05, stroke = 0, shape = 16) +
  scale_fill_gradient2("Predicted change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                       limits = c(-0.051, 0.051), breaks = seq(-0.05, 0.05, by = 0.025), labels = c("<-0.05", 0.025, 0, 0.025, ">0.05"),
                       guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"), ticks.colour = "black")) +
  scale_x_continuous("Stringency Index", expand = c(0,0), breaks = c(-0.5, 20, 40, 60, 80, 100.5), labels = seq(0, 100, by = 20)) +
  scale_y_continuous("Fraction of population\nfully vaccinated (%)", expand = c(0,0), breaks = c(-0.5, 20, 40, 60, 80, 100.5), labels = seq(0, 100, by = 20)) +
  coord_equal() +
  facet_wrap(~ type, ncol = 1) +
  theme(plot.margin = margin(10, 10, 10, 10),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key = element_blank(),
        legend.position = "bottom")

Fig1A + Fig1B + plot_layout(widths = c(1, 0.45), guides = 'collect') & theme(legend.position = "bottom")
ggsave("Figure_1.png", width = 16, height = 14, units = "cm", dpi = 600)
ggsave("Figure_1.pdf", width = 16, height = 14, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure 2
#-------------------------------------------------------------------------------

stanfit_gqs %>%
  spread_draws(preds_country[country,id]) %>%
  right_join(data.frame(id = 1:(31^2),
                        expand.grid(SI = seq(0, 100, length.out = 31),
                                    Vax = seq(0, 100, length.out = 31),
                                    country = 1:38))) -> preds

x = seq(-0.5, 0.5, length.out = 31)
beta_effects <- data.frame(country = NA, int = NA, SI = NA, Vax = NA)[-1,]
for (i in 1:38) {
  beta_effects <- rbind(beta_effects,
                        data.frame(country = i,
                                   .draw = filter(preds, country == i, SI == 0, Vax == 0)$.draw,
                                   int = filter(preds, country == i, SI == 0, Vax == 0)$preds_country,
                                   SI = apply(as.matrix(dplyr::select(pivot_wider(filter(preds, country == i, Vax == 0), values_from = preds_country, names_from = SI, id_cols = .draw), -.draw)), 1, function(y) lm(y ~ x)$coef[2])/100,
                                   Vax = apply(as.matrix(dplyr::select(pivot_wider(filter(preds, country == i, SI == seq(0, 100, length.out = 31)[1]), values_from = preds_country, names_from = Vax, id_cols = .draw), -.draw)), 1, function(y) lm(y ~ x)$coef[2])/100))
  print(i)
}
beta_effects %>%
  pivot_longer(!c(country, .draw)) %>%
  left_join(data.frame(name = c("int", "SI", "Vax"),
                       Name = factor(c("<b>Intercept (adjusted)</b><br><sub>Predicted change in<br>an unvaccinated population<br>without any stringency</sub>", "<b>Stringency Index</b><br><sub>Effect of an increase by<br>10 units (in an<br>unvaccinated population)</sub>", "<b>Fraction of population<br>fully vaccinated</b><br><sub>Effect of an increase by<br>10% (without any stringency)</sub>"),
                                     levels = c("<b>Intercept (adjusted)</b><br><sub>Predicted change in<br>an unvaccinated population<br>without any stringency</sub>", "<b>Stringency Index</b><br><sub>Effect of an increase by<br>10 units (in an<br>unvaccinated population)</sub>", "<b>Fraction of population<br>fully vaccinated</b><br><sub>Effect of an increase by<br>10% (without any stringency)</sub>")))) %>%
  left_join(metadata_countries) %>%
  mutate(value = case_when(name != "int" ~ value * 10,
                           T ~ value)) %>%
  ggplot() +
  geom_segment(data = data.frame(y = levels(metadata_countries$Country_rev)[seq(2, 38, by = 2)]), aes(x = -Inf, xend = +Inf, y = y, yend = y), linewidth = 4, color = "black", alpha = 0.05) +
  stat_interval(aes(x = value, y = Country_rev), .width = c(0.5, 0.8, 0.95, 0.99), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  scale_color_brewer("Credible interval", guide = guide_legend(override.aes = list(linewidth = 2))) +
  facet_wrap(~ Name, scales = "free_x") +
  scale_x_continuous("Predicted change in daily confirmed cases") +
  scale_y_discrete(limits = levels(metadata_countries$Country_rev)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_markdown(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))
ggsave("Figure_2.png", width = 16, height = 16, units = "cm", dpi = 600)
ggsave("Figure_2.pdf", width = 16, height = 16, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure 3
#-------------------------------------------------------------------------------

stanfit %>%
  spread_draws(delay_distribution[lag]) %>%
  ggplot() +
  stat_lineribbon(aes(x = -(lag - 11), y = delay_distribution), color = NA, .width = c(0.5, 0.8, 0.95, 0.99)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous("Lag (days prior to confirmed case)", expand = c(0,0), breaks = seq(-100, 100, by = 5)) +
  scale_y_continuous("Relative lag weight", expand = c(0,0), limits = c(0, NA)) +
  scale_fill_brewer("Credible\ninterval") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey95", size = 0.4),
        panel.grid.minor = element_line(colour = "grey95", size = 0.2),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key = element_blank(),
        legend.key.size = unit(0.5, "cm"))
ggsave("Figure_3.png", width = 16, height = 8, units = "cm", dpi = 600)
ggsave("Figure_3.pdf", width = 16, height = 8, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Obtaining the Numerical values presented in the manuscript
#-------------------------------------------------------------------------------

# Scale coefficient of the random intercept and slopes
stanfit %>%
  spread_draws(tau[cov]) %>%
  group_by(cov) %>%
  summarise(posterior_mean = mean(tau),
            CrI_lower = quantile(tau, 0.025),
            CrI_upper = quantile(tau, 0.975))

# Estimated degrees of freedom for the multivariate t distributed random effects
stanfit %>% spread_draws(nu_mvt) %>% pull(nu_mvt) -> nu_mvt
mean(nu_mvt)
hist(nu_mvt)
quantile(nu_mvt, c(0.025,0.975))

# Posterior probability of more inter-country variation for the effect of SI compared to Vax 
stanfit %>%
  spread_draws(tau[cov]) %>%
  pivot_wider(id_cols = c(.chain, .iteration, .draw), names_from = cov, values_from = tau) %>%
  mutate(difference = `3` - `2`) %>%
  summarise(posterior_prob = mean(difference < 0))

# Posterior mean SI effects
beta_effects %>%
  pivot_longer(!c(country, .draw)) %>%
  left_join(data.frame(name = c("int", "SI", "Vax"),
                       Name = factor(c("Intercept (adjusted)", "Stringency Index", "Fraction of population\nfully vaccinated (%)"),
                                     levels = c("Intercept (adjusted)", "Stringency Index", "Fraction of population\nfully vaccinated (%)")))) %>%
  left_join(metadata_countries) %>%
  filter(name == "SI") %>%
  group_by(Country) %>%
  summarise(value = mean(value)*10) %>%
  arrange(value) %>%
  View(.)

# Calculating posterior effect probabilities for Belarus
beta_effects %>%
  pivot_longer(!c(country, .draw)) %>%
  left_join(data.frame(name = c("int", "SI", "Vax"),
                       Name = factor(c("Intercept (adjusted)", "Stringency Index", "Fraction of population\nfully vaccinated (%)"),
                                     levels = c("Intercept (adjusted)", "Stringency Index", "Fraction of population\nfully vaccinated (%)")))) %>%
  left_join(metadata_countries) %>%
  filter(Country == "Belarus") %>%
  group_by(Name) %>%
  summarise(prob = mean(value < 0))

# Summarising the delay distribution
stanfit %>%
  spread_draws(delay_distribution[lag]) -> delay_distribution
delay_distribution %>%
  group_by(.draw) %>%
  filter(delay_distribution == max(delay_distribution)) %>%
  pull(lag) -> delay_distribution_mode
mean(delay_distribution_mode)
quantile(delay_distribution_mode, c(0.025, 0.975))
delay_distribution %>%
  group_by(.draw) %>%
  mutate(csum = cumsum(delay_distribution)) %>%
  filter(csum > 0.95) %>%
  filter(csum == min(csum)) %>%
  pull(lag) -> delay_distribution_95lag
hist(delay_distribution_95lag)
mean(delay_distribution_95lag)
quantile(delay_distribution_95lag, c(0.025,0.975))

positive_lags <- delay_distribution %>%
  filter(lag <= 10) %>%
  group_by(.draw) %>%
  summarize(delay_distribution = sum(delay_distribution)) %>%
  pull(delay_distribution)
median(positive_lags*100)
quantile(positive_lags*100, 0.025)
quantile(positive_lags*100, 0.975)

#-------------------------------------------------------------------------------
# Figure S1
#-------------------------------------------------------------------------------

data %>%
  mutate(CountryName = case_when(CountryName == "Bosnia and Herzegovina" ~ "Bosn. and Herz.",
                                 TRUE ~ CountryName)) %>%
  ggplot() +
  geom_point(aes(x = Date, y = CountryName), size = 0.1, color = "#0f94a8") +
  scale_x_date("Date", date_breaks = "2 months", date_labels = "%Y-%m", expand = c(0,0)) +
  scale_y_discrete(limits = levels(metadata_countries$Country_rev)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        axis.text.x = element_text(size = 8, angle = -45, hjust = 0),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_blank())
ggsave("Figure_S1.png", width = 16, height = 16, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S2
#-------------------------------------------------------------------------------

data %>%
  mutate(CountryName = case_when(CountryName == "Bosnia and Herzegovina" ~ "Bosn. and Herz.",
                                 TRUE ~ CountryName),
         Y_truncated = case_when(Y > 0.1 ~ 0.1,
                                 Y < -0.1 ~ -0.1,
                                 T ~ Y))  %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = Date, ymin = 0, ymax = X1_lag0/100), fill = "#2567a8", alpha = 0.95) +
  geom_ribbon(aes(x = Date, ymin = 0, ymax = -X2_lag0/100), fill = "#bf17ac", alpha = 0.95) +
  geom_point(aes(x = Date, y = Y, color = Y_truncated), size = 0.1) +
  scale_x_date("Date", date_breaks = "6 months", date_labels = "%Y-%m", expand = c(0,0)) +
  scale_y_continuous("Observed change in daily confirmed cases", expand = c(0,0)) +
  scale_color_gradient2("Observed change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                        limits = c(-0.101, 0.101), breaks = seq(-0.1, 0.1, by = 0.05), labels = c("<-0.1", 0.05, 0, 0.05, ">0.1"),
                        guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"), ticks.colour = "black")) +
  coord_cartesian(ylim = c(-1,1)) +
  facet_wrap(~ CountryName, ncol = 5) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(colour = "grey93", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 8, angle = -45, hjust = 0),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        axis.ticks = element_line(size = 0.3),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))
ggsave("Figure_S2.png", width = 16, height = 21, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S3
#-------------------------------------------------------------------------------

# Fetch the posterior samples for another run with a Gaussian likelihood
stanfit_normal <- readRDS("covid_stringency_effect_fit_Gaussian.rds")

rbind(cbind(model = "Student's t likelihood (nu = 3)", melt(stanfit[seq(1, nrow(stanfit), by = 10),startsWith(colnames(stanfit), "y_rep")])),
      cbind(model = "Gaussian likelihood", melt(stanfit_normal[seq(1, nrow(stanfit_normal), by = 10),startsWith(colnames(stanfit_normal), "y_rep")]))) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_line(data = data, aes(x = Y, color = "Observed data", size = "Observed data", alpha = "Observed data"), stat = "density", adjust = 2, bounds = c(-1,1)) +
  geom_line(aes(x = value, group = draw, color = "Posterior predictive simulation", size = "Posterior predictive simulation", alpha = "Posterior predictive simulation"), stat = "density", adjust = 2, bounds = c(-1,1)) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_x_continuous("Change in daily confirmed cases") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_manual("", values = c("#0f94a8", "black")) +
  scale_alpha_manual("", values = c(1, 0.1)) +
  scale_size_manual("", values = c(1, 0.1)) +
  facet_wrap(~ model) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        axis.ticks = element_line(size = 0.3),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))
ggsave("Figure_S3.png", width = 16, height = 8, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S4
#-------------------------------------------------------------------------------

posterior_preds <- stanfit_gq %>%
  spread_draws(preds[id]) %>%
  right_join(data.frame(id = 1:(101^2),
                        expand.grid(SI = 0:100,
                                    Vax = 0:100)))

FigS4A <- posterior_preds %>%
  group_by(SI, Vax) %>%
  summarise(pred = mean(preds)) %>%
  mutate(pred = case_when(pred > 0.051 ~ 0.051,
                          pred < -0.051 ~ -0.051,
                          T ~ pred))  %>%
  ggplot() +
  geom_tile(aes(x = SI, y = Vax, color = pred, fill = pred)) +
  geom_vline(xintercept = seq(20, 90, by = 20), color = "white", size = 0.4, alpha = 0.35) +
  geom_hline(yintercept = seq(20, 90, by = 20), color = "white", size = 0.4, alpha = 0.35) +
  geom_vline(xintercept = seq(10, 90, by = 20), color = "white", size = 0.2, alpha = 0.35) +
  geom_hline(yintercept = seq(10, 90, by = 20), color = "white", size = 0.2, alpha = 0.35) +
  geom_tile(data = subset(melt(mask$z), value < quantile(mask$z, 0.33)), aes(x = ((Var1-1)/10.09)*1.01-0.5, y = ((Var2-1)/10.09)*1.01-0.5), fill = "white", color = "white") +
  geom_point(data = data, aes(x = X1_lag14, y = X2_lag14), size = 0.35, alpha = 0.05, stroke = 0, shape = 16) +
  scale_color_gradient2("Predicted change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                        limits = c(-0.051, 0.051), breaks = seq(-0.05, 0.05, by = 0.025), labels = c("<-0.05", 0.025, 0, 0.025, ">0.05"),
                        guide = "none") +
  scale_fill_gradient2("Predicted change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                       limits = c(-0.051, 0.051), breaks = seq(-0.05, 0.05, by = 0.025), labels = c("<-0.05", 0.025, 0, 0.025, ">0.05"),
                       guide = "none") +
  scale_x_continuous("Stringency Index", expand = c(0,0), breaks = c(-0.5, 20, 40, 60, 80, 100.5), labels = seq(0, 100, by = 20)) +
  scale_y_continuous("Fraction of population\nfully vaccinated (%)", expand = c(0,0), breaks = c(-0.5, 20, 40, 60, 80, 100.5), labels = seq(0, 100, by = 20)) +
  coord_equal() +
  theme(plot.margin = margin(10, 10, 10, 10),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.ticks = element_line(size = 0.2),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key = element_blank(),
        legend.position = "bottom")
FigS4B <-rbind(cbind(Variable = "Stringency\nIndex", filter(posterior_preds, Vax == 0), value = filter(posterior_preds, Vax == 0)$SI),
               cbind(Variable = "Fraction of population\nfully vaccinated (%)", filter(posterior_preds, SI == 0), value = filter(posterior_preds, SI == 0)$Vax)) %>%
  mutate(Variable = factor(Variable, levels = c("Stringency\nIndex", "Fraction of population\nfully vaccinated (%)"))) %>%
  ggplot(aes(x = value, y = preds, group = .draw)) +
  geom_path(linewidth = 0.1, alpha = 0.02, color = "black") +
  scale_x_continuous("", expand = c(0,0)) +
  scale_y_continuous("Confirmed case change") +
  facet_wrap(~ Variable, ncol = 2, strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey93"),
        panel.spacing.x = unit(1, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        strip.placement = "outside",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

FigS4C <-rbind(cbind(Variable = "Stringency\nIndex", filter(posterior_preds, Vax == 0), value = filter(posterior_preds, Vax == 0)$SI),
               cbind(Variable = "Fraction of population\nfully vaccinated (%)", filter(posterior_preds, SI == 0), value = filter(posterior_preds, SI == 0)$Vax)) %>%
  mutate(Variable = factor(Variable, levels = c("Stringency\nIndex", "Fraction of population\nfully vaccinated (%)"))) %>%
  ggplot(aes(x = value, y = preds, group = .draw)) +
  geom_line(stat="smooth",method = "lm", linewidth = 0.1, alpha = 0.02, se = F, color = "black") +
  scale_x_continuous("", expand = c(0,0)) +
  scale_y_continuous("Confirmed case change") +
  facet_wrap(~ Variable, ncol = 2, strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey93"),
        panel.spacing.x = unit(1, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        strip.placement = "outside",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

x = seq(-0.5, 0.5, length.out = 101)
FigS4D <- rbind(cbind(Variable = "Stringency\nIndex", filter(posterior_preds, Vax == 0), value = filter(posterior_preds, Vax == 0)$SI),
                cbind(Variable = "Fraction of population\nfully vaccinated (%)", filter(posterior_preds, SI == 0), value = filter(posterior_preds, SI == 0)$Vax)) %>%
  dplyr::select(-SI, -Vax, -id, -.chain, -.iteration) %>%
  pivot_wider(names_from = value, values_from = preds) %>%
  dplyr::select(-.draw) %>%
  data.frame(slope = apply(as.matrix(dplyr::select(., -Variable)), 1, function(y) lm(y ~ x)$coef[2])/100) %>%
  dplyr::select(Variable, slope) %>%
  mutate(Variable = factor(Variable, levels = c("Stringency\nIndex", "Fraction of population\nfully vaccinated (%)"))) %>%
  ggplot(aes(x = slope)) +
  geom_density(linewidth = 0.25, alpha = 0.5, color = "black", fill = "black", adjust = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous("Linearized effect on daily confirmed case change") +
  scale_y_continuous("Density", expand = c(0,0)) +
  facet_wrap(~ Variable, ncol = 2, scales = "free_x") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey93"),
        panel.spacing.x = unit(1, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

FigS4A / FigS4B / FigS4C / FigS4D +
  plot_layout(heights = c(1.5, 1, 1, 1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 10, face = "bold"))
ggsave("Figure_S4.png", width = 16, height = 24, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S5
#-------------------------------------------------------------------------------

parameters <- do.call(cbind, lapply(c("lp__","intercept","beta_country","Omega","tau","mu_mvt","alpha_smooth","eta_smooth","rho_lags","alpha_lags","eta_lags","sigma"), function(i) stanfit[,startsWith(colnames(stanfit), i)]))

# Select 27 random (relevant) parameters (repeat multiple times to comprehensively check for any non-convergence)
data.frame(.chain = rep(1:attributes(stanfit)$nchains, each = attributes(stanfit)$dim[1]/attributes(stanfit)$nchains),
           .iteration = (as.numeric(rownames(stanfit))-1) %% (attributes(stanfit)$dim[1]/attributes(stanfit)$nchains),
           parameters[,sample(1:ncol(parameters), 27, replace = F)]) %>%
  pivot_longer(!c(.chain, .iteration)) %>%
  ggplot() +
  geom_path(aes(x = .iteration, y = value, color = factor(.chain)), linewidth = 0.1) +
  scale_x_continuous("Iteration", expand = c(0,0)) +
  scale_y_continuous("Value") +
  scale_color_brewer("Chain", palette = "Set1", guide = guide_legend(override.aes = list(linewidth = 1))) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        axis.ticks = element_line(size = 0.3),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))
ggsave("Figure_S5.png", width = 16, height = 21, units = "cm", dpi = 600)

rhats <- summarise_draws(stanfit[,colnames(parameters)], "rhat")
hist(rhats$rhat)

#-------------------------------------------------------------------------------
# Figure S6
#-------------------------------------------------------------------------------

# Perform country-level posterior predictive checks
cbind(CountryName = levels(metadata_countries$Country)[rep(datalist$country, each = length(seq(1, nrow(stanfit_gqs), by = 10)))], draw = seq(1, nrow(stanfit_gqs), by = 10), melt(stanfit_gqs[seq(1, nrow(stanfit_gqs), by = 10),startsWith(colnames(stanfit_gqs), "y_rep")])) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_line(data = mutate(data,
                          CountryName = case_when(CountryName == "Bosnia and Herzegovina" ~ "Bosn. and Herz.",
                                                  TRUE ~ CountryName)), aes(x = Y, color = "Observed data", size = "Observed data", alpha = "Observed data"), stat = "density", adjust = 2, bounds = c(-1,1)) +
  geom_line(aes(x = value, group = draw, color = "Posterior predictive simulation", size = "Posterior predictive simulation", alpha = "Posterior predictive simulation"), stat = "density", adjust = 2, bounds = c(-1,1)) +
  coord_cartesian(xlim = c(-0.35,0.35)) +
  scale_x_continuous("Change in daily confirmed cases") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_manual("", values = c("#0f94a8", "black")) +
  scale_alpha_manual("", values = c(1, 0.1)) +
  scale_size_manual("", values = c(1, 0.1)) +
  facet_wrap(~ CountryName, ncol = 5) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        axis.ticks = element_line(size = 0.3),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))
ggsave("Figure_S6.png", width = 16, height = 21, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S7
#-------------------------------------------------------------------------------

# Fetch the posterior samples for runs with an alternative prior specification
stanfit_tightpriors <- readRDS("covid_stringency_effect_fit_tightpriors.rds")
stanfit_widepriors <- readRDS("covid_stringency_effect_fit_widepriors.rds")

# Select 27 random (relevant) parameters
selected_parameters <- do.call(c, sapply(c("lp__","intercept","beta_country","Omega","tau","mu_mvt","alpha_smooth","eta_smooth","rho_lags","alpha_lags","eta_lags","sigma"), function(i) which(startsWith(colnames(stanfit), i))))
selected_parameters <- selected_parameters[sample(1:length(selected_parameters), 27, replace = F)]

rbind(cbind(priors = "Regular", melt(stanfit[,selected_parameters])),
      cbind(priors = "Tighter", melt(stanfit_tightpriors[,selected_parameters])),
      cbind(priors = "Wider", melt(stanfit_widepriors[,selected_parameters]))) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_density(aes(x = value, fill = priors), alpha = 0.3, adjust = 2) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_x_continuous("Value") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_fill_brewer("Prior specification", palette = "Set1") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        axis.ticks = element_line(size = 0.3),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
ggsave("Figure_S7.png", width = 16, height = 21, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S8-10: temporal evolution
#-------------------------------------------------------------------------------

# Rerun script without Belarus, or with temporal effect (temporal = 1)

#-------------------------------------------------------------------------------
# Figure S11: temporal evolution
#-------------------------------------------------------------------------------

Sys.setlocale("LC_ALL", "English")
colMeans(stanfit_gqs[,substr(colnames(stanfit_gqs), 1, 4) == "temp"]) %>%
  melt() %>%
  rownames_to_column() %>%
  separate(rowname, c("day","country"), ",") %>%
  mutate(day = min(data$Date) - 1 + as.numeric(gsub("temporal_effect\\[", "", day)),
         country = as.numeric(gsub("\\]", "", country))) %>%
  mutate(value = case_when(value > 0.1 ~ 0.1,
                           value < -0.1 ~ -0.1,
                           T ~ value)) %>%
  left_join(metadata_countries) %>%
  ggplot() +
  geom_tile(aes(x = day, y = Country_rev, fill = value, color = value)) +
  scale_x_date("Date", expand = c(0,0), date_breaks = "2 months", date_labels = "%b %y") +
  scale_color_gradient2("Predicted change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                        breaks = seq(-0.1, 0.1, by = 0.05), labels = c("<-0.1", 0.05, 0, 0.05, ">0.1"),
                        guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"), ticks.colour = "black")) +
  scale_fill_gradient2("Predicted change in daily confirmed cases", low = "#00756c", mid = "#ede7c2", high = "#c23252",
                       breaks = seq(-0.1, 0.1, by = 0.05), labels = c("<-0.1", 0.05, 0, 0.05, ">0.1"),
                       guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"), ticks.colour = "black")) +
  theme(panel.background = element_blank(),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key = element_blank(),
        legend.position = "bottom")
ggsave("Figure_S11.png", width = 16, height = 18, units = "cm", dpi = 600)
