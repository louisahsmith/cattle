# library(rstan)
library(tidyverse)
library(brms)
library(broom)
library(broom.mixed)
options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

dat <- read_csv(here::here("data/PFOS_Dairy_Mice_SHARE_02_28_v2_Revised_v8.csv"), na = c("", ".")) 
# J1_dat <- filter(dat, pfos_plasma_ns == 1) |> 
#   select(pfos_plasma_mean, pfos_plasma_sd, pfos_plasma_ns)
# J2_dat <- filter(dat, pfos_plasma_ns > 1, !is.na(pfos_plasma_sd)) |> 
#   select(pfos_plasma_mean, pfos_plasma_sd, pfos_plasma_ns) |> 
#   uncount(pfos_plasma_ns)
# 
# stan_dat <- list(J1 = nrow(J1_dat),
#                  J2 = nrow(J2_dat),
#                  ii1 = 1:nrow(J1_dat),
#                  ii2 = nrow(J1_dat) + 1:nrow(J2_dat),
#                  y1 = J1_dat$pfos_plasma_mean,
#                  alpha2 = J2_dat$pfos_plasma_mean,
#                  sigma2 = J2_dat$pfos_plasma_sd)
# 
# fit <- stan(file = here::here("stan/simple-model.stan"), data = stan_dat)
# 
# summary(fit)
# beta         3.951090e+00 4.680504e-03 2.749458e-01 
# tau          7.623863e+00 3.348756e-03 2.000751e-01 

brms_dat <- dat |> 
  mutate(pfos_plasma_se = 
           case_when(pfos_plasma_ns == 1 ~ 0,
                     pfos_plasma_ns > 1 & !is.na(pfos_plasma_sd) ~ pfos_plasma_sd / sqrt(pfos_plasma_ns))) |>
  filter(!is.na(pfos_plasma_se)) |> 
  rename_with(~str_replace_all(.x, "_", "."))

count(brms_dat, animal, treatment)

brms_dat |> 
  group_by(animal, treatment) |> 
  summarise(mean(pfos.plasma.mean))

fit1 <- brm(
  pfos.plasma.mean | weights(pfos.plasma.ns) + se(pfos.plasma.se, sigma = TRUE) ~
    treatment + animal + sex + potential.exposure.day + 
    pfos.exposure.level.mean + bw_mean_kg + (1 | study.id),
  family = gaussian(link = "log"),
  data = brms_dat,
 control = list(adapt_delta = 0.9,
               max_treedepth = 15)
)

# compare CLANEP and CAELP
tidy(fit1) |> 
  mutate(across(c(estimate, conf.low, conf.high), ~ifelse(effect == "fixed", exp(.x), .x)))

plot(conditional_effects(fit1, "treatment", 
                         conditions = make_conditions(fit1, "animal")),
     ncol = 2,
     theme = theme_minimal(),
     plot = FALSE)[[1]] +
  labs(y = "mean PFOS (plasma)")

plot(conditional_effects(fit1, "sex"),
     ncol = 2,
     theme = theme_minimal(),
     plot = FALSE)[[1]] +
  labs(y = "mean PFOS (plasma)")
