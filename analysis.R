library(tidyverse)
library(brms)
library(broom)
library(broom.mixed)
options(mc.cores = parallel::detectCores())
options(brms.threads = 4)

dat <- read_csv(here::here("data/PFOS_Dairy_Mice_SHARE_02_28_v2_Revised_v11_CLEANED.csv"),
                na = c("", "."))

brms_dat <- dat |>
  mutate(
    pfos_plasma_se =
      case_when(
        pfos_plasma_ns == 1 ~ 0,
        pfos_plasma_ns > 1 & !is.na(pfos_plasma_sd) ~ pfos_plasma_sd / sqrt(pfos_plasma_ns)
      ),
    treatment = ifelse(treatment %in% c("CONEP", "CONNE"), "CLAEP", treatment),
    treatment = fct_relevel(treatment, "CLANEP"),
    animal2 = ifelse(animal %in% c("rat", "mice"), "micerat", animal),
    sex2 = case_when(
      study_id %in% c(4, 13, 29, 30, 54) ~ "Both",
      sex %in% c("IM", "IF") ~ "Intact",
      .default = sex
    ),
    cens = ifelse(pfos_plasma_mean <= 0.000025, "left", "none"),
    log.pfos.exposure.level.mean = log(pfos_exposure_level_mean),
    log.pfos.exposure.level.mean.lt10 = as.numeric(log.pfos.exposure.level.mean < -10)
  ) |>
  filter(!is.na(pfos_plasma_se)) |>
  rename_with(~ str_replace_all(.x, "_", "."))

# count(brms_dat, sex2)
# 
# brms_dat |>
#   group_by(treatment) |>
#   summarise(mean(pfos.exposure.level.mean), mean(pfos.plasma.mean), n())

# try
# sigma = FALSE
# stronger priors
# log transform first
# fit1 <- brm(
#   pfos.plasma.mean | weights(pfos.plasma.ns) + se(pfos.plasma.se, sigma = TRUE)
#   + cens(cens)
#   ~ 0 + Intercept +
#     treatment +
#     animal +
#     # # sex +
#     potential.exposure.day +
#     s(pfos.exposure.level.mean)  +
#     (1 | study.id),
#   family = gaussian(link = "log"),
#   prior = c(
#     prior(normal(5, 15), class = "b", coef = "treatmentCLAEP"),
#     prior(normal(0, 3), class = "b", coef = "animalmice"),
#     prior(normal(0, 3), class = "b", coef = "animalrat"),
#     prior(normal(0, 3), class = "b", coef = "animalSheep"),
#     rmal(0, .1), class = "b", coef = "potential.exposure.day"),
#     prior(normal(-5, 15), class = "b", coef = "Intercept")
#
#   #           prior(normal(0, 3), class = "b", coef = "sex2IF"),
#   #           prior(normal(0, 3), class = "b", coef = "sex2IM"),
#   # prior(normal(0, 3), class = "b", coef = "sex2CM")
#   ),
#
#   #           prior(normal(0, 50), class = "Intercept")),
#   data = brms_dat,
#   init = 0,
#   # iter = 1000,
#   file = "fit1",
#   file_refit = "on_change"
#   ,
#   control = list(max_treedepth = 18)
# )
#
# pairs(fit1)


# plot(conditional_effects(fit1, "treatment"),
#   ncol = 2,
#   theme = theme_minimal(),
#   plot = FALSE
# )[[1]] +
#   labs(y = "mean PFOS (plasma)")
# 
# tidy(fit1) |>
#   mutate(across(c(estimate, conf.low, conf.high), ~ ifelse(effect == "fixed", exp(.x), .x)))
# 

fit1 <- brm(
  pfos.plasma.mean | 
    weights(pfos.plasma.ns) + 
    se(pfos.plasma.se, sigma = TRUE) + 
    cens(cens) ~ 
    0 + Intercept +
    treatment +
    animal +
    sex +
    potential.exposure.day +
    log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
    # s(pfos.exposure.level.mean) +
    (1 | study.id),
  family = gaussian(link = "log"),
  prior = c(
    prior(normal(5, 15), class = "b", coef = "treatmentCLAEP"),
    prior(normal(0, 3), class = "b", coef = "animalmice"),
    prior(normal(0, 3), class = "b", coef = "animalrat"),
    prior(normal(0, 3), class = "b", coef = "animalSheep"),
    prior(normal(0, .1), class = "b", coef = "potential.exposure.day"),
    prior(normal(-5, 15), class = "b", coef = "Intercept"),
    prior(normal(0, 3), class = "b", coef = "sexIF"),
    prior(normal(0, 3), class = "b", coef = "sexIM")
  ),
  data = brms_dat,
  init = 0,
  control = list(max_treedepth = 18)
)

write_rds(fit1, "fit1_res.rds")



# pairs(fit1)
# 
# 
# fit1 <- brm(
#   pfos.plasma.mean | 
#     weights(pfos.plasma.ns) + 
#     se(pfos.plasma.se, sigma = TRUE) + 
#     cens(cens) ~ 
#     0 + Intercept +
#     # treatment +
#     # animal +
#     # sex +
#     # potential.exposure.day +
#     log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
#     # s(pfos.exposure.level.mean) +
#     (1 | study.id),
#   family = gaussian(link = "log"),
#   prior = c(
#     prior(normal(5, 15), class = "b", coef = "treatmentCLAEP"),
#     prior(normal(0, 3), class = "b", coef = "animalmice"),
#     prior(normal(0, 3), class = "b", coef = "animalrat"),
#     prior(normal(0, 3), class = "b", coef = "animalSheep"),
#     prior(normal(0, .1), class = "b", coef = "potential.exposure.day"),
#     prior(normal(-5, 15), class = "b", coef = "Intercept"),
#     prior(normal(0, 3), class = "b", coef = "sexIF"),
#     prior(normal(0, 3), class = "b", coef = "sexIM")
#   ),
#   data = brms_dat,
#   init = 0
#   # ,
#   # iter = 200,
#   # control = list(max_treedepth = 18)
# )
# 
# pairs(fit1)
# 
# 
# # 
# # brms_dat |>
# #   group_by(sex) |>
# #   summarise(weighted.mean(pfos.plasma.mean, w = pfos.plasma.ns))
# # 
# # select(
# #   brms_dat,
# #   pfos.plasma.mean, pfos.plasma.ns, pfos.plasma.se, treatment, animal,
# #   bw.mean.kg, sex, potential.exposure.day, pfos.exposure.level.mean
# # ) |>
# #   na.omit()
# # 
# 
# # pairs(fit1)
# # stancode(fit1)
# 
# fit2 <- brm(
#     pfos.plasma.mean | 
#       weights(pfos.plasma.ns) + 
#       se(pfos.plasma.se, sigma = TRUE) + 
#       cens(cens) ~ 
#       0 + Intercept +
#       treatment +
#       animal +
#       sex +
#       potential.exposure.day +
#       s(pfos.exposure.level.mean) +
#       (1 | study.id),
#     family = gaussian(link = "log"),
#   prior = c(prior(horseshoe(1, par_ratio = .5), class = "b")),
#   data = brms_dat,
#   init = 0,
#   iter = 5000,
#   control = list(max_treedepth = 20,
#                  adapt_delta = 0.95)
# )
# 
# pairs(fit2)
# 
# bind_rows(tidy(fit1), tidy(fit2), .id = "fit") |>
#   mutate(across(c(estimate, conf.low, conf.high), ~ ifelse(effect == "fixed", exp(.x), .x))) |>
#   mutate(across(where(is.numeric), scales::number)) |>
#   arrange(term)
# 
# plot(
#   conditional_effects(fit2, "treatment",
#     conditions = make_conditions(fit1, c("animal", "sex"))
#   ),
#   ncol = 2,
#   theme = theme_minimal(),
#   plot = FALSE
# )[[1]] +
#   labs(y = "mean PFOS (plasma)")
# 
# plot(conditional_effects(fit2, "treatment"),
#   ncol = 2,
#   theme = theme_minimal(),
#   plot = FALSE
# )[[1]] +
#   labs(y = "mean PFOS (plasma)")
