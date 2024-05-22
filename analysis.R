library(tidyverse)
library(brms)
library(broom)
library(broom.mixed)
options(mc.cores = parallel::detectCores())

dat <- read_csv(here::here("data/PFOS_Dairy_Mice_SHARE_02_28_v2_Revised_v11_CLEANED..csv"),
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
    cens = ifelse(pfos_plasma_mean <= 0.000025, "left", "none"),
    log.pfos.exposure.level.mean = log(pfos_exposure_level_mean),
    log.pfos.exposure.level.mean.lt10 = as.numeric(log.pfos.exposure.level.mean < -10),
    animal2 = fct_collapse(animal, "Rodent" = c("mice", "rat")),
    animal.sex.lactate = fct_cross(animal2, sex, lactating),
    animal.sex.lactate = fct_collapse(animal.sex.lactate, 
                                      "lactating cow/sheep" = c("Cattle:IF:Y", "Sheep:IF:Y"))
  ) |>
  filter(!is.na(pfos_plasma_se)) |>
  rename_with(~ str_replace_all(.x, "_", "."))


# fit1_lactation <- brm(
#   pfos.plasma.mean | 
#     weights(pfos.plasma.ns) + 
#     se(pfos.plasma.se, sigma = TRUE) + 
#     cens(cens) ~ 
#     0 + Intercept +
#     treatment +
#     animal.sex.lactate +
#     potential.exposure.day +
#     # log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
#     s(pfos.exposure.level.mean) +
#     (1 | study.id),
#   family = gaussian(link = "log"),
#   prior = c(
#     prior(normal(5, 15), class = "b", coef = "treatmentCLAEP"),
#     # prior(normal(0, 3), class = "b", coef = "animalmice"),
#     # prior(normal(0, 3), class = "b", coef = "animalrat"),
#     # prior(normal(0, 3), class = "b", coef = "animalSheep"),
#     prior(normal(0, .1), class = "b", coef = "potential.exposure.day"),
#     prior(normal(-5, 15), class = "b", coef = "Intercept")
#     # ,
#     # prior(normal(0, 3), class = "b", coef = "sexIF"),
#     # prior(normal(0, 3), class = "b", coef = "sexIM")
#   ),
#   data = brms_dat,
#   init = 0,
#   iter = 4000,
#   control = list(max_treedepth = 20)
# )
# 
# pairs(fit1_lactation, variable = variables(fit1_lactation)[c(1, 3:9, 14:19)])


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
    # log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
    s(pfos.exposure.level.mean) +
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
  # iter = 4000,
  control = list(max_treedepth = 18)
)

write_rds(fit1, "fit1_res_splines.rds")
# pairs(fit1, variable = variables(fit1)[c(1, 6:7, 20:29)])

# fit1_nosigma <- brm(
#   pfos.plasma.mean | 
#     weights(pfos.plasma.ns) + 
#     # se(pfos.plasma.se) + 
#     cens(cens) ~ 
#     0 + Intercept +
#     treatment +
#     animal +
#     sex +
#     potential.exposure.day +
#     # log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
#     s(pfos.exposure.level.mean) +
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
#   # init = 0,
#   # iter = 4000,
#   control = list(max_treedepth = 18)
# )

fit2 <- brm(
  pfos.plasma.mean | 
    weights(pfos.plasma.ns) + 
    se(pfos.plasma.se, sigma = TRUE) + 
    cens(cens) ~ 
    # 0 + Intercept +
    treatment +
    animal +
    sex +
    potential.exposure.day +
    # log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
    s(pfos.exposure.level.mean) +
    (1 | study.id),
  family = gaussian(link = "log"),
  prior = c(
    prior(horseshoe(1, par_ratio = .5), class = "b")
  ),
  data = brms_dat,
  init = 0,
  iter = 4000,
  control = list(max_treedepth = 18)
)

write_rds(fit2, "fit2_res_splines.rds")

# bind_rows(tidy(fit1), tidy(fit1_nosigma),.id = "fit") |>
#   mutate(across(c(estimate, conf.low, conf.high), ~ ifelse(effect == "fixed", exp(.x), .x))) |>
#   mutate(across(where(is.numeric), scales::number, accuracy = .01, big.mark = ","),
#          est = str_glue("{estimate} ({conf.low}, {conf.high})"))

fit3 <- brm(
  pfos.plasma.mean | 
    weights(pfos.plasma.ns) + 
    se(pfos.plasma.se, sigma = TRUE) + 
    cens(cens) ~ 
    # 0 + Intercept +
    treatment +
    animal +
    sex +
    potential.exposure.day +
    # log.pfos.exposure.level.mean*log.pfos.exposure.level.mean.lt10 +
    s(pfos.exposure.level.mean) +
    (1 | study.id),
  family = gaussian(link = "log"),
  prior = c(
    prior(horseshoe(1, par_ratio = .1), class = "b")
  ),
  data = brms_dat,
  init = 0,
  iter = 4000,
  control = list(max_treedepth = 18)
)

write_rds(fit3, "fit3_res_splines.rds")



