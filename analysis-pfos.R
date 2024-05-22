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
    cens = ifelse(pfos_plasma_mean <= 0.000025, "left", "none")
  ) |>
  filter(!is.na(pfos_plasma_se)) |>
  rename_with(~ str_replace_all(.x, "_", "."))


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
  control = list(max_treedepth = 18)
)

write_rds(fit1, "fit1_res_pfos.rds")
# pairs(fit1, variable = variables(fit1)[c(1, 6:7, 20:29)])


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

write_rds(fit2, "fit2_res_pfos.rds")


