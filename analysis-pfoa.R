library(tidyverse)
library(brms)
library(broom)
library(broom.mixed)
options(mc.cores = parallel::detectCores())

dat <- read_csv(here::here("data/PFOA_Mice_Dataset.csv"),
                na = c("", "."))

brms_dat <- dat |>
  mutate(
    pfoa_plasma_se =
      case_when(
        pfoa_plasma_ns == 1 ~ 0,
        pfoa_plasma_ns > 1 & !is.na(pfoa_plasma_sd) ~ pfoa_plasma_sd / sqrt(pfoa_plasma_ns)
      ),
    animal = str_to_lower(animal),
    treatment = ifelse(treatment %in% c("CONEP", "CONNE"), "CLAEP", treatment),
    treatment = fct_relevel(treatment, "CLANEP"),
    cens = ifelse(pfoa_plasma_mean <= 0.000025, "left", "none")
  ) |>
  filter(!is.na(pfoa_plasma_se)) |>
  rename_with(~ str_replace_all(.x, "_", "."))


fit1 <- brm(
  pfoa.plasma.mean | 
    weights(pfoa.plasma.ns) + 
    se(pfoa.plasma.se, sigma = TRUE) + 
    cens(cens) ~ 
    0 + Intercept +
    treatment +
    animal +
    sex +
    potential.exposure.day +
    s(pfoa.exposure.level.mean) +
    (1 | study.id),
  family = gaussian(link = "log"),
  prior = c(
    prior(normal(5, 15), class = "b", coef = "treatmentCLAEP"),
    prior(normal(0, 3), class = "b", coef = "animalrat"),
    prior(normal(0, .1), class = "b", coef = "potential.exposure.day"),
    prior(normal(-5, 15), class = "b", coef = "Intercept"),
    prior(normal(0, 3), class = "b", coef = "sexIM")
  ),
  data = brms_dat,
  init = 0,
  control = list(max_treedepth = 18)
)

write_rds(fit1, "fit1_res_pfoa.rds")
# pairs(fit1, variable = variables(fit1)[c(1, 6:7, 20:29)])


fit2 <- brm(
  pfoa.plasma.mean | 
    weights(pfoa.plasma.ns) + 
    se(pfoa.plasma.se, sigma = TRUE) + 
    cens(cens) ~ 
    treatment +
    animal +
    sex +
    potential.exposure.day +
    s(pfoa.exposure.level.mean) +
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

write_rds(fit2, "fit2_res_pfoa.rds")

