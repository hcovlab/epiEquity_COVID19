## ---- Libraries --------------------------------------------------------------
suppressWarnings(suppressMessages({
  library(tidyverse)
  library(data.table)
  library(lubridate)
  library(readxl)
  library(optimx)
  library(writexl)
}))

## ---- Configuration ----------------------------------------------------------
CFG <- list(
  data_dir   = "Data",
  out_dir    = "Results",
  pop_code   = "GBR_NP",   # demographics + life expectancy (Population.csv, as in the analytical model)
  pop_year   = 2018,
  lag_weeks  = 2,          # infection -> death reporting lag for the dose decomposition (eq. 17)
  age_groups = c("40-","40","45","50","55","60","65","70","75","80"),  # vaccination-model groups
  # the UK death data are reported in coarse, irregular age bands (not 5-year);
  # the reconstruction runs at these bands, with coverage/IFR/population aggregated to them
  death_bands = list(band = c("45-64","65-74","75-84","85+"),
                     lo = c(45, 65, 75, 85), hi = c(64, 74, 84, 104)),
  # variant mortality relative risk (death | infection), vs non-VOC baseline
  variant_mort = c(nonVOC = 1.00, Alpha = 1.55),       # Davies 2021 (CI 1.39-1.72)
  # sex-prioritization grid for the counterfactual sweep (s = male relative weight)
  s_grid = seq(0.5, 3.0, by = 0.1),
  # default redistribution (attack-rate) scenario for the main analysis
  incidence_scenario = "quasiEmpirical",
  # vaccination model optimiser: "saved" reuses stored (optimal) parameters -- fast
  # and deterministic; "best" re-fits with L-BFGS-B but is slow through the greedy
  # allocator (Phase 5 would vectorise it). Saved params give the same fit (SSE 0.62).
  vacc_algorithm = "saved"
)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# Literature relative-risk parameters (medians + 95% CI for later CI sampling).
# Vaccine effectiveness is encoded as a relative risk RR = 1 - VE (PHE Week 38
# report, Alpha-dominant period; values unchanged from the original code).
#   condition 1/2 = one/two doses, vs unvaccinated (0)
#   Alpha = B.1.1.7 mortality vs non-VOC (Davies 2021)
LitData <- tibble(
  observedCondition = c("1","1","2","2","Alpha"),
  baseCondition     = c("0","0","0","0","nonVOC"),
  Outcome           = c("Infection","Death","Infection","Death","Death"),
  lower             = c(0.30, 0.15, 0.20, 0.01, 1.39),
  median            = c(0.375,0.275,0.225,0.05, 1.55),
  higher            = c(0.45, 0.40, 0.25, 0.10, 1.72)
)
# convenience accessor for a median relative risk
RR <- function(obs, out, base = "0", col = "median")
  LitData[[col]][LitData$observedCondition == obs & LitData$Outcome == out & LitData$baseCondition == base][1]

## ---- Data loading -----------------------------------------------------------

# Weekly COVID-19 deaths (COVerAGE-DB, UK/GB), by sex. The data come in coarse
# irregular age bands; we keep ages 45+ and map the raw lower-bound labels
# (45, 65, 75, 85) to the named death bands.
load_deaths <- function() {
  db <- CFG$death_bands
  band_of <- setNames(db$band, db$lo)              # 45->"45-64", 65->"65-74", 75->"75-84", 85->"85+"
  fread(file.path(CFG$data_dir, "inputDB_UK.csv")) %>%
    filter(Measure == "Deaths", Sex != "b") %>%
    tibble() %>%
    mutate(Age = as.integer(Age), Date = as.Date(Date, format = "%d/%m/%Y")) %>%
    filter(Age %in% db$lo, Date > as.Date("2020-08-25"), Date < as.Date("2021-05-01")) %>%
    mutate(AgeGroup = factor(band_of[as.character(Age)], levels = db$band)) %>%
    select(Date, Sex, AgeGroup, Value) %>%
    arrange(Sex, AgeGroup, Date) %>%
    group_by(Sex, AgeGroup) %>%
    mutate(Value = Value - lag(Value)) %>%            # cumulative -> weekly increments
    ungroup() %>%
    filter(!is.na(Value), Value >= 0) %>%             # drop the seed week and any revision dips
    mutate(epiYear = epiyear(Date), epiWeek = epiweek(Date)) %>%
    group_by(AgeGroup, epiYear, epiWeek, Sex) %>%
    summarise(Value = sum(Value), .groups = "drop")
}

# O'Driscoll IFR, 5-year groups (drop 0/5/10; 80+ -> 80 label).
load_ifr_raw <- function() {
  read_xlsx(file.path(CFG$data_dir, "IFRdata_ODriscoll.xlsx")) %>%
    filter(!AgeGroup %in% c("0","5","10")) %>%
    mutate(AgeGroup = as.factor(ifelse(AgeGroup == "80+", 80, AgeGroup)))
}

# Demographics and remaining life expectancy from Population.csv (GBR_NP, 2018),
# single-year. Returns single-year and 5-year-group summaries (by sex).
load_population <- function() {
  pop <- fread(file.path(CFG$data_dir, "Population.csv")) %>%
    tibble() %>%
    filter(c_code == CFG$pop_code, Year == CFG$pop_year) %>%
    transmute(Age = as.integer(Age),
              Male, Female, ex.Male, ex.Female)
  # map single years to the 10 vaccination/death age groups
  grp <- function(a) cut(a, breaks = c(-Inf,40,45,50,55,60,65,70,75,80,Inf),
                         labels = CFG$age_groups, right = FALSE)
  pop <- pop %>% mutate(AgeGroup = grp(Age))
  list(single = pop, grp_fun = grp)
}

# Vaccination coverage (dose 1 / dose 2 fractions), and the once-only-vaccinated
# share needed by the allocation model.
load_vaccination <- function() {
  v0 <- read_xlsx(file.path(CFG$data_dir, "VaccinationReports.xlsx")) %>%
    mutate(Dose = as.factor(Dose), Fraction = Fraction / 100)
  vr <- v0 %>%
    filter(Dose == 1) %>%
    left_join(filter(v0, Dose == 2), by = c("epiWeek","AgeGroup")) %>%
    mutate(Dose1 = Fraction.x - Fraction.y) %>%        # exactly-one-dose share
    select(epiYear = epiYear.x, epiWeek, AgeGroup, Dose1, Dose2 = Fraction.y) %>%
    mutate(NotVaccinated = 1 - Dose1 - Dose2)
  list(raw = v0, rates = vr)
}

# ONS weekly incidence (per 10k per day), keyed by infection week.
load_incidence <- function()
  read_xlsx(file.path(CFG$data_dir, "Incidence.xlsx"))

# SGTF proportion among deaths -> Alpha share per epi week (non-SGTF = non-VOC).
load_sgtf <- function(death_weeks) {
  read_xlsx(file.path(CFG$data_dir, "SGTFratiosDeaths.xlsx")) %>%
    right_join(tibble(specimen_date = rep(unique(.$specimen_date), 2),
                      SGTF_label = rep(c("SGTF","Non-SGTF"), length(unique(.$specimen_date)))),
               by = c("specimen_date","SGTF_label")) %>%
    replace(is.na(.), 0) %>%
    mutate(epiWeek = epiweek(specimen_date), epiYear = epiyear(specimen_date)) %>%
    group_by(epiYear, epiWeek, SGTF_label) %>% summarize(N = sum(N), .groups = "drop") %>%
    group_by(epiYear, epiWeek) %>% mutate(AlphaRatio = N / sum(N)) %>% ungroup() %>%
    filter(SGTF_label == "SGTF") %>% select(epiYear, epiWeek, AlphaRatio) %>%
    right_join(death_weeks, by = c("epiYear","epiWeek")) %>%
    mutate(AlphaRatio = replace(AlphaRatio, is.na(AlphaRatio), 0),                 # before data: non-VOC only
           AlphaRatio = replace(AlphaRatio, AlphaRatio == 0 & epiYear == 2021, 1)) %>% # after data: Alpha only
    arrange(epiYear, epiWeek)
}

## ---- Interval IFR (C2): single-year exponential fit -> band averages --------
# Following the analytical model and the manuscript, we fit a 2-parameter
# exponential IFR(A) = exp(a + b*A) per sex over ages 20-104, choosing (a,b) so
# that the population-weighted average of the fitted curve over each O'Driscoll
# 5-year band matches the published band IFR on the log scale. Each target band
# (supplied via `bands`) is then assigned the population-weighted mean of the
# fitted single-year curve over its ages -- a proper open-ended-band average.
# Called for both the coarse death bands (reconstruction) and the finer
# vaccination groups (the IFR-weighted allocation in `ef`).
build_interval_ifr <- function(ifr_raw, pop_single, bands) {
  age_min <- 20L; age_max <- 104L
  bind_rows(lapply(c("m","f"), function(sx) {
    popcol <- if (sx == "m") "Male" else "Female"
    pw <- pop_single %>% filter(Age >= age_min, Age <= age_max) %>%
      transmute(Age, w = .data[[popcol]])
    band_avg <- function(a, b, lo, hi) {
      ww <- pw %>% filter(Age >= lo, Age <= hi); sum(exp(a + b * ww$Age) * ww$w) / sum(ww$w)
    }
    d <- ifr_raw %>% filter(Sex == sx) %>%
      mutate(a0 = as.integer(as.character(AgeGroup))) %>%
      filter(!is.na(a0), a0 >= age_min) %>% arrange(a0)
    band_lo <- d$a0; band_hi <- ifelse(d$a0 >= 80, age_max, d$a0 + 4L); odr <- d$IFRmedian
    sse <- function(p) sum((log(mapply(band_avg, p[1], p[2], band_lo, band_hi)) - log(odr))^2)
    fit <- optim(c(log(odr[1]) - 0.1 * band_lo[1], 0.1), sse)   # Nelder-Mead
    grp_ifr <- mapply(function(lo, hi) band_avg(fit$par[1], fit$par[2], lo, hi), bands$lo, bands$hi)
    tibble(AgeGroup = factor(bands$band, levels = bands$band), Sex = sx, IFR = as.numeric(grp_ifr))
  }))
}
## vaccination-group age ranges (single years), for the IFR-weighted allocation
VACC_BANDS <- list(band = CFG$age_groups,
                   lo = c(20,40,45,50,55,60,65,70,75,80), hi = c(39,44,49,54,59,64,69,74,79,104))

## ---- Vaccination campaign model (manuscript eq. 10-16) ------------------------
# Greedy weekly dose allocation across the 10 age groups, driven by age-group
# weights and a sigmoidal uptake-saturation term. Calibrated to observed weekly
# coverage. predSex modes: "parameter" (scalar male weight s), "step" (males one
# band earlier), "IFR" (male/female IFR ratio as the weight).
ef <- function(par, data, pop, algMode = "optim", predSex = FALSE,
               predSexMethod = "parameter", predSexPar = 1) {
  if (algMode == "optim" && predSex) stop("Do not optimise with predSex=TRUE.")
  par_NV <- par[1:4]
  if (predSexMethod != "IFR") {
    par_agegroups <- par[5:14]
  } else {
    male_IFR   <- IFR_GROUP_VEC$m
    female_IFR <- IFR_GROUP_VEC$f
    par_agegroups <- female_IFR
  }
  if (predSex) {
    if (predSexMethod == "parameter")  s <- c(rep(1,10), rep(predSexPar,10))
    else if (predSexMethod == "IFR")   s <- c(rep(1,10), IFR_GROUP_VEC$m / IFR_GROUP_VEC$f)
    basepop <- arrange(pop, Sex, AgeGroup)$Population
  } else {
    pop <- pop %>% group_by(AgeGroup) %>% summarize(Population = sum(Population), .groups = "drop")
    basepop <- arrange(pop, AgeGroup)$Population
  }
  nv1 <- basepop; nv2 <- basepop
  if (algMode == "predict") {
    VaccinationModeled <- data %>% select(epiYear, epiWeek, Dose, DosesAdministeredWeek)
    vacc_modeled <- data %>% ungroup() %>%
      select(-epiWeek, -epiYear, -Dose, -DosesAdministeredWeek) %>% filter(`40-` < 0)
  }
  sse <- 0
  for (i in 1:nrow(data)) {
    avail_week <- data$DosesAdministeredWeek[i]
    for (d in 1:7) {
      availableDoses <- avail_week / 7
      while (availableDoses > 1) {
        nv <- if (data$Dose[i] == "1") nv1 else nv2
        NV <- par_NV[3] + (-par_NV[3]) / (1 + ((nv / basepop) / par_NV[2])^par_NV[1])^par_NV[4]
        if (predSex) {
          if (predSexMethod %in% c("parameter","IFR")) {
            rel <- (s * par_agegroups * NV * basepop) / sum(s * par_agegroups * NV * basepop)
          } else if (predSexMethod == "step") {
            pa <- c(par_agegroups, par_agegroups[2:10], par_agegroups[10])  # males one band earlier
            rel <- (pa * NV * basepop) / sum(pa * NV * basepop)
          }
        } else {
          rel <- (par_agegroups * NV * basepop) / sum(par_agegroups * NV * basepop)
        }
        wk <- numeric(length(rel))
        for (j in seq_along(rel)) {
          give <- min(nv[j], rel[j] * availableDoses)
          wk[j] <- give; nv[j] <- nv[j] - give
        }
        if (data$Dose[i] == "1") nv1 <- nv else nv2 <- nv
        if (algMode == "optim" && d == 7)
          sse <- sse + sum(((1 - (nv / basepop)) - data[i, 4:13])^2)
        availableDoses <- availableDoses - sum(wk)
        if (sum(wk) <= 1e-9) break   # all eligible groups saturated; stop (prevents infinite loop)
      }
    }
    if (algMode == "predict") {
      done <- if (data$Dose[i] == "1") (basepop - nv1) / basepop else (basepop - nv2) / basepop
      vacc_modeled <- rbind(vacc_modeled, done)
    }
  }
  if (algMode == "optim") return(sse)
  cn <- if (predSex)
    c(paste0("f.", CFG$age_groups), paste0("m.", CFG$age_groups)) else CFG$age_groups
  colnames(vacc_modeled) <- cn
  cbind(VaccinationModeled, vacc_modeled)
}

## ---- Build all inputs -------------------------------------------------------
message("Loading data ...")
deathCounts <- load_deaths()
IFRraw      <- load_ifr_raw()
POP         <- load_population()
VAC         <- load_vaccination()
IncidenceRate <- load_incidence()
death_weeks <- deathCounts %>% distinct(epiYear, epiWeek)
SGTFdeaths  <- load_sgtf(death_weeks)

# population by VACCINATION group and sex (for the campaign model)
popCountsSex <- POP$single %>%
  group_by(AgeGroup) %>%
  summarize(m = sum(Male), f = sum(Female), .groups = "drop") %>%
  pivot_longer(c(m, f), names_to = "Sex", values_to = "Population") %>%
  mutate(AgeGroup = factor(AgeGroup, levels = CFG$age_groups))

# DEATH-BAND demographics (4 coarse bands) for the reconstruction + outcomes
db <- CFG$death_bands
vacc_of_age <- function(a) ifelse(a < 40, "40-", ifelse(a < 45, "40",
                                                        ifelse(a >= 80, "80", as.character((a %/% 5) * 5))))
POP_death <- POP$single %>%
  mutate(DBand = cut(Age, breaks = c(db$lo, 105), labels = db$band, right = FALSE)) %>%
  filter(!is.na(DBand))
popCountsDeath <- POP_death %>% group_by(DBand) %>%
  summarize(m = sum(Male), f = sum(Female), .groups = "drop") %>%
  pivot_longer(c(m, f), names_to = "Sex", values_to = "Population") %>%
  transmute(AgeGroup = factor(DBand, levels = db$band), Sex, Population)
lifeExp <- POP_death %>% group_by(DBand) %>%
  summarize(m = sum(ex.Male * Male) / sum(Male), f = sum(ex.Female * Female) / sum(Female), .groups = "drop") %>%
  pivot_longer(c(m, f), names_to = "Sex", values_to = "Lx") %>%
  transmute(AgeGroup = factor(DBand, levels = db$band), Sex, Lx)

# interval IFR at both resolutions: death bands (reconstruction) and vacc groups (ef allocation)
IFR_GROUP <- build_interval_ifr(IFRraw, POP$single, db)
IFR_VACC  <- build_interval_ifr(IFRraw, POP$single, VACC_BANDS)
IFR_GROUP_VEC <- list(
  m = (IFR_VACC %>% filter(Sex == "m") %>% arrange(AgeGroup))$IFR,
  f = (IFR_VACC %>% filter(Sex == "f") %>% arrange(AgeGroup))$IFR
)

# population-weighted weights to aggregate vaccination-group coverage to death bands
.cw_build <- function(by_sex) {
  base <- POP$single %>%
    mutate(DBand = cut(Age, breaks = c(db$lo, 105), labels = db$band, right = FALSE),
           vg = vacc_of_age(Age)) %>% filter(!is.na(DBand))
  if (by_sex) {
    base %>% pivot_longer(c(Male, Female), names_to = "Sex", values_to = "pop") %>%
      mutate(Sex = recode(Sex, Male = "m", Female = "f")) %>%
      group_by(Sex, DBand, vg) %>% summarize(pop = sum(pop), .groups = "drop") %>%
      group_by(Sex, DBand) %>% mutate(weight = pop / sum(pop)) %>% ungroup() %>%
      transmute(AgeGroup = factor(DBand, levels = db$band), vacc_group = vg, Sex, weight)
  } else {
    base %>% mutate(pop = Male + Female) %>%
      group_by(DBand, vg) %>% summarize(pop = sum(pop), .groups = "drop") %>%
      group_by(DBand) %>% mutate(weight = pop / sum(pop)) %>% ungroup() %>%
      transmute(AgeGroup = factor(DBand, levels = db$band), vacc_group = vg, weight)
  }
}
COV_WEIGHTS <- .cw_build(TRUE); COV_WEIGHTS_ALL <- .cw_build(FALSE)

## ---- Vaccination model inputs (doses administered per week) ------------------
popCountsNosex <- popCountsSex %>% group_by(AgeGroup) %>% summarize(Population = sum(Population), .groups = "drop")
VaccinationReal2 <- popCountsNosex %>%
  right_join(VAC$rates, by = "AgeGroup") %>%
  mutate(D1 = round(Population * Dose1 + Population * Dose2),
         D2 = round(Population * Dose2),
         AgeGroup = factor(AgeGroup, levels = CFG$age_groups)) %>%
  group_by(AgeGroup) %>%
  mutate(D1w = D1 - lag(D1, default = 0), D2w = D2 - lag(D2, default = 0)) %>%
  group_by(epiYear, epiWeek) %>%
  summarize(`1` = sum(D1w), `2` = sum(D2w), .groups = "drop") %>%
  pivot_longer(c(`1`,`2`), names_to = "Dose", values_to = "DosesAdministeredWeek")
VaccinationReal <- VAC$raw %>%
  mutate(AgeGroup = factor(AgeGroup, levels = CFG$age_groups)) %>%
  spread(AgeGroup, Fraction) %>%                                  # one column per age group
  left_join(VaccinationReal2, by = c("epiYear","epiWeek","Dose")) %>%
  arrange(Dose, epiYear, epiWeek)

## ---- Vaccination model fit --------------------------------------------------
firstguess <- c(6.88,0.189,1.06,1.822, 0.18,0.57,0.75,1.11,1.23,1.60,2.68,4.80,8.28,20.61)
saved_par  <- c(5.473907,0.2154104,1.06,1.962781,0.1808279,0.5863161,0.7705115,
                1.140127,1.281896,1.656613,2.777422,5.079028,8.469866,20.45464)
if (CFG$vacc_algorithm == "best") {
  message("Fitting vaccination model (L-BFGS-B) ...")
  fit <- optimx(firstguess, ef, method = "L-BFGS-B",
                data = VaccinationReal, pop = popCountsSex, lower = 1e-100)
  vacc_par <- as.numeric(fit[1, 1:14])
} else {
  vacc_par <- saved_par
}

# cache modeled coverage (sex-neutral) once; sex-specific is parameter-dependent
predict_coverage <- function(s_param = 1, method = "parameter", sex = FALSE)
  ef(vacc_par, data = VaccinationReal, pop = popCountsSex, algMode = "predict",
     predSex = sex, predSexMethod = method, predSexPar = s_param)

cat("Phase 0 complete: data, interval IFR, life expectancy and vaccination model ready.\n")

# =============================================================================
# Phase 1 -- reconstruction engine (manuscript eq. 17-23), non-VOC/Alpha
# =============================================================================

## relative risks (death | infection) and (infection | exposure) by dose
R_death <- c(`0` = 1, `1` = RR("1","Death"),     `2` = RR("2","Death"))      # 1, 0.275, 0.05
R_inf   <- c(`0` = 1, `1` = RR("1","Infection"), `2` = RR("2","Infection"))  # 1, 0.375, 0.225
R_voc   <- CFG$variant_mort                                                   # nonVOC 1, Alpha 1.55

## epi-week serial arithmetic (2020 has 53 epi weeks) for the t-2 lag
wser  <- function(y, w) (y - 2020) * 53 + w
wyear <- function(s) 2020 + (s - 1) %/% 53
wweek <- function(s) s - (wyear(s) - 2020) * 53

## Coverage -> dose-status shares (0/1/2) per DEATH BAND, in long form.
## Vaccination-group coverage is aggregated to the coarse death bands using
## population weights (sex-specific when sex = TRUE).
coverage_shares <- function(vm, sex) {
  ag <- CFG$age_groups
  long <- vm %>%
    select(epiYear, epiWeek, Dose, all_of(if (sex) c(paste0("f.",ag), paste0("m.",ag)) else ag)) %>%
    pivot_longer(-c(epiYear, epiWeek, Dose), names_to = "key", values_to = "cov")
  if (sex) long <- long %>% mutate(Sex = substr(key, 1, 1), vacc_group = substr(key, 3, nchar(key)))
  else     long <- long %>% mutate(Sex = NA_character_, vacc_group = key)
  vw <- long %>% select(epiYear, epiWeek, Sex, vacc_group, Dose, cov) %>%
    pivot_wider(names_from = Dose, values_from = cov, names_prefix = "v")   # v1, v2 per vacc group
  agg <- if (sex) inner_join(vw, COV_WEIGHTS, by = c("vacc_group","Sex"), relationship = "many-to-many")
  else     inner_join(vw, COV_WEIGHTS_ALL, by = "vacc_group", relationship = "many-to-many")
  agg %>% group_by(epiYear, epiWeek, Sex, AgeGroup) %>%
    summarize(v1 = sum(v1 * weight), v2 = sum(v2 * weight), .groups = "drop") %>%
    transmute(epiYear, epiWeek, Sex, AgeGroup,
              `0` = 1 - v1, `1` = v1 - v2, `2` = v2) %>%      # dose-status shares per death band
    pivot_longer(c(`0`,`1`,`2`), names_to = "Dose", values_to = "share") %>%
    mutate(share = pmax(share, 0), wk = wser(epiYear, epiWeek))
}

## default reconstruction parameters (relative risks + interval IFR); sensitivity
## analyses pass alternative values through the same engine
PARAMS <- list(Ri = R_inf, Rd = R_death, Rv = R_voc, ifr = IFR_GROUP)

## ---- Baseline exposure (adjusted infections), reconstructed once ------------
# Decompose observed deaths by dose (eq. 17) and variant (SGTF), then invert
# eq. 18-19 to adjusted infections (exposure). Indexed by infection week =
# death week - lag. cov_offset shifts the vaccination timeline relative to the
# epidemic (weeks) for the timing sensitivity (Figure S10).
build_baseline_exposure <- function(P = PARAMS, cov_offset = 0L) {
  cov_neutral <- coverage_shares(predict_coverage(), sex = FALSE) %>%
    transmute(wk_inf = wk + cov_offset, AgeGroup, Dose, share)   # shares at the infection week
  
  # eq. 17 dose distribution among deaths, evaluated at the infection week.
  # NOTE: deaths arise from the *exposed* population, so the correct dose weight
  # is the death risk *given exposure* = R_inf * R_death (infection-selection then
  # death-selection). Weighting by R_death alone (as eq. 17 is written) omits the
  # infection-selection step and over-counts observed deaths by ~33%; including
  # R_inf makes the back- and forward-calculations consistent so the sex-neutral
  # baseline reproduces observed deaths exactly (see manuscript note on eq. 17).
  vaccRatio_D <- cov_neutral %>%
    mutate(w = share * P$Ri[Dose] * P$Rd[Dose]) %>%
    group_by(wk_inf, AgeGroup) %>%
    mutate(vaccRatio_D = w / sum(w)) %>% ungroup() %>%
    select(wk_inf, AgeGroup, Dose, vaccRatio_D)
  
  variant <- SGTFdeaths %>%
    transmute(wk_death = wser(epiYear, epiWeek), nonVOC = 1 - AlphaRatio, Alpha = AlphaRatio) %>%
    pivot_longer(c(nonVOC, Alpha), names_to = "VOC", values_to = "vRatio")
  
  deathCounts %>%
    transmute(wk_death = wser(epiYear, epiWeek), AgeGroup, Sex, Dobs = Value) %>%
    mutate(wk_inf = wk_death - CFG$lag_weeks) %>%
    inner_join(vaccRatio_D, by = c("wk_inf","AgeGroup"), relationship = "many-to-many") %>%
    inner_join(variant, by = "wk_death", relationship = "many-to-many") %>%
    left_join(P$ifr, by = c("AgeGroup","Sex")) %>%
    mutate(
      D   = Dobs * vaccRatio_D * vRatio,                                   # deaths by dose x variant
      adI = D / (P$Rd[Dose] * P$Rv[VOC] * IFR * P$Ri[Dose])               # eq. 18 + 19 -> exposure
    ) %>%
    group_by(wk_inf, AgeGroup, Sex, VOC, IFR) %>%
    summarise(exposure = sum(adI), .groups = "drop")                      # sum over dose
}

## ---- Project deaths under an allocation -------------------------------------
# Redistribute the fixed baseline exposure into dose groups using a counterfactual
# coverage allocation, then re-apply eq. 18-19 forward to obtain deaths.
project_deaths <- function(exposure, cov_shares, P = PARAMS) {
  cov <- cov_shares %>% select(wk_inf = wk, Sex, AgeGroup, Dose, share)
  joined <- if (all(is.na(cov$Sex))) {
    exposure %>% inner_join(select(cov, -Sex), by = c("wk_inf","AgeGroup"),
                            relationship = "many-to-many")
  } else {
    exposure %>% inner_join(cov, by = c("wk_inf","AgeGroup","Sex"),
                            relationship = "many-to-many")
  }
  joined %>%
    mutate(D = exposure * share * P$Ri[Dose] * P$Rd[Dose] * P$Rv[VOC] * IFR) %>%
    group_by(AgeGroup, Sex) %>% summarise(Deaths = sum(D), .groups = "drop")
}

## no-vaccination null: everyone unvaccinated (dose 0), so coverage share = 1 at dose 0
null_allocation <- function(exposure) {
  exposure %>% distinct(wk_inf, AgeGroup, Sex) %>%
    crossing(Dose = c("0","1","2")) %>%
    mutate(share = ifelse(Dose == "0", 1, 0), wk = wk_inf, epiYear = wyear(wk), epiWeek = wweek(wk))
}

# =============================================================================
# Phase 2 -- outcome functionals and per-strategy parity solving (eq. 20-23)
# =============================================================================

## per-sex population (modeled adult range = the death bands) and life expectancy
.Nsex  <- popCountsDeath %>% group_by(Sex) %>% summarise(N = sum(Population), .groups = "drop")
.lifeL <- lifeExp

# eq. 20-23 for one allocation's deaths-by-age-sex, given the null deaths D0
outcomes <- function(deaths, D0) {
  deaths %>%
    left_join(.lifeL, by = c("AgeGroup","Sex")) %>%
    left_join(D0 %>% rename(Deaths0 = Deaths), by = c("AgeGroup","Sex")) %>%
    left_join(.Nsex, by = "Sex") %>%
    group_by(Sex) %>%
    summarise(D = sum(Deaths),
              M = sum(Deaths) / first(N),
              Y = sum(Deaths * Lx) / first(N),
              P = 1 - sum(Deaths) / sum(Deaths0),
              .groups = "drop")
}

## sweep the sex-prioritization parameter s and tabulate the metrics + ratios.
## cov_offset shifts the vaccination timeline (Figure S10); P carries RR/IFR.
strategy_sweep <- function(exposure, s_grid = CFG$s_grid, method = "parameter",
                           P = PARAMS, cov_offset = 0L) {
  shift <- function(cv) mutate(cv, wk = wk + cov_offset)
  D0    <- project_deaths(exposure, coverage_shares0(null_allocation(exposure)), P)
  base_cov <- shift(coverage_shares(predict_coverage(), sex = FALSE))
  Dbase <- project_deaths(exposure, base_cov, P)
  ob <- outcomes(Dbase, D0) %>% mutate(s = NA_real_, scenario = "sex-neutral baseline")
  
  rows <- lapply(s_grid, function(s) {
    cov <- shift(coverage_shares(predict_coverage(s_param = s, method = method, sex = TRUE), sex = TRUE))
    outcomes(project_deaths(exposure, cov, P), D0) %>% mutate(s = s, scenario = "counterfactual")
  })
  res <- bind_rows(ob, rows)
  
  ratios <- res %>% filter(scenario == "counterfactual") %>%
    select(s, Sex, D, M, Y, P) %>%
    pivot_wider(names_from = Sex, values_from = c(D, M, Y, P)) %>%
    transmute(s,
              EqualMortality_ratio       = M_m / M_f,   # EM:  M_m = M_f
              EqualYLL_ratio             = Y_m / Y_f,   # EYS: Y_m = Y_f
              EqualDeathReduction_ratio  = P_m / P_f,   # EDR: P_m = P_f
              dDeath_overall = (D_m + D_f) / (ob$D[ob$Sex=="m"] + ob$D[ob$Sex=="f"]) - 1,
              dDeath_m = D_m / ob$D[ob$Sex=="m"] - 1,
              dDeath_f = D_f / ob$D[ob$Sex=="f"] - 1)
  list(outcomes = res, ratios = ratios, baseline = ob)
}

## allocation helper that bypasses the v1/v2 reshape (null is already long shares)
coverage_shares0 <- function(alloc_long)
  alloc_long %>% transmute(epiYear, epiWeek, Sex, AgeGroup = factor(AgeGroup, levels = CFG$death_bands$band),
                           Dose, share, wk)

## solve for the prioritization s* achieving each parity (linear interpolation of the crossing)
solve_strategies <- function(ratios) {
  cross <- function(x, y) {           # first sign change of (y - 1) in x
    f <- y - 1; idx <- which(f[-1] * f[-length(f)] < 0)
    if (!length(idx)) return(NA_real_)
    i <- idx[1]; x[i] - f[i] * (x[i+1] - x[i]) / (f[i+1] - f[i])
  }
  tibble(
    strategy = c("Equal Mortality","Equal YLL","Equal Death Reduction"),
    s_star   = c(cross(ratios$s, ratios$EqualMortality_ratio),
                 cross(ratios$s, ratios$EqualYLL_ratio),
                 cross(ratios$s, ratios$EqualDeathReduction_ratio))
  )
}

cat("Phase 1-2 functions loaded.\n")


# =============================================================================
# Phase 3-5 -- figures (SVG, analytical-model styling), sensitivity, orchestration
# =============================================================================
suppressWarnings(suppressMessages({ library(patchwork); library(scales) }))

## ---- shared styling (identical palette to the analytical model) -------------
IFR_SEX_LIGHT <- c(m = "#9CD0EE", f = "#F6B6CE")
IFR_SEX_DARK  <- c(m = "#15539e", f = "#b32219")
HEAT_BLUE <- "#2166AC"; HEAT_RED <- "#B2182B"; STRAT_ACCENT <- "#222222"; NEUTRAL_GREY <- "#555555"
# strategy colours match the analytical STRAT_TINT (eqDR/eqYPLL/eqPD/eqRisk)
STRAT_COL <- c("Equal Mortality" = "#1A7DB6", "Equal YLL" = "#6A51A3",
               "Equal Death Reduction" = "#D7642A", "Equal IFR Threshold" = "#2C8C5A")
pub_theme <- function(base_size = 19) {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = "grey15"),
          plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
          axis.title = element_text(face = "bold", size = base_size + 1, colour = "grey20"),
          axis.text  = element_text(size = base_size - 2, colour = "grey35"),
          axis.line  = element_line(colour = "grey40", linewidth = 0.7),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.7),
          strip.text = element_text(face = "bold", size = base_size, colour = "grey15"),
          strip.background = element_rect(fill = "grey93", colour = NA),
          legend.title = element_text(face = "bold", size = base_size - 3),
          legend.text  = element_text(size = base_size - 4),
          legend.key = element_blank(), legend.position = "top",
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          plot.margin = margin(10, 12, 10, 10))
}
age_scale <- function(name = "Age") scale_colour_viridis_d(option = "C", end = 0.92, name = name)
SEX_LAB <- c(f = "Female", m = "Male")

## ---- date + vaccination-age-difference helpers ------------------------------
wser_to_date <- function(s) as.Date("2019-12-29") + (s - 1) * 7
ANALYSIS_BOX <- tibble(xmin = wser_to_date(wser(2020, 50)), xmax = wser_to_date(wser(2021, 17)))

# all-age weekly deaths by sex (whole-population panel of Figure 5)
deathsTotal <- fread(file.path(CFG$data_dir, "inputDB_UK.csv")) %>%
  filter(Measure == "Deaths", Sex != "b") %>% tibble() %>%
  mutate(Age = as.integer(Age), Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  filter(Date > as.Date("2020-08-25"), Date < as.Date("2021-05-01")) %>%
  arrange(Sex, Age, Date) %>% group_by(Sex, Age) %>% mutate(V = Value - lag(Value)) %>% ungroup() %>%
  filter(!is.na(V), V >= 0) %>%
  group_by(epiYear = epiyear(Date), epiWeek = epiweek(Date), Sex) %>% summarise(Deaths = sum(V), .groups = "drop")

# mean age within each vaccination group (per sex) and population, for age-of-vaccination
AGE_MID <- POP$single %>% group_by(AgeGroup) %>%
  summarise(m = sum(Age * Male) / sum(Male), f = sum(Age * Female) / sum(Female), .groups = "drop") %>%
  pivot_longer(c(m, f), names_to = "Sex", values_to = "mid") %>% mutate(AgeGroup = as.character(AgeGroup))
POPL <- popCountsSex %>% mutate(AgeGroup = as.character(AgeGroup)) %>% left_join(AGE_MID, by = c("AgeGroup","Sex"))

# weekly mean age of the vaccinated population by sex AND dose, for prioritization s
mean_vacc_age <- function(s = 1, method = "parameter") {
  ag <- CFG$age_groups
  predict_coverage(s_param = s, method = method, sex = TRUE) %>% filter(Dose %in% c("1","2")) %>%
    select(epiYear, epiWeek, Dose, all_of(c(paste0("m.",ag), paste0("f.",ag)))) %>%
    pivot_longer(-c(epiYear, epiWeek, Dose), names_to = "key", values_to = "cov") %>%
    mutate(Sex = substr(key, 1, 1), AgeGroup = substr(key, 3, nchar(key))) %>%
    left_join(POPL, by = c("AgeGroup","Sex")) %>% mutate(nv = cov * Population) %>%
    group_by(epiYear, epiWeek, Dose, Sex) %>% summarise(mage = sum(mid * nv) / sum(nv), .groups = "drop") %>%
    pivot_wider(names_from = Sex, values_from = mage) %>%
    transmute(wk = wser(epiYear, epiWeek), date = wser_to_date(wk),
              Dose = paste0("Dose ", Dose), age_diff = m - f)
}
CAMPAIGN_WKS <- wser(2020, 50):wser(2021, 20)
# campaign-mean age difference for a given s, per dose ("1"/"2")
age_diff_at <- function(s, dose = "1", method = "parameter") {
  p <- mean_vacc_age(s, method) %>%
    filter(Dose == paste0("Dose ", dose), wk %in% CAMPAIGN_WKS, !is.nan(age_diff))
  if (!nrow(p)) NA_real_ else mean(p$age_diff)
}
# fast lookup over a precomputed s -> age-difference grid (both doses; filled in the run)
AGE_DIFF_GRID <- NULL
age_diff_lookup <- function(s, dose = "1")
  approx(AGE_DIFF_GRID$s, AGE_DIFF_GRID[[if (dose == "2") "ad2" else "ad1"]], xout = s, rule = 2)$y

## first s where a m/f ratio curve crosses 1 (linear interpolation)
.cross <- function(x, y) {
  f <- y - 1; idx <- which(f[-1] * f[-length(f)] < 0)
  if (!length(idx)) return(NA_real_)
  i <- idx[1]; x[i] - f[i] * (x[i+1] - x[i]) / (f[i+1] - f[i])
}

## ---- S9 attack-rate redistribution (equalAgeSexONS removed) ------------------
redistribute_exposure <- function(EXP, scenario) {
  if (scenario == "quasiEmpirical") return(EXP)
  popA <- popCountsDeath %>% mutate(AgeGroup = as.character(AgeGroup))
  E <- EXP %>% mutate(AgeGroup = as.character(AgeGroup))
  if (scenario == "equalSex") {
    out <- E %>% group_by(wk_inf, AgeGroup, VOC) %>% mutate(tot = sum(exposure)) %>% ungroup() %>%
      left_join(popA, by = c("AgeGroup","Sex")) %>%
      group_by(wk_inf, AgeGroup, VOC) %>% mutate(exposure = tot * Population / sum(Population)) %>% ungroup()
  } else if (scenario == "equalAgeSex") {
    out <- E %>% group_by(wk_inf, VOC) %>% mutate(tot = sum(exposure)) %>% ungroup() %>%
      left_join(popA, by = c("AgeGroup","Sex")) %>%
      group_by(wk_inf, VOC) %>% mutate(exposure = tot * Population / sum(Population)) %>% ungroup()
  }
  out %>% transmute(wk_inf, AgeGroup = factor(AgeGroup, levels = CFG$death_bands$band), Sex, VOC, IFR, exposure)
}

## ---- cached coverage + fast (optionally time-shifted) sweep -----------------
cache_coverage <- function(s_grid = CFG$s_grid, method = "parameter")
  list(neutral = coverage_shares(predict_coverage(), sex = FALSE),
       per_s   = lapply(s_grid, function(s)
         coverage_shares(predict_coverage(s_param = s, method = method, sex = TRUE), sex = TRUE)),
       s_grid  = s_grid)
sweep_cached <- function(exposure, cc, P = PARAMS, offset = 0L) {
  D0 <- project_deaths(exposure, coverage_shares0(null_allocation(exposure)), P)
  shift <- if (offset == 0L) function(cov) cov else function(cov) pad_and_shift(cov, offset)
  rows <- Map(function(cov, s) outcomes(project_deaths(exposure, shift(cov), P), D0) %>% mutate(s = s),
              cc$per_s, cc$s_grid)
  bind_rows(rows) %>% select(s, Sex, M, Y, P) %>%
    pivot_wider(names_from = Sex, values_from = c(M, Y, P)) %>%
    transmute(s, EM = M_m / M_f, EY = Y_m / Y_f, ED = P_m / P_f)
}
solve_from <- function(rat)
  c(`Equal Mortality` = .cross(rat$s, rat$EM), `Equal YLL` = .cross(rat$s, rat$EY),
    `Equal Death Reduction` = .cross(rat$s, rat$ED))

## ---- S8 parameter sampling --------------------------------------------------
draw_RR <- function(obs, out, base = "0") {
  r <- LitData[LitData$observedCondition == obs & LitData$Outcome == out & LitData$baseCondition == base, ]
  runif(1, r$lower, r$higher)
}
sample_params <- function()
  list(Ri  = c(`0` = 1, `1` = draw_RR("1","Infection"), `2` = draw_RR("2","Infection")),
       Rd  = c(`0` = 1, `1` = draw_RR("1","Death"),     `2` = draw_RR("2","Death")),
       Rv  = c(nonVOC = 1, Alpha = draw_RR("Alpha","Death","nonVOC")),
       ifr = build_interval_ifr(IFRraw %>% mutate(IFRmedian = runif(n(), IFRlowCI, IFRhighCI)),
                                POP$single, CFG$death_bands))

## time-shift a coverage table with unvaccinated padding before / held value after
pad_and_shift <- function(cov, offset) {
  wkr <- range(cov$wk); ag <- unique(cov$AgeGroup); sx <- unique(cov$Sex)
  pre  <- expand_grid(wk = (wkr[1]-12):(wkr[1]-1), AgeGroup = ag, Sex = sx, Dose = c("0","1","2")) %>%
    mutate(share = ifelse(Dose == "0", 1, 0))
  post <- expand_grid(wk = (wkr[2]+1):(wkr[2]+12)) %>%
    crossing(cov %>% filter(wk == wkr[2]) %>% select(AgeGroup, Sex, Dose, share))
  bind_rows(select(cov, wk, AgeGroup, Sex, Dose, share), pre, post) %>% mutate(wk = wk + offset)
}


# small helper: a white-background in-plot legend at a corner
.inleg <- function(x, y, dir = "vertical")
  theme(legend.position = c(x, y), legend.justification = c(x, y), legend.direction = dir,
        legend.background = element_rect(fill = "white", colour = NA), legend.margin = margin(2,4,2,4))

## ===== Figure 5: the UK wave + vaccination-model fit (2 x 2) =================
fig5 <- function() {
  xr  <- as.Date(c("2020-09-01", "2021-05-20"))
  xsc <- scale_x_date(limits = xr, date_breaks = "2 months", date_labels = "%b %y")
  box <- geom_rect(data = ANALYSIS_BOX, inherit.aes = FALSE,
                   aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey88", alpha = 0.5)
  dks <- deathsTotal %>% mutate(date = wser_to_date(wser(epiYear, epiWeek)))
  pTL <- ggplot(dks, aes(date, Deaths, colour = Sex)) + box +
    geom_line(linewidth = 1.1) + geom_point(size = 1.3) +
    scale_colour_manual(values = IFR_SEX_DARK, labels = SEX_LAB) +
    labs(x = NULL, y = "Weekly deaths (all ages)", colour = NULL) + xsc + pub_theme(17) + .inleg(0.03, 0.97)
  vb <- SGTFdeaths %>% mutate(date = wser_to_date(wser(epiYear, epiWeek)))
  pBL <- ggplot(vb, aes(date, AlphaRatio)) + box +
    geom_area(fill = HEAT_RED, alpha = 0.16) + geom_line(colour = HEAT_RED, linewidth = 1.1) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(x = NULL, y = "Alpha (SGTF) share of deaths") + xsc + pub_theme(17)
  cov_fit <- function(dose, legend) {
    obs <- VAC$raw %>% filter(Dose == dose) %>%
      mutate(date = wser_to_date(wser(epiYear, epiWeek)), AgeGroup = factor(AgeGroup, levels = CFG$age_groups))
    mod <- predict_coverage() %>% filter(Dose == dose) %>% select(epiYear, epiWeek, all_of(CFG$age_groups)) %>%
      pivot_longer(-c(epiYear, epiWeek), names_to = "AgeGroup", values_to = "cov") %>%
      mutate(date = wser_to_date(wser(epiYear, epiWeek)), AgeGroup = factor(AgeGroup, levels = CFG$age_groups))
    p <- ggplot() +
      geom_line(data = mod, aes(date, cov, colour = AgeGroup), linewidth = 0.9) +
      geom_point(data = obs, aes(date, Fraction, colour = AgeGroup), size = 1.2, shape = 21, fill = "white", stroke = 0.6) +
      scale_y_continuous(labels = percent, limits = c(0, 1)) + age_scale() +
      labs(x = NULL, y = paste0("Dose-", dose, " coverage"), colour = "Age") + xsc + pub_theme(17)
    if (legend) p + .inleg(0.03, 0.97, "vertical") + guides(colour = guide_legend(nrow = 3, byrow = TRUE))
    else p + theme(legend.position = "none")
  }
  (pTL | cov_fit("1", FALSE)) / (pBL | cov_fit("2", TRUE))
}

## ===== Figure 6: sex-specific campaign by dose + age-difference translation =
fig6 <- function(RES, s_rep = 2) {
  ag  <- CFG$age_groups
  xr  <- as.Date(c("2020-12-01", "2021-05-20")); xsc <- scale_x_date(limits = xr, date_breaks = "2 months", date_labels = "%b %y")
  cov_age <- predict_coverage(s_param = s_rep, method = "parameter", sex = TRUE) %>% filter(Dose %in% c("1","2")) %>%
    select(epiYear, epiWeek, Dose, all_of(c(paste0("m.",ag), paste0("f.",ag)))) %>%
    pivot_longer(-c(epiYear, epiWeek, Dose), names_to = "key", values_to = "cov") %>%
    mutate(Sex = substr(key,1,1), AgeGroup = factor(substr(key,3,nchar(key)), levels = ag),
           date = wser_to_date(wser(epiYear, epiWeek)), Dose = paste0("Dose ", Dose))
  cov_panel <- function(dz, legend) {
    p <- ggplot(filter(cov_age, Dose == dz), aes(date, cov, colour = AgeGroup, linetype = Sex)) +
      geom_line(linewidth = 0.85) + scale_y_continuous(labels = percent, limits = c(0, 1)) + age_scale() +
      scale_linetype_manual(values = c(m = "solid", f = "dashed"), labels = SEX_LAB, name = NULL) +
      labs(x = NULL, y = paste0(dz, " coverage (s = ", s_rep, ")"), colour = "Age") + xsc + pub_theme(16)
    if (!legend) p + .inleg(0.13, 0.97, "vertical") +
      guides(colour = guide_legend(nrow = 3, byrow = TRUE, order = 1), linetype = guide_legend(order = 2))
    else p + theme(legend.position = "none")
  }
  prof <- map_dfr(c(1, 1.5, 2, 2.5), function(s) mean_vacc_age(s) %>% mutate(s = factor(s))) %>% filter(!is.nan(age_diff))
  ad_panel <- function(dz, legend) {
    p <- ggplot(filter(prof, Dose == dz), aes(date, age_diff, colour = s)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") + geom_line(linewidth = 1.1) +
      scale_colour_viridis_d(option = "D", end = 0.9, name = "Prioritization s") +
      labs(x = NULL, y = paste0(dz, " age difference (yrs)")) + xsc + pub_theme(16)
    if (legend) p + .inleg(0.03, 0.93, "horizontal") else p + theme(legend.position = "none")
  }
  (cov_panel("Dose 1", TRUE) | cov_panel("Dose 2", FALSE)) /
    (ad_panel("Dose 1", TRUE)  | ad_panel("Dose 2", FALSE))
}

## ===== Figure 7: equity strategies, age-difference translation, trade-off ===
fig7 <- function(RES) {
  st <- RES$strategies %>% filter(!is.na(s_star)); smax <- max(RES$ratios$s)
  rl <- RES$ratios %>% select(s, `Equal Mortality` = EqualMortality_ratio, `Equal YLL` = EqualYLL_ratio,
                              `Equal Death Reduction` = EqualDeathReduction_ratio) %>%
    pivot_longer(-s, names_to = "Strategy", values_to = "ratio")
  pA <- ggplot(rl, aes(s, ratio, colour = Strategy)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(data = st, aes(xintercept = s_star, colour = strategy), linetype = "dotted", linewidth = 0.8, show.legend = FALSE) +
    geom_line(linewidth = 1.3) + scale_colour_manual(values = STRAT_COL) +
    lims(x=c(0.4,3.1))+
    labs(x = NULL, y = "Male / female outcome ratio", colour = NULL) + pub_theme(17) + .inleg(0.97, 0.97)
  ad_s <- AGE_DIFF_GRID %>% filter(s <= smax) %>%
    pivot_longer(c(ad1, ad2), names_to = "Dose", values_to = "ad") %>%
    mutate(Dose = ifelse(Dose == "ad1", "Dose 1", "Dose 2"))
  pMid <- ggplot(ad_s, aes(s, ad, linetype = Dose)) +
    geom_vline(data = st, aes(xintercept = s_star, colour = strategy), linetype = "dotted", linewidth = 0.8) +
    geom_line(linewidth = 1.2, colour = STRAT_ACCENT) +
    scale_linetype_manual(values = c(`Dose 1` = "solid", `Dose 2` = "dashed"), name = NULL) +
    scale_colour_manual(values = STRAT_COL, guide = "none") +
    lims(x=c(0.4,3.1))+
    labs(x = NULL, y = "Vaccination age difference,\nmale - female (yrs)") + pub_theme(17) + .inleg(0.97, 0.97)
  dl <- RES$ratios %>% select(s, Overall = dDeath_overall, Male = dDeath_m, Female = dDeath_f) %>%
    pivot_longer(-s, names_to = "Group", values_to = "d") %>% mutate(d = 100 * d)
  pB <- ggplot(dl, aes(s, d, colour = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(data = st, aes(xintercept = s_star), linetype = "dotted", colour = "grey55", linewidth = 0.7) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = c(Overall = NEUTRAL_GREY, Male = IFR_SEX_DARK[["m"]], Female = IFR_SEX_DARK[["f"]])) +
    lims(x=c(0.4,3.1))+
    labs(x = "Male prioritization weight s", y = "Death change vs baseline (%)", colour = NULL) + pub_theme(17) +
    theme(legend.position = c(0.3, 0.95), legend.justification = c(0.5, 1), legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", colour = NA))
  pA / pMid / pB
}

## ===== Figure S8: robustness in s* and in age-difference units ===============
## S8 Monte-Carlo draws (shared by Figure S8 and the results workbook so the
## reported numbers match the figure exactly). Resamples OR/IFR parameters and
## solves for the parity prioritization s* per strategy on each draw.
compute_s8_draws <- function(EXP, cc, n = 100) {
  map_dfr(seq_len(n), function(i) {
    P <- sample_params(); rat <- sweep_cached(build_baseline_exposure(P), cc, P)
    as_tibble(as.list(solve_from(rat))) %>% mutate(draw = i)
  }) %>% pivot_longer(-draw, names_to = "Strategy", values_to = "s_star") %>% filter(!is.na(s_star)) %>%
    mutate(age_diff = age_diff_lookup(s_star, "1"))
}
figS8 <- function(EXP, cc, n = 100, draws = NULL) {
  if (is.null(draws)) draws <- compute_s8_draws(EXP, cc, n)
  draws <- draws %>% mutate(Strategy = factor(Strategy, levels = c("Equal Death Reduction","Equal Mortality","Equal YLL")))
  vio <- function(xvar, xlab, ydrop) {
    p <- ggplot(draws, aes({{ xvar }}, Strategy, fill = Strategy)) +
      geom_violin(alpha = 0.5, colour = NA) + geom_boxplot(width = 0.16, outlier.size = 0.5) +
      scale_fill_manual(values = STRAT_COL, guide = "none") + labs(x = xlab, y = NULL) + pub_theme(17)
    if (ydrop) p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) else p
  }
  vio(s_star, "Parity prioritization s*", FALSE) | vio(age_diff, "Vaccination age difference (years)", TRUE)
}

## ===== Figure S9: attack-rate robustness (exact s*) + age-difference =========
figS9 <- function(EXP, cc) {
  scen <- c("quasiEmpirical","equalSex","equalAgeSex")
  res <- map_dfr(scen, function(sc)
    as_tibble(as.list(solve_from(sweep_cached(redistribute_exposure(EXP, sc), cc)))) %>% mutate(scenario = sc)) %>%
    pivot_longer(-scenario, names_to = "Strategy", values_to = "s_star") %>%
    mutate(scenario = factor(scenario, levels = scen), age_diff = age_diff_lookup(s_star, "1"))
  dg <- position_dodge(0.8)
  pL <- ggplot(res, aes(scenario, s_star, fill = Strategy)) +
    geom_col(position = dg, width = 0.72) +
    geom_text(aes(label = ifelse(is.na(s_star), "n/a", sprintf("%.1f", s_star))),
              position = dg, vjust = -0.3, size = 4, colour = "grey25") +
    scale_fill_manual(values = STRAT_COL) + scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Parity prioritization s*", fill = NULL) + pub_theme(16) +
    theme(axis.text.x = element_text(angle = 12, hjust = 0.5,vjust = 0.75)) + .inleg(0.03, 0.97)
  pR <- ggplot(res, aes(scenario, age_diff, fill = Strategy)) +
    geom_col(position = dg, width = 0.72) +
    geom_text(aes(label = ifelse(is.na(age_diff), "n/a", sprintf("%.1f", age_diff))),
              position = dg, vjust = 1.2, size = 4, colour = "grey25") +
    scale_fill_manual(values = STRAT_COL, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.08, 0))) +
    labs(x = NULL, y = "Vaccination age difference,\nmale - female (yrs)") + pub_theme(16) +
    theme(axis.text.x = element_text(angle = 12, hjust = 0.5,vjust = 0.75))
  pL | pR
}

## ===== Figure S10: wave-vaccination timing sensitivity (3 panels) ============
figS10 <- function(EXP, cc) {
  bc0  <- coverage_shares(predict_coverage(), sex = FALSE)
  cc15 <- coverage_shares(predict_coverage(s_param = 1.5, method = "parameter", sex = TRUE), sex = TRUE)
  Dbase <- project_deaths(EXP, bc0); gm <- function(D, s) sum(D$Deaths[D$Sex == s])
  ref_m <- gm(Dbase,"m"); ref_f <- gm(Dbase,"f"); ref_t <- ref_m + ref_f
  dchg <- map_dfr(seq(-8, 8, 1), function(o) {
    Dn <- project_deaths(EXP, pad_and_shift(bc0, o)); Dc <- project_deaths(EXP, pad_and_shift(cc15, o))
    tibble(offset = o, `Total (neutral)` = (gm(Dn,"m") + gm(Dn,"f")) / ref_t - 1,
           `Male (neutral)` = gm(Dn,"m")/ref_m - 1, `Female (neutral)` = gm(Dn,"f")/ref_f - 1,
           `Male (s=1.5)` = gm(Dc,"m")/ref_m - 1, `Female (s=1.5)` = gm(Dc,"f")/ref_f - 1)
  }) %>% pivot_longer(-offset, names_to = "Series", values_to = "d") %>%
    mutate(d = 100 * d, Sex = ifelse(grepl("Male", Series), "m", ifelse(grepl("Female", Series), "f", "t")),
           strat = ifelse(grepl("neutral", Series), "sex-neutral", "sex-specific (s=1.5)"))
  pA <- ggplot(dchg, aes(offset, d, colour = Sex, linetype = strat)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") + geom_line(linewidth = 1.15) +
    scale_colour_manual(values = c(m = IFR_SEX_DARK[["m"]], f = IFR_SEX_DARK[["f"]], t = NEUTRAL_GREY),
                        labels = c(f = "Female", m = "Male", t = "Total")) +
    scale_linetype_manual(values = c(`sex-neutral` = "dashed", `sex-specific (s=1.5)` = "solid"), name = NULL) +
    labs(x = NULL, y = "Death change vs observed (%)", colour = NULL) + pub_theme(16) +
    theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1), legend.box = "horizontal",
          legend.background = element_rect(fill = "white", colour = NA))
  ss <- map_dfr(seq(-8, 8, 1), function(o)
    as_tibble(as.list(solve_from(sweep_cached(EXP, cc, offset = o)))) %>% mutate(offset = o)) %>%
    pivot_longer(-offset, names_to = "Strategy", values_to = "s_star")
  pB <- ggplot(ss, aes(offset, s_star, colour = Strategy)) + geom_line(linewidth = 1.2) + geom_point(size = 1.6) +
    scale_colour_manual(values = STRAT_COL) +
    labs(x = NULL, y = "Parity prioritization s*", colour = NULL) + pub_theme(16) + .inleg(0.03, 0.97)
  pC <- ggplot(ss %>% mutate(ad = age_diff_lookup(s_star, "1")), aes(offset, ad, colour = Strategy)) +
    geom_line(linewidth = 1.2) + geom_point(size = 1.6) + scale_colour_manual(values = STRAT_COL, guide = "none") +
    labs(x = "Vaccination timing offset (weeks; + = later)",
         y = "Age difference at s*,\nmale - female (yrs)", colour = NULL) + pub_theme(16)
  pA / pB / pC + plot_layout(heights = c(1.05, 1, 1))
}
run_analysis <- function(s_grid = CFG$s_grid) {
  message("Building baseline exposure (", CFG$incidence_scenario, ") ...")
  EXP   <- build_baseline_exposure()
  D0    <- project_deaths(EXP, coverage_shares0(null_allocation(EXP)))
  Dbase <- project_deaths(EXP, coverage_shares(predict_coverage(), sex = FALSE))
  message("Sweeping prioritization s over ", length(s_grid), " values ...")
  sw  <- strategy_sweep(EXP, s_grid = s_grid)
  st  <- solve_strategies(sw$ratios)
  ets <- outcomes(project_deaths(EXP, coverage_shares(predict_coverage(method = "IFR", sex = TRUE), sex = TRUE)), D0) %>%
    mutate(strategy = "Equal IFR Threshold")
  list(exposure = EXP, null = D0, baseline = Dbase, outcomes = sw$outcomes,
       ratios = sw$ratios, strategies = st, ets = ets)
}
summarise_run <- function(RES) {
  base <- outcomes(RES$baseline, RES$null)
  cat("\n================ UK counterfactual: results ================\n")
  cat("Null (no vaccination) deaths :", round(sum(RES$null$Deaths)),
      "   sex-neutral baseline :", round(sum(RES$baseline$Deaths)), "\n\n")
  print(base %>% mutate(across(c(M, Y, P), ~round(., 4)), D = round(D)))
  cat("\nParity prioritization s* per criterion:\n"); print(RES$strategies %>% mutate(s_star = round(s_star, 3)))
  cat("\nEqual IFR Threshold (IFR-weighted) allocation:\n")
  print(RES$ets %>% mutate(across(c(M, Y, P), ~round(., 4)), D = round(D)))
  invisible(base)
}

## ===== Results workbook: key results, key parameters, and figure data =======
## Mirrors the analytical script's summary export. Reuses the same engines as
## the figures (solve_from / sweep_cached / age_diff_lookup) so the exported
## numbers match the figures; the random S8 draws are passed in (shared with
## Figure S8); the deterministic S9/S10 sweeps are recomputed here.
export_uk_results <- function(RES, EXP, s8_draws, cc_s9, cc_t10,
                              file = file.path(CFG$out_dir, "UK_counterfactual_results.xlsx"),
                              offsets = seq(-8, 8, 1)) {
  message("Building UK results workbook ...")
  base <- outcomes(RES$baseline, RES$null)
  
  ## key model parameters
  key_parameters <- tibble(
    parameter = c("population (demography)", "population year", "infection->death lag (weeks)",
                  "Alpha mortality RR vs non-VOC", "death rate-ratio (1 dose)", "death rate-ratio (2 doses)",
                  "infection rate-ratio (1 dose)", "infection rate-ratio (2 doses)",
                  "main sweep s grid", "attack-rate scenario (main)", "vaccination optimiser"),
    value = c(CFG$pop_code, as.character(CFG$pop_year), as.character(CFG$lag_weeks),
              as.character(unname(CFG$variant_mort["Alpha"])),
              as.character(round(unname(R_death["1"]), 4)), as.character(round(unname(R_death["2"]), 4)),
              as.character(round(unname(R_inf["1"]), 4)),  as.character(round(unname(R_inf["2"]), 4)),
              sprintf("%.1f-%.1f step %.1f", min(CFG$s_grid), max(CFG$s_grid), CFG$s_grid[2] - CFG$s_grid[1]),
              CFG$incidence_scenario, CFG$vacc_algorithm))
  
  ## IFR by death band (per sex)
  IFR_by_band <- IFR_GROUP %>% tidyr::pivot_wider(names_from = Sex, values_from = IFR) %>%
    transmute(death_band = as.character(AgeGroup), IFR_male = m, IFR_female = f,
              male_female_ratio = m / f)
  
  ## death totals
  null_t <- sum(RES$null$Deaths); base_t <- sum(RES$baseline$Deaths)
  death_totals <- tibble(quantity = c("null (no vaccination) deaths", "sex-neutral baseline deaths",
                                      "deaths prevented by vaccination", "fraction of deaths prevented"),
                         value = c(round(null_t), round(base_t), round(null_t - base_t), round(1 - base_t / null_t, 4)))
  
  ## parity prioritization s* per strategy (+ Equal IFR Threshold) with the
  ## campaign-mean vaccination-age difference (years) at s* and the death change
  ## relative to the sex-neutral baseline at s*, in % (Figure 7)
  dchg_at <- function(ss, col) if (is.na(ss)) NA_real_ else round(100 * approx(RES$ratios$s, RES$ratios[[col]], ss, rule = 2)$y, 2)
  parity_summary <- bind_rows(RES$strategies, tibble(strategy = "Equal IFR Threshold", s_star = NA_real_)) %>%
    mutate(age_diff_dose1 = age_diff_lookup(s_star, "1"),
           age_diff_dose2 = age_diff_lookup(s_star, "2"),
           death_change_overall_pct = vapply(s_star, function(s) dchg_at(s, "dDeath_overall"), numeric(1)),
           death_change_male_pct    = vapply(s_star, function(s) dchg_at(s, "dDeath_m"), numeric(1)),
           death_change_female_pct  = vapply(s_star, function(s) dchg_at(s, "dDeath_f"), numeric(1)))
  
  ## S8: OR/IFR Monte-Carlo robustness summary (from the shared draws; Figure S8)
  q <- function(x, p) unname(quantile(x, p, names = FALSE))
  or_ifr_robustness <- s8_draws %>% group_by(strategy = as.character(Strategy)) %>%
    summarise(n_draws = n(),
              s_star_median = round(median(s_star), 3), s_star_q025 = round(q(s_star, .025), 3),
              s_star_q975 = round(q(s_star, .975), 3), age_diff_median = round(median(age_diff), 2),
              age_diff_q025 = round(q(age_diff, .025), 2), age_diff_q975 = round(q(age_diff, .975), 2),
              .groups = "drop")
  
  ## S9: attack-rate (infection redistribution) scenarios (Figure S9)
  attackrate_sensitivity <- map_dfr(c("quasiEmpirical", "equalSex", "equalAgeSex"), function(sc)
    as_tibble(as.list(solve_from(sweep_cached(redistribute_exposure(EXP, sc), cc_s9)))) %>% mutate(scenario = sc)) %>%
    pivot_longer(-scenario, names_to = "strategy", values_to = "s_star") %>%
    mutate(age_diff = age_diff_lookup(s_star, "1"))
  
  ## S10: wave-vaccination timing sensitivity (Figure S10)
  timing_sensitivity <- map_dfr(offsets, function(o)
    as_tibble(as.list(solve_from(sweep_cached(EXP, cc_t10, offset = o)))) %>% mutate(offset_weeks = o)) %>%
    pivot_longer(-offset_weeks, names_to = "strategy", values_to = "s_star") %>%
    mutate(age_diff = age_diff_lookup(s_star, "1"))
  
  ## Epidemic-wave summary (Figure 5, left). Two windows are reported: the full
  ## death series displayed in Figure 5 (the wave and the full non-VOC->Alpha
  ## transition), and the narrower analysis (reconstruction) window, which spans
  ## only the weeks where lagged coverage and variant data overlap the deaths.
  ## Observed deaths in the analysis window are matched on death week x age band
  ## x sex, giving a like-for-like comparison with the baseline reconstruction.
  .ld   <- CFG$lag_weeks
  .dcf  <- deathCounts %>% mutate(wk_death = wser(epiYear, epiWeek))           # full death series
  .wdf  <- sort(unique(.dcf$wk_death))
  .sgf  <- SGTFdeaths %>% mutate(wk_death = wser(epiYear, epiWeek)) %>% filter(wk_death %in% .wdf) %>% arrange(wk_death)
  .keys <- EXP %>% distinct(wk_inf, AgeGroup, Sex) %>% mutate(wk_death = wk_inf + .ld)
  .dcm  <- .dcf %>% semi_join(.keys, by = c("wk_death", "AgeGroup", "Sex"))    # observed, analysis window
  .wd   <- sort(unique(EXP$wk_inf + .ld))
  .sgm  <- SGTFdeaths %>% mutate(wk_death = wser(epiYear, epiWeek)) %>% filter(wk_death %in% .wd) %>% arrange(wk_death)
  .obs_a <- sum(.dcm$Value)
  epidemic_wave <- tibble(
    quantity = c(
      "display period start (full death series)", "display period end (full death series)",
      "observed deaths, full window, total", "observed deaths, full window, male", "observed deaths, full window, female",
      "Alpha share of deaths, first displayed week (%)", "Alpha share of deaths, last displayed week (%)",
      "analysis (reconstruction) period start", "analysis (reconstruction) period end",
      "observed deaths in analysis window, total", "observed deaths in analysis window, male", "observed deaths in analysis window, female",
      "Alpha share of deaths, first analysis week (%)", "Alpha share of deaths, last analysis week (%)",
      "baseline reconstruction vs observed (analysis window), deviation (%)"),
    value = c(
      as.character(wser_to_date(min(.wdf))), as.character(wser_to_date(max(.wdf))),
      as.character(round(sum(.dcf$Value))), as.character(round(sum(.dcf$Value[.dcf$Sex == "m"]))), as.character(round(sum(.dcf$Value[.dcf$Sex == "f"]))),
      as.character(round(100 * .sgf$AlphaRatio[1], 1)), as.character(round(100 * .sgf$AlphaRatio[nrow(.sgf)], 1)),
      as.character(wser_to_date(min(.wd))), as.character(wser_to_date(max(.wd))),
      as.character(round(.obs_a)), as.character(round(sum(.dcm$Value[.dcm$Sex == "m"]))), as.character(round(sum(.dcm$Value[.dcm$Sex == "f"]))),
      as.character(round(100 * .sgm$AlphaRatio[1], 1)), as.character(round(100 * .sgm$AlphaRatio[nrow(.sgm)], 1)),
      as.character(round(100 * abs(base_t - .obs_a) / max(.obs_a, 1), 2))))
  
  ## Baseline (s=1) and Equal IFR Threshold male-to-female outcome ratios (Figure 7A)
  .mk_ratio <- function(df, label) tibble(allocation = label,
                                          Mm_Mf = df$M[df$Sex == "m"] / df$M[df$Sex == "f"],
                                          Ym_Yf = df$Y[df$Sex == "m"] / df$Y[df$Sex == "f"],
                                          Pm_Pf = df$P[df$Sex == "m"] / df$P[df$Sex == "f"])
  outcome_ratios <- bind_rows(.mk_ratio(base, "sex-neutral baseline (s=1)"),
                              .mk_ratio(RES$ets, "Equal IFR Threshold")) %>%
    mutate(across(c(Mm_Mf, Ym_Yf, Pm_Pf), ~ round(.x, 4)))
  
  ## Death change vs the sex-neutral baseline across vaccination-timing offsets
  ## (Figure S10, top): sex-neutral and s = 1.5 campaigns, overall and by sex.
  .bc0  <- coverage_shares(predict_coverage(), sex = FALSE)
  .cc15 <- coverage_shares(predict_coverage(s_param = 1.5, method = "parameter", sex = TRUE), sex = TRUE)
  .Dref <- project_deaths(EXP, .bc0); .gm <- function(D, s) sum(D$Deaths[D$Sex == s])
  .rm <- .gm(.Dref, "m"); .rf <- .gm(.Dref, "f"); .rt <- .rm + .rf
  timing_death_change <- map_dfr(offsets, function(o) {
    Dn <- project_deaths(EXP, pad_and_shift(.bc0, o)); Dc <- project_deaths(EXP, pad_and_shift(.cc15, o))
    tibble(offset_weeks = o,
           neutral_overall_pct = round(100 * ((.gm(Dn,"m") + .gm(Dn,"f")) / .rt - 1), 2),
           neutral_male_pct    = round(100 * (.gm(Dn,"m") / .rm - 1), 2),
           neutral_female_pct  = round(100 * (.gm(Dn,"f") / .rf - 1), 2),
           s1p5_male_pct       = round(100 * (.gm(Dc,"m") / .rm - 1), 2),
           s1p5_female_pct     = round(100 * (.gm(Dc,"f") / .rf - 1), 2)) })
  
  README <- tibble(
    sheet = c("epidemic_wave", "key_parameters", "IFR_by_band", "death_totals", "outcomes_by_s",
              "ratios_by_s", "outcome_ratios", "parity_summary", "baseline_by_sex", "ets_by_sex",
              "age_difference_by_s", "or_ifr_robustness", "attackrate_sensitivity",
              "timing_sensitivity", "timing_death_change", "deaths_null", "deaths_baseline"),
    contents = c("Analysis window, observed deaths by sex, Alpha share of deaths (first/last week), and baseline-vs-observed deviation (Figure 5, left).",
                 "Key model parameters (rate-ratios, variant RR, lag, optimiser).",
                 "Infection fatality ratio by coarse death band, per sex.",
                 "Null vs sex-neutral-baseline death totals and the fraction prevented.",
                 "Outcome metrics (mortality M, YLL Y, prevented fraction P, deaths D) by sex over the prioritization sweep s (Figure 7).",
                 "Male/female outcome ratios over s, with the per-s death change vs baseline (dDeath, fraction); the crossing of 1 gives each strategy's parity s* (Figure 7A, 7C).",
                 "Male-to-female mortality / YLL / prevented-fraction ratios for the sex-neutral baseline (s=1) and the Equal IFR Threshold allocation (Figure 7A).",
                 "Parity prioritization s* per strategy (+ Equal IFR Threshold) with the campaign-mean vaccination-age difference and the death change vs baseline (overall, male, female, %) at s* (Figure 7).",
                 "Sex-neutral baseline outcomes vs the no-vaccination counterfactual.",
                 "Equal IFR Threshold (IFR-weighted) allocation outcomes by sex.",
                 "Campaign-mean vaccination-age difference (male-female, years) vs s for dose 1 (ad1) and dose 2 (ad2) (Figures 6, 7B).",
                 "OR/IFR Monte-Carlo robustness: median and 95% interval of parity s* per strategy, in s and in vaccination-age-difference units (Figure S8).",
                 "Attack-rate robustness: parity s* and age difference per strategy under three infection-redistribution scenarios (Figure S9).",
                 "Timing robustness: parity s* and age difference per strategy across vaccination-timing offsets (Figure S10, middle/bottom).",
                 "Death change vs the sex-neutral baseline across vaccination-timing offsets, for the sex-neutral and s=1.5 campaigns, overall and by sex (Figure S10, top).",
                 "Reconstructed deaths by age/sex under the no-vaccination counterfactual.",
                 "Reconstructed deaths by age/sex under the sex-neutral baseline campaign."))
  
  writexl::write_xlsx(list(
    README = README, epidemic_wave = epidemic_wave, key_parameters = key_parameters,
    IFR_by_band = IFR_by_band, death_totals = death_totals,
    outcomes_by_s = RES$outcomes, ratios_by_s = RES$ratios, outcome_ratios = outcome_ratios,
    parity_summary = parity_summary, baseline_by_sex = base, ets_by_sex = RES$ets,
    age_difference_by_s = AGE_DIFF_GRID, or_ifr_robustness = or_ifr_robustness,
    attackrate_sensitivity = attackrate_sensitivity, timing_sensitivity = timing_sensitivity,
    timing_death_change = timing_death_change, deaths_null = RES$null, deaths_baseline = RES$baseline), file)
  message("  wrote ", file, " (", 18, " sheets)")
  invisible(file)
}

if (!exists("SKIP_MAIN")) {
  RES <- run_analysis(); summarise_run(RES)
  
  message("Precomputing s -> vaccination-age-difference grid (both doses) ...")
  AGE_DIFF_GRID <<- map_dfr(seq(0.5, 15, 0.1), function(s) {
    p <- mean_vacc_age(s) %>% filter(wk %in% CAMPAIGN_WKS, !is.nan(age_diff)) %>%
      group_by(Dose) %>% summarise(ad = mean(age_diff), .groups = "drop")
    tibble(s = s, ad1 = p$ad[p$Dose == "Dose 1"], ad2 = p$ad[p$Dose == "Dose 2"])
  })
  cc_s8   <- cache_coverage(seq(0.5, 4, 0.1))         # OR/IFR robustness: parity s* is modest
  cc_t10  <- cache_coverage(seq(0.5, 15, 0.1))       # timing: s* search constrained to 15
  cc_s9   <- cache_coverage(seq(0.5, 15, 0.1))        # attack-rate: unconstrained search for exact s*
  s8_draws <- compute_s8_draws(RES$exposure, cc_s8, 100)  # shared by the workbook and Figure S8
  
  ## results workbook (key results + key parameters + figure data)
  export_uk_results(RES, RES$exposure, s8_draws = s8_draws, cc_s9 = cc_s9, cc_t10 = cc_t10)
  
  message("Rendering figures (SVG + PNG) ...")
  figdir <- file.path(CFG$out_dir, "figures"); dir.create(figdir, showWarnings = FALSE)
  save_fig <- function(p, name, w, h) tryCatch({
    ggsave(file.path(figdir, paste0(name, ".svg")), p, width = w, height = h)
    ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, dpi = 140)
    cat("  wrote", name, "\n")
  }, error = function(e) cat("  FAILED", name, ":", conditionMessage(e), "\n"))
  
  save_fig(fig5(),                              "Figure5_wave_and_fit",    14, 9.5)
  save_fig(fig6(RES),                           "Figure6_vaccination_sex", 14, 11)
  save_fig(fig7(RES),                           "Figure7_strategies",      10, 13)
  save_fig(figS8(RES$exposure, cc_s8, draws = s8_draws), "FigureS8_robustness", 12, 6.5)
  save_fig(figS9(RES$exposure, cc_s9),          "FigureS9_attack_rate",    12, 7)
  save_fig(figS10(RES$exposure, cc_t10),        "FigureS10_timing",        10, 14)
  cat("\nDone. Outputs in", CFG$out_dir, "\n")
}