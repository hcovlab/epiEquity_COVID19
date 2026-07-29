# Libraries----
library(optimx)
library(tidyverse)
library(scales)
library(gridExtra)
library(ggrepel)
library(png)
library(grid)
library(ggpattern)
library(patchwork)   # composite main figure (shift profile + heatmap) and the sensitivity grids

#load("Results0501.all.4strat_exp2.RData", envir = .GlobalEnv)

# Lists----
functions_optim <- list(
  genlog = list(
    nop = 4,
    body = "parmsvector['par1']/(1+exp(-parmsvector['par2']*(age-(parmsvector['par3']))))+parmsvector['par4']",
    body.m = "pars.m[1]/(1+exp(-pars.m[2]*(age.d-(pars.m[3]))))+pars.m[4]",
    body.f = "pars.f[1]/(1+exp(-pars.f[2]*(age.d-(pars.f[3]))))+pars.f[4]",
    parints = list(c(0, 1), c(0.01, 1), c(10, 200), c(1e-6, 1e-4))
  ),
  exp2 = list(
    nop = 2,
    body = "parmsvector['par1']*exp(age*parmsvector['par2'])",
    body.m = "pars.m[1]*exp(age.d*pars.m[2])",
    body.f = "pars.f[1]*exp(age.d*pars.f[2])",
    parints = list(c(1e-7, 1e-5), c(0.01, 0.5))
  ),
  exp3 = list(
    nop = 3,
    body = "parmsvector['par1']*exp(parmsvector['par2']*age)+parmsvector['par3']",
    body.m = "pars.m[1]*exp(pars.m[2]*age.d)+pars.m[3]",
    body.f = "pars.f[1]*exp(pars.f[2]*age.d)+pars.f[3]",
    parints = list(c(1e-7, 1e-5), c(0.01, 0.05), c(1e-8, 1e-4))
  ),
  G = list(
    nop = 3,
    body = "parmsvector['par1']*exp(-exp(parmsvector['par2']-parmsvector['par3']*age))",
    body.m = "pars.m[1]*exp(-exp(pars.m[2]-pars.m[3]*age.d))",
    body.f = "pars.f[1]*exp(-exp(pars.f[2]-pars.f[3]*age.d))",
    parints = list(c(0.01, 1), c(1,5), c(0.001,0.1))
  )
)

functions_optim_shift <- list(
  genlog = list(
    nop = 1,
    body = "pars.f[1]/(1+exp(-pars.f[2]*(age.d-(pars.f[3])+parmsvector['par1'])))+pars.f[4]",
    parints = list(c(1, 5))
  ),
  exp2 = list(
    nop = 1,
    body = "pars.f[1]*exp((age.d+parmsvector['par1'])*pars.f[2])",
    parints = list(c(1, 5))
  ),
  exp3 = list(
    nop = 1,
    body = "pars.f[1]*exp((age.d+parmsvector['par1'])*pars.f[2])+pars.f[3]",
    parints = list(c(1, 5))
  ),
  G = list(
    nop = 1,
    body = "pars.f[1]*exp(-exp(pars.f[2]-pars.f[3]*(age.d+parmsvector['par1'])))",
    parints = list(c(1, 5))
  )
)

IFR_optim <- list(
  ODriscoll = list(
    male = tibble(
      lower = c(0.002, 0.0005, 0.001, 0.002, 0.007, 0.014, 0.028, 0.048, 0.091, 0.144, 0.250, 0.384, 0.511, 1.244, 1.976, 3.651, 9.276),
      median = c(0.003, 0.001, 0.001, 0.003, 0.008, 0.017, 0.033, 0.056, 0.106, 0.168, 0.291, 0.448, 0.595, 1.452, 2.307, 4.260, 10.825),
      upper = c(0.004, 0.001, 0.002, 0.003, 0.009, 0.020, 0.038, 0.065, 0.123, 0.195, 0.336, 0.518, 0.688, 1.680, 2.668, 4.929, 12.523)
    ) / 100,
    female = tibble(
      lower = c(0.002, 0.0005, 0.0005, 0.002, 0.004, 0.008, 0.013, 0.021, 0.038, 0.063, 0.105, 0.169, 0.273, 0.598, 0.893, 1.838, 4.935),
      median = c(0.003, 0.001, 0.001, 0.002, 0.005, 0.009, 0.015, 0.025, 0.044, 0.073, 0.123, 0.197, 0.318, 0.698, 1.042, 2.145, 5.759),
      upper = c(0.003, 0.001, 0.001, 0.003, 0.006, 0.011, 0.018, 0.029, 0.051, 0.085, 0.142, 0.228, 0.367, 0.807, 1.206, 2.482, 6.662)
    ) / 100,
    ints = c(seq(0, 80, by = 5), 100),
    gender = TRUE
  ),
  Salje = list(
    male = tibble(
      lower = c(0.0005, 0.003, 0.02, 0.03, 0.1, 0.6, 1.5, 7.9),
      median = c(0.001, 0.007, 0.03, 0.06, 0.2, 1, 2.7, 14),
      upper = c(0.003, 0.01, 0.05, 0.1, 0.4, 1.6, 3.4, 22.7)
    ) / 100,
    female = tibble(
      lower = c(0.0005, 0.002, 0.007, 0.02, 0.08, 0.3, 0.7, 2.7),
      median = c(0.001, 0.004, 0.01, 0.03, 0.1, 0.5, 1.3, 4.9),
      upper = c(0.002, 0.007, 0.02, 0.06, 0.2, 0.8, 2.1, 8)
    ) / 100,
    ints = c(0, 20, 30, 40, 50, 60, 70, 80, 100),
    gender = TRUE
  )
)



# Functions----
round_any <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

eval_model_expr <- function(expr, env = parent.frame()) {
  eval(parse(text = expr), envir = env)
}

intervals_from_ints <- function(ints) {
  tibble(
    age = age,
    factor = rep(seq_len(length(ints) - 1), times = diff(ints))
  )
}

weighted_interval_mean <- function(z, values, weights, intervals) {
  idx <- intervals$factor == z
  sum(values[idx] * weights[idx] / sum(weights[idx]))
}

dem.gen <- function(population,year, permillion) {
  pop <- read.table(file = "Data/Population.csv", sep = " ", dec = ".", header = T, stringsAsFactors = F)
  pop <- pop %>%
    filter(Year==year,c_code==population)%>%
    mutate(
      fr.Total = Total / sum(Total),
      fr.Male = Male / sum(Total),
      fr.Female = Female / sum(Total)
    ) %>%
    filter(Age >= age[1], Age <= last(age))
  
  if (permillion) {
    pop %>%
      mutate(
        Total = Total * 1e6 / sum(dem$Total),
        Male = Male * 1e6 / sum(dem$Total),
        Female = Female * 1e6 / sum(dem$Total)
      )
  }
  
  return(pop)
}

agebased <- function(population = population_def,
                     year = year_def,
                     optFunction,
                     optIFR = optIFR_def,
                     rF_grid = rF_grid_def,
                     fit_methods = c("nlminb", "L-BFGS-B"),    # curated, bounded, base-R optimisers
                     shift_methods = c("nlminb", "L-BFGS-B"),   # (replaces optimx all.methods)
                     n_fits = 0) {                              # extra fits to IFR draws ~ U(lower, upper) for the fit-ensemble figure
  #### Demographic data
  
  dem <- dem.gen(population, year, permillion = F)
  optim_def <- functions_optim[[optFunction]]
  shift_def <- functions_optim_shift[[optFunction]]
  ifr_def <- IFR_optim[[optIFR]]
  
  #### Agebased model setup
  
  ef <- function(par, func = optFunction, IFR = optIFR) {
    nop <- optim_def$nop
    
    for (sex in c("m", "f")) {
      start <- if (sex == "m") 1 else nop + 1
      parmsvector <- par[start:(start + nop - 1)] * parscale
      names(parmsvector) <- paste0("par", seq_len(nop))
      
      assign(paste0("fx.", sex, ".age"), eval_model_expr(optim_def$body), envir = .GlobalEnv)
      assign(paste0("pars.", sex), parmsvector, envir = .GlobalEnv)
    }
    dem_fm <- filter(dem, fr.Female >= age[1])
    assign("fx.fm.age", (fx.f.age * dem_fm$fr.Female + fx.m.age * dem_fm$fr.Male) / dem_fm$fr.Total, envir = .GlobalEnv)
    
    intervals <- intervals_from_ints(ints)
    
    interval_seq <- seq_len(length(ints) - 1)
    mort.pred.male.age <<- sapply(interval_seq, weighted_interval_mean, values = fx.m.age, weights = dem$Male, intervals = intervals)
    mort.pred.female.age <<- sapply(interval_seq, weighted_interval_mean, values = fx.f.age, weights = dem$Female, intervals = intervals)
    meanage.male <<- sapply(interval_seq, weighted_interval_mean, values = age, weights = dem$fr.Male, intervals = intervals)
    meanage.female <<- sapply(interval_seq, weighted_interval_mean, values = age, weights = dem$fr.Female, intervals = intervals)
    
    sse <- (sum((log(mort.pred.male.age) - log(mortality.male.age))^2) + sum((log(mort.pred.female.age) - log(mortality.female.age))^2)) / (length(mortality.male.age) + length(mortality.female.age))
    return(sse)
  }
  IFR.Gen <- function(coords = F) {
    if (coords) {
      for (i in 1:3) {
        ints <<- ifr_def$ints
        
        mortality.male.age <- ifr_def$male[, i]
        mortality.male.age <- mortality.male.age[ints[-length(ints)] >= age[1] & ints[-length(ints)] <= (age[length(age)])]
        mortality.female.age <- ifr_def$female[, i]
        mortality.female.age <- mortality.female.age[ints[-length(ints)] >= age[1] & ints[-length(ints)] <= (age[length(age)])]
        
        ints <<- ints[ints >= age[1] & ints <= (age[length(age)] + 1)]
        
        xcoord <<- rep(ints, rep(2, length(ints)))
        xcoord <<- xcoord[c(-1, -length(xcoord))]
        if (i == 1) {
          ycoord.male.lower <<- rep(mortality.male.age, rep(2, length(mortality.male.age)))
          ycoord.female.lower <<- rep(mortality.female.age, rep(2, length(mortality.female.age)))
        } else if (i == 2) {
          mortality.male.age <<- mortality.male.age
          mortality.female.age <<- mortality.female.age
          ycoord.male.median <<- rep(mortality.male.age, rep(2, length(mortality.male.age)))
          ycoord.female.median <<- rep(mortality.female.age, rep(2, length(mortality.female.age)))
        } else if (i == 3) {
          ycoord.male.upper <<- rep(mortality.male.age, rep(2, length(mortality.male.age)))
          ycoord.female.upper <<- rep(mortality.female.age, rep(2, length(mortality.female.age)))
        }
      }
    } else {
      ints <<- ifr_def$ints
      
      maleIFR <- ifr_def$male
      femaleIFR <- ifr_def$female
      neutralIFR <- ifr_def$mort
      
      mortality.male.age <<- runif(nrow(maleIFR), min = maleIFR$lower, max = maleIFR$upper)
      mortality.female.age <<- runif(nrow(femaleIFR), min = femaleIFR$lower, max = femaleIFR$upper)
      mortality.male.age <<- mortality.male.age[ints[-length(ints)] >= age[1] & ints[-length(ints)] <= (age[length(age)])]
      mortality.female.age <<- mortality.female.age[ints[-length(ints)] >= age[1] & ints[-length(ints)] <= (age[length(age)])]
      
      ycoord.male <<- rep(mortality.male.age, rep(2, length(mortality.male.age)))
      ycoord.female <<- rep(mortality.female.age, rep(2, length(mortality.female.age)))
      
      ints <<- ints[ints >= age[1] & ints <= (age[length(age)] + 1)]
    }
  }
  
  parints <- optim_def$parints
  parscale <- sapply(parints, function(x) runif(min = x[1], max = x[2], 1))
  firstguess <- rep(1, 2 * optim_def$nop)
  
  #### Agebased model fit
  
  IFR.Gen(coords = T)
  estimate_all <<- optimx(firstguess, ef, method = fit_methods, lower = 1e-100)
  bestmethod <- rownames(estimate_all)[which.min(estimate_all$value)]
  estimate <- optimx(firstguess, ef, method = bestmethod, lower = 0)
  
  #### Shift model setup
  
  eqRisk.model <- function(par, func = optFunction, IFR = optIFR) {
    age.d <- seq(age[1], last(age), by = 0.01)
    nop <- shift_def$nop
    parmsvector <- par
    names(parmsvector) <- "par1"
    
    fx.m.age.d <- eval_model_expr(optim_def$body.m)
    assign("fx.m.age.d.shift", eval_model_expr(shift_def$body), envir = .GlobalEnv)
    
    # sse<-(sum((log(fx.m.age.d.shift) - log(fx.m.age.d))^2) )
    sse <- (sum(((fx.m.age.d.shift) - (fx.m.age.d))^2))
    return(sse)
  }
  
  #### Shift model fit
  
  estimate_all2 <<- optimx(1, eqRisk.model, method = shift_methods, lower = 1e-100)
  bestmethod2 <- rownames(estimate_all2)[which.min(estimate_all2$value)]
  estimate2 <- optimx(1, eqRisk.model, method = bestmethod2, lower = 0)
  
  #### Setting up tibbles
  
  ycoords.age.results <<- tibble(it = numeric(), agecat = numeric(), gender = character(), value = numeric())
  sse.age.results <<- tibble(it = numeric(), sse = numeric(), shift = numeric(), sse2 = numeric())
  fx.age.age.results <<- tibble(it = numeric(), age = numeric(), gender = character(), value = numeric())
  fx.age.agegroup.results <<- tibble(it = numeric(), agecat = numeric(), gender = character(), value = numeric())
  
  #### Saving results
  
  ycoords.age.results <<- bind_rows(
    ycoords.age.results,
    tibble(it = 0, agecat = xcoord, gender = "m", value = ycoord.male.median),
    tibble(it = 0, agecat = xcoord, gender = "f", value = ycoord.female.median)
  )
  
  fx.age.agegroup.results <<- bind_rows(
    fx.age.agegroup.results,
    tibble(it = 0, agecat = meanage.male, gender = "m", value = mort.pred.male.age),
    tibble(it = 0, agecat = meanage.female, gender = "f", value = mort.pred.female.age)
  )
  
  fx.age.age.results <<- bind_rows(
    fx.age.age.results,
    tibble(it = 0, age = age[1]:(age[1]+length(fx.m.age)-1), gender = "m", value = fx.m.age),
    tibble(it = 0, age = age[1]:(age[1]+length(fx.f.age)-1), gender = "f", value = fx.f.age),
    tibble(it = 0, age = age[1]:(age[1]+length(fx.m.age)-1), gender = "ms", value = fx.m.age.d.shift[seq(age[1], last(age), by = 0.01) %% 1 == 0])
  )
  
  #### Optional fit ensemble (for the IFR fit figure): refit to IFR realisations
  #### sampled uniformly from the published 95% credible intervals, then restore
  #### the canonical (median) fit so nothing downstream is affected.
  if (n_fits > 0) {
    for (k in seq_len(n_fits)) {
      IFR.Gen(coords = FALSE)                                  # mortality.*.age ~ U(lower, upper)
      est_k <- tryCatch(optimx(firstguess, ef, method = bestmethod, lower = 0), error = function(e) NULL)
      if (is.null(est_k)) next
      invisible(ef(as.numeric(est_k[1, seq_len(2 * optim_def$nop)])))   # set fx.m.age / fx.f.age for this draw
      fx.age.age.results <<- bind_rows(
        fx.age.age.results,
        tibble(it = k, age = age[1]:(age[1]+length(fx.m.age)-1), gender = "m", value = fx.m.age),
        tibble(it = k, age = age[1]:(age[1]+length(fx.f.age)-1), gender = "f", value = fx.f.age))
    }
    IFR.Gen(coords = TRUE)                                     # restore consensus median + step coords
    invisible(ef(as.numeric(estimate[1, seq_len(2 * optim_def$nop)])))  # restore canonical fx.* / mort.pred.* / meanage.*
  }
  
  sse.age.results <<- bind_rows(
    sse.age.results,
    tibble(it = 0, sse = estimate$value, shift = estimate2$p1, sse2 = estimate2$value)
  )
  
  #### IFR into dem
  
  dem<-dem%>%
    mutate(
      IFR.Male=fx.m.age,
      IFR.Female=fx.f.age,
      IFR.s.Male=fx.m.age.d.shift[seq(age[1], last(age), by = 0.01) %% 1 == 0]
    )
  
  #### Separable threshold constants (precomputed once; constant across the grid)
  ## Every per-sex mortality metric is one-dimensional in its own threshold, so
  ## the (mMA x fMA) sweep collapses to reverse-cumulative sums. With efficiency
  ## eps = 1 - rF, deaths at ages >= the threshold are scaled by (1 - eps).
  fm_age   <- dem$Age
  fm_Dm    <- dem$Male   * dem$IFR.Male                 # deaths (attack rate a cancels in ratios)
  fm_Df    <- dem$Female * dem$IFR.Female
  fm_DmY   <- fm_Dm * dem$ex.Male                       # YLL-weighted deaths
  fm_DfY   <- fm_Df * dem$ex.Female
  fm_Nm    <- sum(dem$Male);   fm_Nf   <- sum(dem$Female)
  fm_TDm   <- sum(fm_Dm);      fm_TDf  <- sum(fm_Df)
  fm_TDmY  <- sum(fm_DmY);     fm_TDfY <- sum(fm_DfY)
  fm_AbDm  <- rev(cumsum(rev(fm_Dm)))                   # sum over ages >= each threshold
  fm_AbDf  <- rev(cumsum(rev(fm_Df)))
  fm_AbDmY <- rev(cumsum(rev(fm_DmY)))
  fm_AbDfY <- rev(cumsum(rev(fm_DfY)))
  fm_fxm   <- dem$IFR.Male                              # single-year IFR, aligned to fm_age
  fm_fxf   <- dem$IFR.Female
  fm_shift <- round(estimate2$p1, 0)
  
  #### Four strategies for the optimal female threshold (vectorised) ############
  ## For one (mMA, rF): male metrics depend only on mMA, female metrics only on
  ## fMA. Each strategy's optimum fMA is argmin |1 - mMetric/fMetric| over the
  ## candidate female thresholds fMA >= mMA. Same numbers as the original per-cell
  ## sweep (verified to machine precision), without the nested rowwise().
  fourModels <- function(mMA, rF) {
    eps <- 1 - rF
    mDR <- (fm_TDm  - eps * fm_AbDm ) / fm_Nm
    fDR <- (fm_TDf  - eps * fm_AbDf ) / fm_Nf
    mPD <-  eps * fm_AbDm / fm_TDm
    fPD <-  eps * fm_AbDf / fm_TDf
    mY  <- (fm_TDmY - eps * fm_AbDmY) / fm_Nm
    fY  <- (fm_TDfY - eps * fm_AbDfY) / fm_Nf
    
    im   <- match(mMA, fm_age)            # male threshold index
    fidx <- which(fm_age >= mMA)          # candidate female thresholds (ascending)
    fAge <- fm_age[fidx]
    
    pick <- function(mval, fvec) {        # argmin |1 - ratio|; first min -> lowest fMA on ties
      r <- mval / fvec[fidx]
      j <- which.min(abs(1 - r))
      list(fMA = fAge[j], m = mval, f = fvec[fidx][j], ratio = unname(r[j]),
           vis = isTRUE(r[j] >= 0.95 && r[j] <= 1.05))
    }
    dr <- pick(mDR[im], fDR); pd <- pick(mPD[im], fPD); yp <- pick(mY[im], fY)
    
    eqDR   <- tibble(eqDR.fMA = dr$fMA, eqDR.mDR = dr$m, eqDR.fDR = dr$f,
                     eqDR.mDR_fDR = dr$ratio, eqDR.vis = dr$vis)
    eqPD   <- tibble(eqPD.fMA = pd$fMA, eqPD.mPD = pd$m, eqPD.fPD = pd$f,
                     eqPD.mPD_fPD = pd$ratio, eqPD.vis = pd$vis)
    eqYPLL <- tibble(eqYPLL.fMA = yp$fMA, eqYPLL.mYPLL = yp$m, eqYPLL.fYPLL = yp$f,
                     eqYPLL.mYPLL_fYPLL = yp$ratio, eqYPLL.vis = yp$vis)
    
    sN <- tibble(sN.fMA = mMA, sN.fRisk = fm_fxf[im],
                 sN.fYPLL = fY[im], sN.fPD = fPD[im], sN.fDR = fDR[im])
    
    fRiskAge <- mMA + fm_shift
    if (fRiskAge <= last(age)) {
      mr <- fm_fxm[im]; fr <- fm_fxf[match(fRiskAge, fm_age)]
      eqRisk <- tibble(eqRisk.fMA = fRiskAge, eqRisk.mRisk = mr, eqRisk.fRisk = fr,
                       eqRisk.mRisk_fRisk = mr / fr,
                       eqRisk.vis = isTRUE(mr / fr >= 0.95 && mr / fr <= 1.05))
    } else {
      eqRisk <- tibble(eqRisk.fMA = NA_real_, eqRisk.mRisk = NA_real_,
                       eqRisk.fRisk = NA_real_, eqRisk.mRisk_fRisk = NA_real_,
                       eqRisk.vis = FALSE)
    }
    
    bind_cols(eqDR, eqPD, eqRisk, eqYPLL, sN)
  }
  
  # Calculating optimal female protection ages
  
  minAge.results <<- expand.grid(
    mMA = age[-1],
    rF = rF_grid
    #rF = c(0.1, 1 / 3, 0.5, 2 / 3, 0.9)
  ) %>%
    mutate(fMA.neutral = mMA,
           pop=population_def,
           year=year_def) %>%
    rowwise() %>%
    mutate(r = fourModels(mMA, rF))
  
  return(1)
}

# === Matrix (heatmap) plotting for the four strategies ========================
# The four-strategy model (fourModels) evaluates each mortality metric on the
# full (thr_m x thr_f) plane internally but only stores the equality optimum.
# The helpers below rebuild that full surface for plotting WITHOUT touching the
# model. The surface is reconstructed from objects that persist after a
# single-country agebased() run:
#   population_def, year_def : demography, via dem.gen()
#   fx.m.age, fx.f.age       : single-year IFR for males / females
#
# Notation follows the methods section: thr_f (female age threshold) and
# thr_m (male age threshold, with thr_m <= thr_f), and efficiency
# epsilon = 1 - rF of the enhanced protection applied above the threshold.
# The infection attack rate 'a' is uniform across ages and cancels in every
# ratio, so it is omitted here.

.above <- function(v) rev(cumsum(rev(v)))   # sum over all ages >= each threshold

# Per-age demography + IFR for the country currently held in the globals.
.strategy_inputs <- function() {
  dem <- dem.gen(population_def, year_def, permillion = FALSE)
  dem$IFR.Male   <- fx.m.age      # positional alignment, exactly as the model does
  dem$IFR.Female <- fx.f.age
  list(
    age   = dem$Age,
    Nm    = dem$Male,
    Nf    = dem$Female,
    Dm    = dem$Male   * dem$IFR.Male,                  # deaths (a = 1)
    Df    = dem$Female * dem$IFR.Female,
    DmYLL = dem$Male   * dem$IFR.Male   * dem$ex.Male,  # YLL-weighted deaths
    DfYLL = dem$Female * dem$IFR.Female * dem$ex.Female,
    IFRm  = dem$IFR.Male,
    IFRf  = dem$IFR.Female
  )
}

# 1-D threshold functions at one efficiency epsilon, indexed by threshold age.
# Deaths above the threshold are scaled by (1 - eps), hence e.g.
#   mDR(t) = [ total deaths - eps * (deaths at ages >= t) ] / total population.
.threshold_metrics <- function(inp, epsilon) {
  Adm <- .above(inp$Dm);    Adf <- .above(inp$Df)
  Aym <- .above(inp$DmYLL); Ayf <- .above(inp$DfYLL)
  list(
    mDR   = (sum(inp$Dm)    - epsilon * Adm) / sum(inp$Nm),
    fDR   = (sum(inp$Df)    - epsilon * Adf) / sum(inp$Nf),
    mPD   =  epsilon * Adm / sum(inp$Dm),
    fPD   =  epsilon * Adf / sum(inp$Df),
    mYPLL = (sum(inp$DmYLL) - epsilon * Aym) / sum(inp$Nm),
    fYPLL = (sum(inp$DfYLL) - epsilon * Ayf) / sum(inp$Nf),
    mRisk = inp$IFRm,
    fRisk = inp$IFRf
  )
}

# Strategy metadata: title, the male/female threshold-metric whose ratio is
# colour-mapped, a colour-bar label, whether that ratio depends on epsilon, and
# whether the colour direction is inverted. Colour is oriented so blue always
# means "males die more". For death rate / YLL / IFR that is the high-ratio end;
# for prevented deaths (mPD/fPD) a higher ratio means males are protected MORE
# (die less), so that one panel's scale is reversed (invert = TRUE).
.strategy_meta <- list(
  eqRisk = list(title = "Individual Risk Equality",   m = "mRisk", f = "fRisk", lab = "mIFR / fIFR", eps_dep = FALSE, invert = FALSE),
  eqPD   = list(title = "Benefit Equality", m = "mPD",   f = "fPD",   lab = "mPD / fPD",   eps_dep = FALSE, invert = TRUE),
  eqDR   = list(title = "Mortality Equality",       m = "mDR",   f = "fDR",   lab = "mDR / fDR",   eps_dep = TRUE,  invert = FALSE),
  eqYPLL = list(title = "Life-year Equality",             m = "mYPLL", f = "fYPLL", lab = "mYLL / fYLL", eps_dep = TRUE,  invert = FALSE)
)

# Full (thr_m x thr_f) ratio surface for one strategy at one efficiency.
.strategy_surface <- function(inp, strat, epsilon) {
  meta <- .strategy_meta[[strat]]
  mm   <- .threshold_metrics(inp, epsilon)
  M <- outer(mm[[meta$m]], mm[[meta$f]], "/")          # rows = thr_m, cols = thr_f
  dimnames(M) <- list(thr_m = inp$age, thr_f = inp$age)
  df <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  names(df)[3] <- "ratio"
  df$thr_m   <- as.integer(as.character(df$thr_m))
  df$thr_f   <- as.integer(as.character(df$thr_f))
  df$epsilon <- epsilon
  df
}

# Per-thr_f equalising optimum: for each female threshold, the male threshold
# thr_m <= thr_f whose metric ratio is closest to 1 (the same argmin the model
# uses, here resolved per column because thr_f is the x-axis in this version).
# Returns the optimal thr_m, the equalising shift thr_f - thr_m, and whether the
# residual is within the solution tolerance `tol` (+/-5%, as in the manuscript).
.optimal_points <- function(inp, strat, epsilon, tol = 0.05) {
  df <- .strategy_surface(inp, strat, epsilon)
  df <- df[df$thr_m <= df$thr_f & is.finite(df$ratio), , drop = FALSE]
  if (!nrow(df))
    return(data.frame(thr_f = numeric(), thr_m_opt = numeric(), ratio = numeric(),
                      shift = numeric(), within_tol = logical(), epsilon = numeric()))
  do.call(rbind, lapply(split(df, df$thr_f), function(d) {
    j <- which.min(abs(d$ratio - 1))
    data.frame(thr_f = d$thr_f[1], thr_m_opt = d$thr_m[j], ratio = d$ratio[j],
               shift = d$thr_f[1] - d$thr_m[j],
               within_tol = abs(d$ratio[j] - 1) <= tol, epsilon = epsilon)
  }))
}

# ---- Shared visual style -----------------------------------------------------
# Clean, BioRender-inspired look applied to every figure: white "card" panels
# with a light border, bold centred titles, soft axis lines, and a faint
# horizontal guide grid. Soft blue = male, soft pink = female throughout.
col_male   <- "#7FB9DD"
col_female <- "#F3B7D0"

# One shared palette so every figure reads as one analysis.
IFR_SEX_LIGHT <- c(m = "#9CD0EE", f = "#F6B6CE")   # fit ensemble / soft sex colours
IFR_SEX_DARK  <- c(m = "#15539e", f = "#b32219")   # consensus steps, error bars, points
HEAT_BLUE     <- "#2166AC"                          # "males die more" end of the diverging scale
HEAT_RED      <- "#B2182B"                          # "females die more" end
# Accent for the equalising-shift profile (top panels) and the optimal markers.
STRAT_ACCENT  <- "#222222"
# Single colour for every top shift-profile panel: with a shared colour and a
# shared y-axis, the differing panel shapes make the strategies' contrast obvious.
SHIFT_COL     <- "#6A51A3"
# Per-strategy tint, retained for the detailed/efficiency views.
STRAT_TINT <- c(eqDR = "#1A7DB6", eqYPLL = "#6A51A3", eqPD = "#D7642A", eqRisk = "#2C8C5A")

pub_theme <- function(base_size = 22) {
  theme_classic(base_size = base_size) +
    theme(
      text               = element_text(colour = "grey15"),
      plot.title         = element_text(face = "bold", size = base_size + 5, hjust = 0.5,
                                        margin = margin(b = base_size / 2), colour = "grey15"),
      plot.subtitle      = element_blank(),
      plot.caption       = element_blank(),
      axis.title         = element_text(face = "bold", size = base_size + 1, colour = "grey20"),
      axis.text          = element_text(size = base_size - 1, colour = "grey35"),
      axis.line          = element_line(colour = "grey40", linewidth = 0.7),
      axis.ticks         = element_line(colour = "grey40", linewidth = 0.7),
      strip.text         = element_text(face = "bold", size = base_size, colour = "grey15"),
      strip.background   = element_rect(fill = "grey93", colour = NA),
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 2),
      legend.key         = element_blank(),
      legend.key.height  = unit(1.6, "lines"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.background    = element_rect(fill = "white", colour = "grey80", linewidth = 0.6),
      plot.margin        = margin(14, 16, 14, 14)
    )
}

# Every figure is written to its own SVG file (SVG, unlike PDF, is single-page).
# Uses svglite when available (cleanest editable SVG), else base grDevices::svg.
.open_svg <- function(file, width, height) {
  if (requireNamespace("svglite", quietly = TRUE))
    svglite::svglite(file, width = width, height = height)
  else
    svg(file, width = width, height = height)
}

# Render one plot or grob to a standalone SVG file.
.save_svg <- function(p, file, width, height) {
  .open_svg(file, width, height)
  on.exit(dev.off())
  if (inherits(p, c("gtable", "grob", "gList"))) grid::grid.draw(p) else print(p)
}

# Rectilinear outline of the equality band. For every cell within `tol` of
# equality, only the edges shared with a non-band cell (or the grid border) are
# drawn, so the boundary hugs cell edges exactly with no internal mesh -- this
# stays accurate when zoomed in. Returns unit segments, faceted by epsilon.
.band_outline <- function(df, tol) {
  empty <- data.frame(x = numeric(), y = numeric(), xend = numeric(),
                      yend = numeric(), epsilon = numeric())
  w <- df[is.finite(df$ratio) & abs(df$ratio - 1) <= tol, c("epsilon", "thr_f", "thr_m")]
  if (!nrow(w)) return(empty)
  inkey <- paste(w$epsilon, w$thr_f, w$thr_m)
  edge <- function(dx, dy, x1, y1, x2, y2) {
    keep <- !(paste(w$epsilon, w$thr_f + dx, w$thr_m + dy) %in% inkey)
    if (!any(keep)) return(empty)
    data.frame(x = w$thr_f[keep] + x1, y = w$thr_m[keep] + y1,
               xend = w$thr_f[keep] + x2, yend = w$thr_m[keep] + y2,
               epsilon = w$epsilon[keep])
  }
  rbind(edge( 1, 0,  0.5, -0.5,  0.5,  0.5),    # right edge
        edge(-1, 0, -0.5, -0.5, -0.5,  0.5),    # left edge
        edge( 0, 1, -0.5,  0.5,  0.5,  0.5),    # top edge
        edge( 0,-1, -0.5, -0.5,  0.5, -0.5))    # bottom edge
}

# A single strategy heatmap. Cells are coloured by log2(male metric / female
# metric), oriented so blue = males die more, red = females die more. The
# equality band (|ratio - 1| <= tol) is outlined along cell edges, and the
# thr_f = thr_m diagonal is a dashed reference line. Pass `lim` to share a
# symmetric colour scale across facets.
.strategy_heatmap <- function(df, meta, tol = 0.05, domain_only = FALSE,
                              lim = NULL, base_size = 18, optimal = NULL,
                              square = "coord", legend = TRUE, title = NULL) {
  if (domain_only) df <- df[df$thr_f >= df$thr_m, , drop = FALSE]
  df$logr <- log2(df$ratio)
  
  # symmetric colour limit, driven by the modelled domain (thr_f >= thr_m)
  if (is.null(lim)) {
    dom <- df$logr[df$thr_f >= df$thr_m & is.finite(df$logr)]
    lim <- if (length(dom)) as.numeric(stats::quantile(abs(dom), 0.99)) else NA
    if (!is.finite(lim) || lim <= 0) lim <- max(abs(df$logr[is.finite(df$logr)]))
  }
  # five evenly spaced breaks in log2 space (ends, quarter points, centre) so the
  # horizontal colourbar never crowds its tick labels
  brk <- c(-lim, -lim / 2, 0, lim / 2, lim)
  fmt_ratio <- function(r)                         # compact, fixed-width-ish labels
    ifelse(r >= 100, formatC(signif(r, 2), format = "fg", big.mark = ""),
           ifelse(r >= 10,  sprintf("%.0f", r),
                  ifelse(r >= 1,   sprintf("%.1f", r),
                         ifelse(r >= 0.1, sprintf("%.2f", r),
                                formatC(signif(r, 1), format = "fg")))))
  
  # blue = "males die more": high-ratio end for most metrics, low-ratio end for
  # prevented deaths (invert = TRUE), where higher mPD/fPD means males die less.
  lo <- if (isTRUE(meta$invert)) HEAT_BLUE else HEAT_RED
  hi <- if (isTRUE(meta$invert)) HEAT_RED  else HEAT_BLUE
  
  outline <- .band_outline(df, tol)
  
  p <- ggplot(df, aes(thr_f, thr_m)) +
    geom_raster(aes(fill = logr)) +
    geom_segment(data = outline, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "31",
                colour = "grey25", linewidth = 0.4)
  
  # optimal thr_m for each thr_f (within +/-5%): the equalising combinations
  if (!is.null(optimal) && nrow(optimal)) {
    op <- optimal[optimal$within_tol & is.finite(optimal$thr_m_opt), , drop = FALSE]
    if (nrow(op))
      p <- p + geom_point(data = op, aes(thr_f, thr_m_opt),
                          shape = 21, fill = "white", colour = STRAT_ACCENT,
                          size = 1.5, stroke = 0.5, inherit.aes = FALSE)
  }
  
  p <- p +
    scale_fill_gradient2(
      name = meta$lab, low = lo, mid = "white", high = hi,
      midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
      breaks = brk, labels = fmt_ratio(2^brk),
      guide = if (legend) "colourbar" else "none"
    ) +
    labs(x = expression(thr[f]), y = expression(thr[m]),
         title = if (is.null(title)) meta$title else title) +
    pub_theme(base_size) +
    theme(panel.grid = element_blank(), axis.line = element_blank(),
          panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.7))
  
  if (identical(square, "aspect"))
    p + coord_cartesian(xlim = range(df$thr_f), ylim = range(df$thr_m), expand = FALSE) +
    theme(aspect.ratio = 1)
  else if (identical(square, "none"))
    p + coord_cartesian(xlim = range(df$thr_f), ylim = range(df$thr_m), expand = FALSE)
  else
    p + coord_fixed(expand = FALSE)
}

# Top marginal panel: the equalising shift (thr_f - thr_m) versus the female
# threshold thr_f, restricted to the equalising optima. Every panel uses one
# colour and one shared y-axis (ylim), so the differing shapes make the
# strategies' contrast obvious. The x-axis is dropped because the heatmap
# directly below carries it; the panel is the same width as that heatmap.
.shift_profile <- function(opt, xlim, ylim, title = NULL, base_size = 15) {
  o <- opt[opt$within_tol & is.finite(opt$shift), , drop = FALSE]
  ggplot(o, aes(thr_f, shift)) +
    geom_line(colour = SHIFT_COL, linewidth = 1) +
    geom_point(colour = SHIFT_COL, size = 1.4) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = NULL, y = expression(thr[f] - thr[m]), title = title) +
    pub_theme(base_size) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.7),
          plot.background = element_blank(),
          plot.margin = margin(6, 6, 0, 6))
}

# One strategy panel: the shift profile sitting directly on top of its heatmap,
# sharing the thr_f axis and width, with the per-thr_f optimum marked on the
# heatmap and the colourbar collected to the right of the pair. heights and
# margins are tuned so the two read as a single plot.
.strategy_composite <- function(inp, strat, epsilon, opt, ylim,
                                tol_band = 0.05, base_size = 15) {
  meta <- .strategy_meta[[strat]]
  surf <- .strategy_surface(inp, strat, epsilon)
  xlim <- range(inp$age)
  p_top  <- .shift_profile(opt, xlim, ylim, title = meta$title, base_size = base_size)
  p_heat <- .strategy_heatmap(surf, meta, tol = tol_band, optimal = opt, title = "",
                              base_size = base_size, square = "none", legend = TRUE) +
    theme(plot.title   = element_blank(),
          legend.title = element_text(size = base_size - 2),
          legend.text  = element_text(size = base_size - 3),
          plot.margin  = margin(0, 6, 6, 6)) +
    guides(fill = guide_colourbar(barwidth = unit(0.6, "lines"), barheight = unit(7, "lines")))
  (p_top / p_heat) +
    plot_layout(heights = c(1, 4), guides = "collect") &
    theme(legend.position = "right")
}

plots.models <- function(detailed = FALSE, main_file = NULL,
                         eqDR_file = NULL, eqYPLL_file = NULL,
                         rF.vis = rF_vis_def) {
  inp     <- .strategy_inputs()
  eps_vis <- 1 - rF.vis
  tol     <- 0.05   # equality band drawn on the heatmaps: +/-5% (as in the manuscript)
  
  if (!detailed) {
    # ---- Main results (Figure 4): 2 x 2; each strategy = shift profile (top)
    #      over its heatmap, at efficiency eps_vis, with the equalising optimum
    #      marked. All top panels share one colour and one y-axis. --------------
    order4 <- c("eqRisk","eqPD","eqDR", "eqYPLL")
    opts   <- lapply(order4, function(s) .optimal_points(inp, s, eps_vis, tol = 0.05))
    names(opts) <- order4
    ymax <- max(c(0, unlist(lapply(opts, function(o) {
      v <- o$shift[o$within_tol & is.finite(o$shift)]; if (length(v)) max(v) else 0 }))))
    ylim <- c(0, ymax * 1.08)
    
    comps <- lapply(order4, function(s)
      .strategy_composite(inp, s, eps_vis, opt = opts[[s]], ylim = ylim,
                          tol_band = tol, base_size = 15))
    
    G <- wrap_plots(comps, ncol = 2) +
      plot_annotation(
        title = sprintf("Equality objective outcome metrics   (efficiency \u03b5 = %.2f)", eps_vis),
        theme = theme(plot.title = element_text(face = "bold", size = 24, hjust = 0.5,
                                                margin = margin(b = 6)),
                      plot.background = element_rect(fill = "white", colour = NA)))
    .save_svg(G, main_file, width = 15, height = 16)
    
  } else {
    # ---- Efficiency sensitivity (Figures S4/S5): the epsilon-dependent
    #      strategies faceted over efficiency. ---------------------------------
    eps_grid <- sort(unique(1 - rF_grid_def))
    
    facet_strategy <- function(s) {
      meta <- .strategy_meta[[s]]
      df   <- do.call(rbind, lapply(eps_grid, function(e) .strategy_surface(inp, s, e)))
      lr   <- log2(df$ratio)
      dom  <- lr[df$thr_f >= df$thr_m & is.finite(lr)]
      lim  <- as.numeric(stats::quantile(abs(dom), 0.99))
      .strategy_heatmap(df, meta, tol = tol, lim = lim, base_size = 16) +
        facet_wrap(~ epsilon,
                   labeller = labeller(epsilon = function(x) paste0("\u03b5 = ", x)))
    }
    
    .save_svg(facet_strategy("eqDR"),   eqDR_file,   width = 16, height = 15)
    .save_svg(facet_strategy("eqYPLL"), eqYPLL_file, width = 16, height = 15)
  }
}

plots.all <- function(profiles_file, efficiency_file,
                      rF_summary = plots_all_rF_def,
                      covariate  = plots_all_covariate_def) {  # cross-country colour covariate
  covariate <- match.arg(covariate, c("sex_ratio", "dr_ratio", "e0_gap"))
  order4       <- c("eqRisk","eqPD","eqDR", "eqYPLL")
  strat_titles <- vapply(order4, function(s) .strategy_meta[[s]]$title, character(1))
  eps_vis      <- 1 - rF_summary
  
  cov_lab <- c(sex_ratio = "elderly sex ratio  (F:M, 65+)",
               dr_ratio  = "male / female death-rate ratio",
               e0_gap    = "life-expectancy gap  E0(f) - E0(m)  [yrs]")[[covariate]]
  
  ## ---- per-country covariates -------------------------------------------------
  ## e0_gap comes from Population.csv; dr_ratio (crude no-protection male/female
  ## death-rate ratio) is captured during the multi-country loop into country.info.
  e0 <- read.table("Data/Population.csv", sep = " ", dec = ".", header = TRUE,
                   stringsAsFactors = FALSE) %>%
    mutate(e0_gap = ex.Female - ex.Male) %>%
    filter(Age == 0) %>%
    transmute(pop = c_code, year = Year, e0_gap)
  attribs <- if (exists("country.info") && nrow(country.info))
    suppressMessages(dplyr::full_join(country.info, e0)) else e0
  if (!covariate %in% names(attribs))
    stop(sprintf("covariate '%s' is unavailable; run the multi-country loop so country.info is populated.", covariate))
  
  ## ---- tidy optimum table from the stacked multi-country results --------------
  collect <- function(strategies, rF_values) {
    base <- minAge.results.all %>% filter(rF %in% rF_values) %>% tidyr::unnest(r)
    do.call(rbind, lapply(strategies, function(s) {
      m <- .strategy_meta[[s]]
      data.frame(
        pop      = base$pop,
        year     = base$year,
        epsilon  = 1 - base$rF,
        mMA      = base$mMA,
        strategy = factor(m$title, levels = strat_titles),
        fMA      = base[[paste0(s, ".fMA")]],
        ratio    = base[[paste0(s, ".", m$m, "_", m$f)]],
        vis      = base[[paste0(s, ".vis")]],
        shift    = base[[paste0(s, ".fMA")]] - base$mMA,
        stringsAsFactors = FALSE
      )
    }))
  }
  add_cov <- function(d) {
    d <- suppressMessages(dplyr::left_join(d, attribs[, c("pop", "year", covariate)]))
    d$cov <- d[[covariate]]
    d
  }
  
  ## Equity is only meaningful where it is actually attainable: keep cells whose
  ## optimum ratio lies inside the +/-5% equality band (matching the heatmaps) and drop
  ## the rest, so no misleading "closest argmin" shift values are drawn.
  tol    <- 0.05
  attain <- function(d) d[is.finite(d$ratio) & abs(d$ratio - 1) <= tol, , drop = FALSE]
  
  ## ---- Figure S6: shift profiles, 2 x 2 by strategy, at eps_vis -------------
  dfp <- attain(add_cov(collect(order4, rF_summary)))
  pA <- ggplot(dfp, aes(mMA, shift, group = pop, colour = cov)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ strategy, ncol = 2) +
    scale_colour_viridis_c(name = cov_lab) +
    labs(x = expression(thr[m]), y = expression(thr[f] - thr[m])) +
    pub_theme(20)
  .save_svg(pA, profiles_file, width = 16, height = 13)
  
  ## ---- Figure S7: efficiency small-multiples for the eps-dependent strategies
  dfe <- attain(add_cov(collect(c("eqDR", "eqYPLL"), sort(unique(rF_grid_def)))))
  pC <- ggplot(dfe, aes(mMA, shift, group = pop, colour = cov)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.1) +
    facet_grid(strategy ~ epsilon,
               labeller = labeller(epsilon = function(x) paste0("\u03b5 = ", x))) +
    scale_colour_viridis_c(name = cov_lab) +
    labs(x = expression(thr[m]), y = expression(thr[f] - thr[m])) +
    pub_theme(15)
  .save_svg(pC, efficiency_file, width = 25, height = 10)
  
  invisible(list(profiles = dfp, efficiency = dfe))
}

# ---- IFR fit figure ----------------------------------------------------------
# Conceptually early: a single figure comparing the published "consensus" IFR
# estimates (median +/- 95% credible interval) with the age-based model -- both
# the continuous parametric fit (soft bands) and the per-interval model
# estimates that match the consensus. Requires a completed single-country
# agebased() run (fx.*, mort.pred.*, meanage.*, ints, IFR_optim).
# Build the IFR fit figure from the globals left by a single agebased() run.
# callouts = TRUE draws the four explanatory annotations + gender pictograms
# (for the standalone Figure); FALSE keeps a clean panel for the sensitivity
# grid, with an optional `title`. The model estimate per 5-year category is the
# population-weighted mean IFR of the fit; the error bar is the FULL RANGE of
# that weighted mean across the stored fit replicates (it = 0 plus the IFR-draw
# refits), i.e. the spread of the quantity the cost function compares.
.ifr_fit_plot <- function(base_size = 18, callouts = TRUE, title = NULL) {
  n_fit <- length(ints) - 1
  
  steps <- rbind(
    data.frame(x = xcoord, lower = ycoord.male.lower,   median = ycoord.male.median,   upper = ycoord.male.upper,   sex = "m"),
    data.frame(x = xcoord, lower = ycoord.female.lower, median = ycoord.female.median, upper = ycoord.female.upper, sex = "f"))
  steps$col <- IFR_SEX_DARK[steps$sex]
  
  fits <- dplyr::filter(fx.age.age.results, gender %in% c("m", "f"))
  fits$col <- IFR_SEX_LIGHT[fits$gender]
  
  ## per-interval, per-replicate weighted-mean IFR; error bar = min..max over replicates
  dem       <- dem.gen(population_def, year_def, FALSE)
  intervals <- intervals_from_ints(ints)
  iseq      <- seq_len(n_fit)
  wmean <- function(v, w) sapply(iseq, weighted_interval_mean, values = v, weights = w, intervals = intervals)
  ens <- do.call(rbind, lapply(split(fits, list(fits$gender, fits$it), drop = TRUE), function(d) {
    g <- d$gender[1]; w <- if (g == "m") dem$Male else dem$Female
    data.frame(cat = iseq, gender = g, it = d$it[1], pred = wmean(d$value[order(d$age)], w))
  }))
  est <- ens %>%
    dplyr::group_by(gender, cat) %>%
    dplyr::summarise(pred  = pred[it == 0],
                     lower = min(pred), upper = max(pred), .groups = "drop") %>%
    dplyr::mutate(age = ifelse(gender == "m", meanage.male[cat], meanage.female[cat]),
                  col = IFR_SEX_DARK[gender])
  
  p <- ggplot() +
    geom_line(data = fits, aes(age, value, colour = col, group = interaction(gender, it)),
              linewidth = 0.35, alpha = 0.5) +
    geom_line(data = steps, aes(x, lower,  colour = col, group = sex), linewidth = 0.45, linetype = 2) +
    geom_line(data = steps, aes(x, upper,  colour = col, group = sex), linewidth = 0.45, linetype = 2) +
    geom_line(data = steps, aes(x, median, colour = col, group = sex), linewidth = 0.95) +
    geom_errorbar(data = est, aes(age, ymin = lower, ymax = upper, colour = col), width = 1.7, linewidth = 0.6) +
    geom_point(data = est, aes(age, pred, colour = col), size = 1.8) +
    scale_colour_identity() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    coord_cartesian(xlim = c(18, 100), ylim = c(3e-5, 5e-1)) +
    labs(x = "Age", y = "IFR", title = title) +
    pub_theme(base_size)
  
  if (callouts) {
    p <- p +
      annotate("text", x = 42, y = 4e-3, label = "\u2642", colour = "#7FB9DD", size = 18) +
      annotate("text", x = 60, y = 7e-4, label = "\u2640", colour = "#F3B7D0", size = 18) +
      annotate("text", x = 50, y = 0.34, label = "age-based model:\nparametric IFR estimates",
               fontface = "bold", size = 5, hjust = 0.5, lineheight = 0.95) +
      annotate("segment", x = 63, y = 0.30, xend = 86, yend = 0.16, arrow = grid::arrow(length = grid::unit(0.02, "npc"))) +
      annotate("text", x = 31, y = 0.05, label = "age-based model:\nIFR model estimates\nthat match the consensus",
               fontface = "bold", size = 5, hjust = 0.5, lineheight = 0.95) +
      annotate("segment", x = 49, y = 0.05, xend = 72, yend = 0.045, arrow = grid::arrow(length = grid::unit(0.02, "npc"))) +
      annotate("text", x = 70, y = 5.5e-3, label = "consensus IFR estimate", size = 4.4, hjust = 0.5) +
      annotate("segment", x = 77, y = 7.5e-3, xend = 81, yend = 0.052, arrow = grid::arrow(length = grid::unit(0.02, "npc"))) +
      annotate("text", x = 93, y = 1.4e-2, label = "95% credible\nintervals", size = 4.4, hjust = 0.5, lineheight = 0.95) +
      annotate("segment", x = 91, y = 2.2e-2, xend = 89, yend = 0.085, arrow = grid::arrow(length = grid::unit(0.02, "npc")))
  }
  p
}

plots.diagnostics <- function(fit_file, ratio_file) {
  .save_svg(.ifr_fit_plot(base_size = 18, callouts = TRUE), fit_file, width = 11, height = 13)
  
  ## raw sex-specific IFR ratio by age for both published sources (step plot)
  ifr_ratio_df <- function(src) {
    d <- IFR_optim[[src]]
    r <- d$male$median / d$female$median
    data.frame(age = d$ints, ratio = c(r, r[length(r)]), source = src)
  }
  ratio_dat <- rbind(ifr_ratio_df("ODriscoll"), ifr_ratio_df("Salje"))
  p_ratio <- ggplot(ratio_dat, aes(age, ratio, colour = source)) +
    geom_step(linewidth = 1.5, direction = "hv") +
    scale_colour_manual(NULL,
                        values = c(ODriscoll = "royalblue", Salje = "firebrick"),
                        labels = c(ODriscoll = "O'Driscoll et al.", Salje = "Salje et al.")) +
    coord_cartesian(xlim = c(5, 90), ylim = c(1, 3)) +
    labs(x = "Age", y = "male IFR / female IFR",
         title = "Ratio of sex-specific IFR estimates by age") +
    pub_theme(18)
  .save_svg(p_ratio, ratio_file, width = 12, height = 9)
  
  invisible(NULL)
}

# ---- Sensitivity grids (three figures from one function) ---------------------
# Re-runs the single-country fit for every (IFR source x fit function) scenario
# and assembles, from the captured panels:
#   <prefix>1_fit_sensitivity.svg : IFR fit for each scenario (sources x fits)
#   <prefix>2_analytical_ODriscoll.svg : 4 strategies x fits, ODriscoll IFR
#   <prefix>3_analytical_Salje.svg     : 4 strategies x fits, Salje IFR
# This loop overwrites the fit globals (and is slow, refitting every scenario), so it
# is meant to run after the main outputs are written. Failed fits get a
# placeholder panel so the grids always assemble.
.blank_panel <- function(title = "") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "fit unavailable", colour = "grey55", size = 4.5) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13, colour = "grey15"),
          plot.background = element_rect(fill = "white", colour = "grey85"),
          plot.margin = margin(8, 10, 8, 10))
}

plots.sensitivity <- function(fit_file, od_file, salje_file,
                              fits    = c("exp2", "exp3", "G"),
                              sources = c("ODriscoll", "Salje"),
                              n_fits  = 10) {
  fit_label <- c(exp2 = "exp2", exp3 = "exp3", G = "Gompertz")
  order4    <- c("eqDR", "eqYPLL", "eqPD", "eqRisk")
  ncol_f    <- length(fits)
  eps_vis   <- 1 - rF_vis_def
  title_th  <- theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5, colour = "grey15"),
                     plot.background = element_rect(fill = "white", colour = NA))
  
  ## ---- run every (source, fit) once; capture panels + surfaces --------------
  ifr_panels <- list(); heat_surf <- list()
  for (src in sources) for (fit in fits) {
    key <- paste(src, fit, sep = "|")
    ok <- FALSE; att <- 1
    while (!ok && att <= 3) {
      att <- att + 1
      ok <- tryCatch({ agebased(optFunction = fit, optIFR = src,
                                rF_grid = rF_grid_def, n_fits = n_fits); TRUE },
                     error = function(e) FALSE)
    }
    ttl <- sprintf("%s \u2014 %s", src, fit_label[[fit]])
    if (!ok) { ifr_panels[[key]] <- .blank_panel(ttl); heat_surf[[key]] <- NULL; next }
    ifr_panels[[key]] <- .ifr_fit_plot(base_size = 12, callouts = FALSE, title = ttl)
    inp <- .strategy_inputs()
    sl  <- lapply(order4, function(s) .strategy_surface(inp, s, eps_vis))
    names(sl) <- order4; heat_surf[[key]] <- sl
  }
  
  ## ---- Figure S1: IFR fit sensitivity (rows = sources, cols = fits) ----------
  s1 <- do.call(c, lapply(sources, function(src)
    lapply(fits, function(fit) ifr_panels[[paste(src, fit, sep = "|")]])))
  S1 <- wrap_plots(s1, ncol = ncol_f)
  .save_svg(S1, fit_file, width = 4.6 * ncol_f, height = 6.0 * length(sources))
  
  ## ---- Figures S2/S3: analytical sensitivity (rows = strategies, cols = fits)-
  analytical_grid <- function(src, out_file) {
    panels <- list(); k <- 0
    for (si in seq_along(order4)) {
      s <- order4[si]
      lrs <- unlist(lapply(fits, function(fit) {
        sf <- heat_surf[[paste(src, fit, sep = "|")]]; if (is.null(sf)) return(numeric())
        lr <- log2(sf[[s]]$ratio); lr[sf[[s]]$thr_f >= sf[[s]]$thr_m & is.finite(lr)] }))
      lim <- if (length(lrs)) as.numeric(stats::quantile(abs(lrs), 0.99)) else 1
      for (ci in seq_along(fits)) {
        fit <- fits[ci]; sf <- heat_surf[[paste(src, fit, sep = "|")]]
        ttl <- if (si == 1) fit_label[[fit]] else ""
        if (is.null(sf)) { k <- k + 1; panels[[k]] <- .blank_panel(ttl); next }
        p <- .strategy_heatmap(sf[[s]], .strategy_meta[[s]], tol = 0.05, lim = lim,
                               base_size = 12, square = "aspect", legend = FALSE, title = ttl)
        if (ci == 1) p <- p + labs(y = .strategy_meta[[s]]$title)
        else        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
        if (si != length(order4))
          p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                         axis.ticks.x = element_blank())
        k <- k + 1; panels[[k]] <- p
      }
    }
    G <- wrap_plots(panels, ncol = ncol_f)
    .save_svg(G, out_file, width = 4.3 * ncol_f, height = 4.8 * length(order4))
  }
  analytical_grid("ODriscoll", od_file)
  analytical_grid("Salje",     salje_file)
  invisible(NULL)
}

# Commands----
age <- 20:99               # Analysis age range. Use a contiguous integer vector; IFR intervals are filtered to this range.
bestMethod <- T            # Kept for compatibility; currently not used elsewhere in this script.
optIFR_def <- "ODriscoll"  # IFR source options: names(IFR_optim) -> "ODriscoll", "Salje".
population_def <- "GBR_NP" # Single-country run. Options: any c_code in Data/Population.csv, e.g. "GBR_NP", "CZE", "BGR".
year_def<-2018             # Single-country run year. Must exist for population_def in Data/Population.csv.
optf="exp2"              # Fit options: names(functions_optim) -> "exp2", "exp3", "G" (Gompertz).

# Protection levels in detailed plots
rF_grid_def <- seq(0.1, 0.9, by = 0.1) # Protection-level grid tested by agebased(); alternatives: c(0.1, 1 / 3, 0.5, 2 / 3, 0.9), any numeric vector in [0, 1].
rF_vis_def <- 0.5        # Protection level shown by plots.models(); should usually be one of rF_grid_def.
plots_all_rF_def <- 0.5  # Protection level used by plots.all() cross-country summaries; should usually be one of rF_grid_def.
n_fits_def <- 100         # IFR fit figure: number of extra fits to sampled IFR draws (the thin fit-ensemble lines). 0 disables.
# Cross-country summary plots
plots_all_covariate_def <- "sex_ratio" # plots.all() colour covariate: "sex_ratio" (elderly F:M, 65+), "dr_ratio" (male/female death-rate ratio) or "e0_gap".
# Sensitivity grids (Figures S1-S3). SLOW: refits every IFR-source x fit-function scenario.
run_sensitivity <- TRUE  # FALSE skips the three sensitivity grids entirely.
n_fits_sens     <- 100    # IFR-draw refits per scenario panel in the S1 fit-sensitivity grid.
sens_fits       <- c("exp2", "exp3", "G")             # fit functions compared (columns): exp2, exp3, Gompertz.
sens_sources    <- c("ODriscoll", "Salje")            # IFR sources compared.

# ---- Output file names, numbered to match the manuscript figure list --------
# Figure 2 (UK counterfactual model design) is a BioRender schematic and is not
# produced by this script.
fig_src <- optIFR_def; fig_fit <- optf
fig1_ratio_svg     <- "Results/Figure1_IFR_ratio_ODriscoll_Salje.svg"                       # Fig 1: O'Driscoll vs Salje IFR ratio
fig3_fit_svg       <- paste0("Results/Figure3_IFR_fit_", fig_src, ".svg")                    # Fig 3: IFR fit
fig4_main_svg      <- paste0("Results/Figure4_analytical_main_", fig_src, "_", fig_fit, ".svg")            # Fig 4: analytical main results
figS1_fitsens_svg  <- "Results/FigureS1_IFR_fit_sensitivity.svg"                             # Fig S1: IFR fit sensitivity
figS2_od_svg       <- "Results/FigureS2_analytical_sensitivity_ODriscoll.svg"                # Fig S2: analytical sensitivity (ODriscoll)
figS3_salje_svg    <- "Results/FigureS3_analytical_sensitivity_Salje.svg"                    # Fig S3: analytical sensitivity (Salje)
figS4_eqDR_svg     <- paste0("Results/FigureS4_efficiency_sensitivity_eqDR_",  fig_src, "_", fig_fit, ".svg")   # Fig S4: efficiency sensitivity, Mortality Equality
figS5_eqYLL_svg    <- paste0("Results/FigureS5_efficiency_sensitivity_eqYLL_", fig_src, "_", fig_fit, ".svg")   # Fig S5: efficiency sensitivity, Life-year Equality
figS6_demog_svg    <- paste0("Results/FigureS6_demographic_sensitivity_", fig_src, "_", fig_fit, ".svg")        # Fig S6: demographic sensitivity
figS7_demogeff_svg <- paste0("Results/FigureS7_demographic_efficiency_sensitivity_", fig_src, "_", fig_fit, ".svg")  # Fig S7: demographics x efficiency
workspace_file     <- paste0("Results/4strategies_", fig_src, "_", fig_fit, "_20260722.RData")
results_xlsx       <- paste0("Results/analytical_results_", fig_src, "_", fig_fit, ".xlsx")  # curated results / figure-data workbook

# ============================================================================
# Results export -- curated workbook of the key analytical results and the
# data underlying the figures (mirrors the counterfactual script's summary).
# Re-fits the reference IFR source x fit so it is robust to placement after the
# multi-country loop / sensitivity grids (which overwrite the fit globals); the
# cross-country demographic summary is read from minAge.results.all in memory.
# ============================================================================
export_analytical_results <- function(file = results_xlsx,
                                      src = optIFR_def, fit = optf, eps_ref = rF_vis_def,
                                      ref_pop = "GBR_NP", ref_year = 2018,
                                      thr_show = seq(50, 90, 5), thr_m_ref = 65,
                                      thr_m_ill = 65, thr_unif = 65) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) { message("openxlsx unavailable; skipping results export."); return(invisible()) }
  message("Building analytical results workbook ...")
  population_def <<- ref_pop; year_def <<- ref_year
  agebased(optFunction = fit, optIFR = src, rF_grid = rF_grid_def, n_fits = 0)  # restore reference fit globals
  inp       <- .strategy_inputs()
  ets_shift <- as.numeric(sse.age.results$shift[1])     # Individual Risk Equality constant age shift (years)
  a0        <- age[1]
  idx       <- function(A) A - a0 + 1
  smap      <- c(eqDR = "Mortality Equality", eqYPLL = "Life-year Equality",
                 eqPD = "Benefit Equality", eqRisk = "Individual Risk Equality")
  loc <- function(strat, eps) { op <- .optimal_points(inp, strat, eps, tol = 0.05)
  op[op$within_tol & is.finite(op$shift), c("thr_f", "thr_m_opt", "shift"), drop = FALSE] }
  
  ## ---- config -------------------------------------------------------------
  config <- tibble(parameter = c("IFR source", "IFR fit function", "analysis age range",
                                 "reference population", "population year", "reference efficiency (epsilon)",
                                 "efficiency grid", "tolerance (|ratio-1|)", "strategy label = code"),
                   value = c(src, fit, paste0(min(age), "-", max(age)), ref_pop, ref_year, eps_ref,
                             paste(rF_grid_def, collapse = ", "), "0.05",
                             "Mortality Equality=eqDR; Life-year Equality=eqYPLL; Benefit Equality=eqPD; Individual Risk Equality=eqRisk"))
  
  ## ---- helpers ------------------------------------------------------------
  ## Summarise the equalising shift over the within-tolerance parity locus.
  ## EMS / EYS solutions exist only at high female thresholds, so a range is
  ## reported rather than a single value, and lookups never extrapolate.
  ssum <- function(strat, eps) { op <- loc(strat, eps)
  if (!nrow(op)) return(tibble(shift_min = NA_real_, shift_median = NA_real_, shift_max = NA_real_,
                               n_threshold_solutions = 0L, thr_f_from = NA_real_, thr_f_to = NA_real_))
  tibble(shift_min = round(min(op$shift), 2), shift_median = round(median(op$shift), 2),
         shift_max = round(max(op$shift), 2), n_threshold_solutions = nrow(op),
         thr_f_from = min(op$thr_f), thr_f_to = max(op$thr_f)) }
  at_thr <- function(strat, eps, x) { op <- loc(strat, eps)
  if (nrow(op) < 2) return(NA_real_); approx(op$thr_f, op$shift, x, rule = 1)$y }  # NA outside the locus
  
  ## ---- key scalar estimates ----------------------------------------------
  ages_all  <- age[age >= 30 & age <= 90]
  ratio_all <- fx.m.age[idx(ages_all)] / fx.f.age[idx(ages_all)]
  ## IFR exponential slope b = d ln(IFR)/d age (exact for the exp2 baseline)
  ifr_slope <- function(fx) (log(fx[idx(90)]) - log(fx[idx(30)])) / 60
  b_f <- ifr_slope(fx.f.age); b_m <- ifr_slope(fx.m.age)
  key_estimates <- tibble(
    quantity = c("Individual Risk Equality age shift (years, males earlier)",
                 "IFR exponential slope b, female (per year of age)",
                 "IFR exponential slope b, male (per year of age)",
                 "male/female IFR ratio (geometric mean, ages 30-90)",
                 "male/female IFR ratio at age 30", "male/female IFR ratio at age 50",
                 "male/female IFR ratio at age 70", "male/female IFR ratio at age 90"),
    value = round(c(ets_shift, b_f, b_m, exp(mean(log(ratio_all))),
                    fx.m.age[idx(30)]/fx.f.age[idx(30)], fx.m.age[idx(50)]/fx.f.age[idx(50)],
                    fx.m.age[idx(70)]/fx.f.age[idx(70)], fx.m.age[idx(90)]/fx.f.age[idx(90)]), 4))
  
  ## ---- per-strategy equalising-shift summary at reference epsilon (Figure 4)
  strategy_shift_summary <- map_dfr(names(smap), function(s) {
    if (s == "eqRisk")
      tibble(strategy = smap[[s]], shift_min = ets_shift, shift_median = ets_shift,
             shift_max = ets_shift, n_threshold_solutions = NA_integer_,
             thr_f_from = NA_real_, thr_f_to = NA_real_, epsilon = eps_ref)
    else bind_cols(tibble(strategy = smap[[s]]), ssum(s, eps_ref)) %>% mutate(epsilon = eps_ref) })
  
  ## ---- IFR fit (Figures 1, 3) --------------------------------------------
  ages_show <- seq(max(25, min(age)), min(95, max(age)), by = 5)
  IFR_fit <- tibble(age = ages_show,
                    IFR_male   = fx.m.age[idx(ages_show)],
                    IFR_female = fx.f.age[idx(ages_show)],
                    male_female_ratio = fx.m.age[idx(ages_show)] / fx.f.age[idx(ages_show)])
  
  ## ---- strategy shifts at representative female thresholds (Figure 4) -----
  strategy_shifts <- map_dfr(names(smap), function(s) {
    if (s == "eqRisk")
      tibble(strategy = smap[[s]], thr_f = thr_show, shift = ets_shift, epsilon = eps_ref)
    else tibble(strategy = smap[[s]], thr_f = thr_show,
                shift = vapply(thr_show, function(x) at_thr(s, eps_ref, x), numeric(1)),
                epsilon = eps_ref) }) %>% mutate(thr_m_opt = thr_f - shift)
  
  ## ---- full parity loci (Figure 4 underlying data) -----------------------
  parity_loci <- map_dfr(c("eqDR", "eqYPLL", "eqPD"), function(s) {
    op <- .optimal_points(inp, s, eps_ref, tol = 0.05)
    tibble(strategy = smap[[s]], thr_f = op$thr_f, thr_m_opt = op$thr_m_opt,
           shift = op$shift, ratio = op$ratio, within_tol = op$within_tol) })
  
  ## ---- efficiency sensitivity (Figures S4, S5; eqDR & eqYPLL depend on eps)
  efficiency_sensitivity <- map_dfr(rF_grid_def, function(e)
    map_dfr(c("eqDR", "eqYPLL"), function(s)
      bind_cols(tibble(epsilon = e, strategy = smap[[s]]),
                ssum(s, e) %>% select(shift_min, shift_median, shift_max, n_threshold_solutions))))
  
  ## ---- demographic summary across countries (Figures S6, S7) -------------
  demographic_summary <- NULL
  if (exists("minAge.results.all") && nrow(minAge.results.all) &&
      exists("country.info") && nrow(country.info)) {
    mr <- minAge.results.all %>% filter(rF == eps_ref) %>% tidyr::unnest(r) %>%
      filter(mMA == thr_m_ref) %>%
      transmute(pop, year,
                `Mortality Equality`       = eqDR.fMA   - mMA,
                `Life-year Equality`             = eqYPLL.fMA - mMA,
                `Benefit Equality` = eqPD.fMA   - mMA,
                `Individual Risk Equality`   = eqRisk.fMA - mMA)
    demographic_summary <- suppressMessages(dplyr::left_join(
      country.info %>% transmute(pop, year, elderly_sex_ratio = sex_ratio, death_rate_ratio = dr_ratio), mr)) %>%
      mutate(thr_m = thr_m_ref, epsilon = eps_ref) %>% arrange(elderly_sex_ratio)
  }
  
  ## ---- illustration paragraph: fixed male threshold + uniform threshold ---
  ## (a) optimal female threshold per strategy at the fixed male threshold
  ##     (model solver output for the reference population, set by agebased() above)
  illust_thresholds <- tibble(strategy = character(), male_threshold = numeric(),
                              optimal_female_threshold = numeric())
  if (exists("minAge.results") && nrow(minAge.results)) {
    rr <- minAge.results %>% filter(rF == eps_ref, mMA == thr_m_ill) %>% tidyr::unnest(r)
    if (nrow(rr)) illust_thresholds <- tibble(
      strategy = c("Individual Risk Equality", "Benefit Equality", "Mortality Equality", "Life-year Equality"),
      male_threshold = thr_m_ill,
      optimal_female_threshold = c(rr$eqRisk.fMA[1], rr$eqPD.fMA[1], rr$eqDR.fMA[1], rr$eqYPLL.fMA[1]))
  }
  ## (b) per-sex outcomes under a uniform threshold applied to both sexes
  mm_u <- .threshold_metrics(inp, eps_ref); ju <- which(inp$age == thr_unif)
  illust_uniform <- tibble(
    metric = c("IFR at threshold", "mortality rate (per capita)",
               "proportion of deaths averted", "years of life lost (per capita)"),
    male   = c(mm_u$mRisk[ju], mm_u$mDR[ju], mm_u$mPD[ju], mm_u$mYPLL[ju]),
    female = c(mm_u$fRisk[ju], mm_u$fDR[ju], mm_u$fPD[ju], mm_u$fYPLL[ju])) %>%
    mutate(male_female_ratio = male / female)
  
  ## ---- assemble workbook --------------------------------------------------
  README <- tibble(sheet = c("config", "key_estimates", "IFR_fit", "strategy_shift_summary",
                             "strategy_shifts", "parity_loci", "efficiency_sensitivity",
                             "demographic_summary", "illustration"),
                   contents = c("Run parameters and the strategy label-to-code map.",
                                "Headline scalar estimates: Individual Risk Equality age shift, the male/female IFR exponential slopes (b), and IFR ratios.",
                                "Fitted single-year IFR by sex and the male/female ratio at representative ages (Figures 1, 3).",
                                "Per-strategy equalising shift summarised over the parity locus (min / median / max, with the female-threshold range), reference efficiency (Figure 4). EMS/EYS solutions exist only at high thresholds, hence a range.",
                                "Equalising shift (thr_f - thr_m) per strategy at representative female thresholds, reference efficiency; NA where no parity solution exists at that threshold (Figure 4).",
                                "Full parity locus per outcome-based strategy at reference efficiency (data underlying Figure 4).",
                                "Equalising-shift range (min / median / max over the locus) across efficiency for the efficiency-dependent strategies (Figures S4, S5).",
                                paste0("Per-country equalising shift per strategy at thr_m = ", thr_m_ref,
                                       ", reference efficiency, with the elderly (65+) sex ratio (Figures S6, S7)."),
                                paste0("Results-illustration paragraph: optimal female thresholds at male threshold = ", thr_m_ill,
                                       ", and the per-sex IFR / mortality rate / deaths averted / YLL under a uniform threshold = ", thr_unif,
                                       " (reference efficiency).")))
  wb <- openxlsx::createWorkbook(); add <- function(n, d) { openxlsx::addWorksheet(wb, n); openxlsx::writeData(wb, n, d) }
  add("README", README); add("config", config); add("key_estimates", key_estimates)
  add("IFR_fit", IFR_fit); add("strategy_shift_summary", strategy_shift_summary)
  add("strategy_shifts", strategy_shifts); add("parity_loci", parity_loci)
  add("efficiency_sensitivity", efficiency_sensitivity)
  if (!is.null(demographic_summary)) add("demographic_summary", demographic_summary)
  ## illustration sheet: two labelled blocks mirroring the Results paragraph
  openxlsx::addWorksheet(wb, "illustration")
  openxlsx::writeData(wb, "illustration",
                      paste0("Optimal female thresholds at male threshold = ", thr_m_ill, " (epsilon = ", eps_ref, ")"), startRow = 1)
  openxlsx::writeData(wb, "illustration", illust_thresholds, startRow = 2)
  r2 <- nrow(illust_thresholds) + 4
  openxlsx::writeData(wb, "illustration",
                      paste0("Per-sex outcomes under a uniform threshold = ", thr_unif, " applied to both sexes (epsilon = ", eps_ref, ")"), startRow = r2)
  openxlsx::writeData(wb, "illustration", illust_uniform, startRow = r2 + 1)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  message("  wrote ", file, " (", length(names(wb)), " sheets)")
  invisible(file)
}

agebased(optFunction = optf, rF_grid = rF_grid_def, n_fits = n_fits_def) ## long step; uses age, population_def, year_def, optIFR_def, optf, rF_grid_def. n_fits adds the IFR-figure fit ensemble.
plots.diagnostics(fit_file = fig3_fit_svg, ratio_file = fig1_ratio_svg)                       # Figure 3 (IFR fit) and Figure 1 (IFR ratio).
plots.models(detailed = F, main_file = fig4_main_svg, rF.vis = rF_vis_def)                    # Figure 4 (analytical main results).
plots.models(detailed = T, eqDR_file = figS4_eqDR_svg, eqYPLL_file = figS5_eqYLL_svg, rF.vis = rF_vis_def) # Figures S4 / S5 (efficiency sensitivity).

popvector <- read.table(file = "Data/Population.csv", sep = " ", dec = ".", header = T, stringsAsFactors = F)%>% # Multi-country run source.
  drop_na()%>%
  group_by(c_code)%>%
  summarize(max(Year))#%>% # Uses latest available year per selected country.
#filter(c_code %in% c("CZE","GBR_NP","BGR")) # Multi-country options: replace with any Data/Population.csv c_code values.

minAge.results.all<-tibble() # Stores combined multi-country age-threshold results.
country.info<-tibble()       # Stores per-country covariates for plots.all() (e.g. crude death-rate ratio).
for (P in 1:nrow(popvector)) {
  population_def<<-as.character(popvector[P,1])
  year_def<<-as.integer(popvector[P,2])
  
  r <- NULL
  attempt <- 1 # Retry control for stochastic/optimizer failures; loop allows up to 3 attempts.
  while( is.null(r) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      r <- agebased(optFunction = optf, rF_grid = rF_grid_def)
    )
  } 
  
  minAge.results.all<<-bind_rows(minAge.results.all,minAge.results)
  .ci_dem <- dem.gen(population_def, year_def, FALSE) # per-country demographic covariates
  .old    <- .ci_dem$Age >= 65                        # elderly cut-off for the population sex ratio
  country.info<<-bind_rows(country.info, tibble(
    pop = population_def, year = year_def,
    sex_ratio = sum(.ci_dem$Female[.old]) / sum(.ci_dem$Male[.old]),   # elderly (65+) female:male population ratio
    dr_ratio  = (sum(.ci_dem$Male   * fx.m.age) / sum(.ci_dem$Male)) /  # crude (no-protection) male/female death-rate ratio
      (sum(.ci_dem$Female * fx.f.age) / sum(.ci_dem$Female))))
}

plots.all(
  profiles_file   = figS6_demog_svg,
  efficiency_file = figS7_demogeff_svg,
  rF_summary = plots_all_rF_def,
  covariate  = plots_all_covariate_def
) # Figure S6 (demographic sensitivity) and Figure S7 (demographics x efficiency).

# Save workspace----
save.image(file = workspace_file) # Change file name or optf to keep separate analysis runs.

# Sensitivity grids (Figures S1-S3)----
# Run last: this refits every IFR-source x fit-function scenario (slow, and it
# overwrites the single-country fit globals). The single-country state has
# already been saved above. Reset to the canonical diagnostic country first.
if (isTRUE(run_sensitivity)) {
  population_def <- "GBR_NP"; year_def <- 2018
  plots.sensitivity(fit_file = figS1_fitsens_svg, od_file = figS2_od_svg, salje_file = figS3_salje_svg,
                    fits = sens_fits, sources = sens_sources, n_fits = n_fits_sens)
}

# Results workbook (key results + figure data)----
# Placed last so the cross-country demographic summary (minAge.results.all) is
# complete; the function re-fits the reference IFR source x fit internally.
export_analytical_results()
