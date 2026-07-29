# Towards equity-based public health measures in epidemics: a case study of COVID-19

Code accompanying **Zsichla L. & Müller V.**, *Towards equity-based public health measures in epidemics: a case study of COVID-19* (Institute of Biology, ELTE Eötvös Loránd University, Budapest). [preprint / journal — **ADD link/DOI**]

The repository contains the two analyses reported in the paper. The **threshold-shift framework** derives sex-specific age thresholds under four equity objectives from age- and sex-specific COVID-19 infection fatality ratios (IFRs). The **counterfactual reconstruction** reproduces the United Kingdom's 2020–2021 winter epidemic wave and vaccination campaign and tests those predictions under realistic epidemic conditions.

## Scripts

- **`COVID19_SexEqualityStrategies.R`** — Threshold-shift framework. Fits single-year IFR–age curves by sex (two-parameter exponential in the main analysis; three-parameter exponential and Gompertz as sensitivity), then solves for the female/male age thresholds that equalize each of four equity objectives — *Individual Risk Equality* (equal individual IFR), *Benefit Equality* (equal proportion of deaths prevented), *Mortality Equality* (equal per-capita mortality), and *Life-year Equality* (equal per-capita years of life lost). Includes the intervention-effectiveness and demographic sensitivity analyses. Reproduces Figures 1, 3, 4 and the associated supplementary figures, and writes a results workbook (`.xlsx`).

- **`COVID19_counterfactualAnalyses.R`** — Counterfactual UK vaccination model. Calibrates an age-structured, sex-neutral vaccination model to the observed rollout, then reconstructs infections and deaths under counterfactual sex-prioritized allocation (continuous prioritization parameter *s*) by back-calculation from observed weekly deaths, accounting for vaccine effectiveness, age structure, incidence, and the non-VOC→Alpha transition. Computes the prioritization weight at which each equity objective reaches parity, with Monte-Carlo, attack-rate, and rollout-timing sensitivity analyses. Reproduces Figure 5 and the associated supplementary figures, and writes a results workbook (`.xlsx`).

## Data

Input files are read from the paths defined near the top of each script (e.g. a `Data/` subdirectory). Sources (full references in the paper):

- Population counts and remaining life expectancy — Human Mortality Database
- Age- and sex-specific IFRs — O'Driscoll et al. (main), Salje et al. (sensitivity)
- Weekly UK COVID-19 deaths by age and sex — COVerAGE-DB
- Weekly UK vaccination coverage — UK COVID-19 Vaccine Surveillance Report (Week 38, 2021)
- SGTF (Alpha share) among deaths — Davies et al.
- Weekly incidence — UK Office for National Statistics

## Requirements

R (≥ 4.3.1) with:

```r
install.packages(c("tidyverse","data.table","optimx","patchwork","ggpattern",
                   "ggrepel","gridExtra","svglite","openxlsx","writexl","readxl",
                   "lubridate","scales","png"))
```

## Citation

> Zsichla L, Müller V. Towards equity-based public health measures in epidemics: a case study of COVID-19. [journal, year — **ADD**].

## Contact

Levente Zsichla — Institute of Biology, ELTE Eötvös Loránd University, Budapest, Hungary.
