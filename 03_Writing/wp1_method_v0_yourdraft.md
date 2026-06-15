# WP1 — Method (v0, original draft by HK)

> Status: **v0 (HK first-draft, 2026-05-13)**
> Author: HK
> Note: This is the first written draft straight from the head, intentionally
> rough. Subsequent versions (v1, v2, ...) will iterate on this base.

---

## What is currently lacking (self-identified)

- Model validation plan
- Specific dataset
- Geographic scope

## Draft text

I will develop an inference model with two layers reconstructing the past
history of chikungunya (2014–2025) and forwardly simulating future
outbreaks over the next 10 years to identify the regions suitable for
vaccine effectiveness evaluation.

First, I will use a renewal-based model within a Bayesian framework. In
this model framework, I will use Dataset 1 to select sites suitable for
chikungunya vaccine effectiveness evaluation during the next 10 years.

To overcome the unidentifiability problem of surveillance bias, I will
use Dataset 2 that provides cumulative seroprevalence which will be used
as an anchor of the immune population to date. I will inform this to
likelihood estimation and will infer state-level time-varying reporting
rates.

Using the set of inferred parameters (effective reproductive number, true
observed cases, time-varying reporting rates), I will do forward
simulation using the age-structured SEIR model. Based on the previously
developed model from my work, I will use future projected climatic
covariates over the next 10 years based on the Shared Socioeconomic
Pathway (SSP) scenarios from CMIP.

The outputs will be driven from Model 1 and Model 2, respectively. From
the Model 1, I will infer effective reproductive number, true observed
cases, time-varying reporting rates. From Model 2, the probability of
outbreaks based on outbreak decision threshold (see supplementary) over
the next 1, 3, 5, and 10 years, the total number of lab-confirmed
symptomatic cases that are required to conduct vaccine effectiveness
study, and rankings of the states will be derived.

### Datasets

**Dataset 1**: weekly lab-confirmed symptomatic cases, climatic
covariates (temperature, precipitation, and evaporation), chikungunya
environmental suitability map, serology survey data, and demographic
data.

**Dataset 2**: Seroprevalence data.

---

*end of v0 draft (HK)*
