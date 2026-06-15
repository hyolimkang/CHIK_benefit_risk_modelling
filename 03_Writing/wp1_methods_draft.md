# WP1 — Methods (DRAFT v0)

> Status: **v0 — intentionally imperfect, written to be reacted to.**
> Author: HK · Date: 2026-05-13
> Read once, mark up freely, then revise. Do not optimise on first read.

---

## 0. Working title

> A modular, climate-informed Bayesian framework for chikungunya outbreak
> reconstruction, vaccination strategy evaluation, and post-licensure
> effectiveness study design in Brazil — with transferability to Oropouche
> and other emerging arboviruses.

*(react: title too long / acronym needed / "framework" overused / etc.)*

---

## 1. Overall framework

We propose a **modular two-layer Bayesian framework** designed around four
operational deliverables:

1. **Reconstruction** of past chikungunya transmission across Brazil at
   municipality resolution, integrating fragmented surveillance and
   serological evidence;
2. **Probabilistic projection** of where and when future outbreaks are most
   likely to occur, under explicit climate scenarios and with calibrated
   uncertainty;
3. **Comparative evaluation** of alternative vaccination deployment
   strategies (mass, age-targeted, reactive stockpile, hybrid) in terms of
   infections, cases, DALYs averted, and number needed to vaccinate;
4. **Prioritisation** of post-licensure vaccine-effectiveness (VE)
   evaluation sites, defined as municipalities where laboratory-confirmed
   symptomatic cases are most likely to accrue with calibrated uncertainty.

The framework comprises:

- **Layer 1 — Bayesian Reconstruction Engine**: a renewal-equation–based
  semi-mechanistic model that jointly estimates time-varying transmission
  potential, latent incidence, susceptible pool dynamics, and reporting
  probability, with serology providing the identifiability anchor.
- **Layer 2 — Forward Simulation Engine**: an age-structured discrete-time
  SEIR model parameterised by Layer 1 posteriors, used to compare
  vaccination strategies under climate scenarios.

The two layers share inferred parameters but are computationally decoupled:
Layer 1 is fit once to historical data; Layer 2 runs forward simulations
conditional on Layer 1 posterior draws. This decoupling makes scenario
comparison tractable at municipality scale and is essential to the
*executable, interpretable* design constraint of the framework.

### 1.1 Why this architecture (vs alternatives)

Existing approaches typically belong to one of two paradigms:

(a) **Statistical / spatio-temporal hierarchical models** (Lowe et al,
2013 *Stat Med*; 2016 *PLoS Med*; 2021 *Lancet Planet Health*): excellent
at reconstructing past transmission with rich covariate effects, but do
not represent interventions mechanistically.

(b) **Mechanistic compartmental models** (Bhatt et al, 2013 *Nature*;
Ribeiro Dos Santos et al, 2025 *Nat Med*): represent interventions
naturally, but are typically calibrated to country-aggregated burden, miss
sub-national heterogeneity, and rarely propagate surveillance bias.

Recent precedents from COVID-19 (Abbott et al, 2020 *Wellcome Open Res*
[EpiNow2]; Flaxman et al, 2020 *Nature*) demonstrate that
**semi-mechanistic Bayesian renewal frameworks** can bridge these paradigms.
We adapt this approach to chikungunya, addressing its specific features
(vector-borne transmission, sparse serology, sporadic outbreak dynamics,
high sub-clinical fraction). Our key methodological contribution is the
**explicit coupling of serology likelihood with the renewal equation**,
which resolves the otherwise irreducible incidence-reporting confounding
identified as a major limitation in Ribeiro Dos Santos et al (2025).

*(react: too long? cut to one paragraph? expand on Brazil specifics?)*

---

## 2. Data sources

| Stream | Source | Coverage | Use |
|---|---|---|---|
| Notified + lab-confirmed cases | SINAN (SUS / Ministry of Health) | 5,570 municipalities, weekly, 2014–2024 | Layer 1 likelihood `Y_it` |
| Age-stratified seroprevalence | Curated meta-dataset: Cunha 2017, Brito 2020, Dias 2018, Donalisio 2017, PROVIT, CHIK-VAC | ~25 surveys in ~15 municipalities, 2014–2024 | Layer 1 anchoring likelihood |
| Temperature | ERA5 reanalysis, 0.25° | All Brazil, daily, 1990–present | Layer 1 R_it covariate |
| Precipitation | CHIRPS, 0.05° | All Brazil, daily, 1981–present | Layer 1 R_it covariate, lagged |
| Aedes suitability | Kraemer et al, 2019 | 5 km × 5 km global | Layer 1 R_it covariate |
| Population by age | IBGE censuses (2010, 2022) + intercensal interpolation | Municipality × age group, annual | Layer 2 denominator |
| Birth rates | DataSUS / IBGE | State × year | Susceptible pool replenishment |
| LACEN testing capacity | State health secretariats + Ministry of Health archives | State × year | Layer 1 ρ_it covariate |
| Hospitalisation rates | SINAN-SIH linkage (existing in-house) | Municipality × age | Layer 2 outcome conversion |

*(react: which data sources are not yet accessible? what permissions
needed? FioCruz / LACEN collaboration formalised?)*

---

## 3. Layer 1 — Bayesian Reconstruction Engine

### 3.1 Latent transmission process

For municipality `i ∈ {1,...,N}` and week `t ∈ {1,...,T}`, the latent
expected weekly symptomatic incidence `C_{i,t}` follows a discrete
renewal equation:

$$
\mathbb{E}[C_{i,t}] = R_{i,t} \cdot \sum_{s=1}^{S} C_{i,t-s} \, g_s
$$

where `g_s` is the discretised generation interval (Gamma, mean 14 days,
informed by Riswari et al, Salje et al) and `S` is the truncation lag
(8 weeks). The renewal equation operates on the latent expected scale;
stochasticity is carried by the observation layer (Section 3.4).

### 3.2 Transmission potential `R_{i,t}`

$$
\log R_{i,t} = \alpha_i + \beta^\top X_{i,t} + \log\!\left(\frac{S_{i,t-1}}{N_i}\right) + \varphi_t
$$

with:

- `α_i ~ Normal(μ_α, σ_α²)`: municipality random intercept
- `X_{i,t}`: standardised climate × vector covariates
- `S_{i,t-1} / N_i`: susceptible proportion (mass-action term)
- `φ_t`: AR(1) temporal random effect on log-R scale

`β^⊤ X_{i,t}` is parameterised following Ross-Macdonald
climate-driven `R_0` (Mordecai et al, 2017; Liu-Helmersson et al, 2014):

$$
\log R_0(T, R, V) = 2\log a(T) + \log b + \log c + \log(V/N) - \mu_v(T)\,\mathrm{EIP}(T) - \log \mu_v(T) - \log \gamma_h
$$

with biting rate `a(T)`, mosquito mortality `μ_v(T)`, extrinsic incubation
period `EIP(T)` as temperature-response functions fixed from literature;
`V/N` as a logistic function of Aedes suitability and lagged precipitation;
`b`, `c`, `γ_h` fixed with informative priors permitting modest deviation.
This grounds the transmission regression in vector biology while keeping
the inferential problem identifiable.

### 3.3 Susceptible pool dynamics

$$
S_{i,t+1} = S_{i,t} - C_{i,t} + b_i N_i \cdot \Delta t
$$

with `b_i` the municipality-specific weekly birth rate. Initial susceptibility:

$$
S_{i,0} = N_i \cdot (1 - \pi_i^{(0)})
$$

where `π_i^{(0)}` is informed by direct serology where available, spatial
smoothing otherwise, and pre-2014 outbreak history as a prior.

### 3.4 Observation process

Reported lab-confirmed cases:

$$
Y_{i,t} \sim \mathrm{NegBin}(\rho_{i,t} \cdot C_{i,t}, \phi_Y)
$$

$$
\mathrm{logit}\,\rho_{i,t} = \gamma_i + \eta^\top Z_{i,t}
$$

with `Z_{i,t}` including time-varying LACEN testing capacity, state-level
health system indicators, and outbreak-relative position (testing
saturation during peaks). `γ_i` is a municipality random effect.

### 3.5 Serology likelihood — methodological core

For each survey at municipality `i`, age group `a`, calendar year `τ_k`,
with `n` subjects and `k_pos` seropositive:

$$
k_\text{pos} \sim \mathrm{Binomial}(n, \pi_{i,a,\tau_k})
$$

$$
\pi_{i,a,\tau_k} = 1 - \exp\!\left(-\sum_{u \le \tau_k} \lambda_{i,a,u}\right)
$$

with per-susceptible force of infection from the renewal process:

$$
\lambda_{i,a,u} = \frac{C_{i,u}}{S_{i,u-1}/N_i} \cdot \mathbb{1}\{\text{individual was age } a \text{ at time } u\}
$$

**This explicit coupling of renewal incidence with serology likelihood is
the central methodological contribution of Layer 1.** It resolves the
otherwise unidentifiable ρ–C confounding by anchoring cumulative latent
incidence to direct serological evidence wherever the latter is available,
and propagates that anchoring spatially through the hierarchical structure.

### 3.6 Inference

Two implementations will be developed and compared:

- **Stan / cmdstanr**: full Bayesian via NUTS; used for state-level
  (27 units) reconstruction and for methodological validation on simulated
  data.
- **R-INLA**: integrated nested Laplace approximation; used for full
  municipality resolution (5,570 units) following Lowe et al precedents.

Consistency between implementations will be validated on a shared
state-level subset before deploying INLA for full-resolution analysis.

### 3.7 Layer 1 outputs

Posterior draws for: `R_{i,t}`, `C_{i,t}`, `S_{i,t}`, `ρ_{i,t}`, `β`, `η`,
`σ_α`, `σ_γ`, climate response coefficients, and forward-projected
analogues under three climate scenarios (RCP 4.5, RCP 8.5, observational
continuation).

*(react: too dense for grant? move into supplement? simplify equations?)*

---

## 4. Layer 2 — Age-Structured SEIR Forward Simulation

### 4.1 Compartmental structure

Discrete-time (weekly) by age group `a ∈ {<15, 15–39, 40–64, 65+}`:

$$
\begin{aligned}
S_{i,a,t+1} &= S_{i,a,t} - \lambda_{i,a,t} S_{i,a,t} - V_{i,a,t} + B_{i,a,t} \\
E_{i,a,t+1} &= E_{i,a,t} + \lambda_{i,a,t} S_{i,a,t} - \sigma E_{i,a,t} \\
I_{i,a,t+1} &= I_{i,a,t} + \sigma E_{i,a,t} - \gamma_h I_{i,a,t} \\
R_{i,a,t+1} &= R_{i,a,t} + \gamma_h I_{i,a,t} + V_{i,a,t} \cdot \mathrm{VE}_\text{inf}
\end{aligned}
$$

with FOI:

$$
\lambda_{i,a,t} = \beta_{i,t} \cdot \sum_{a'} M_{a,a'} \frac{I_{i,a',t}}{N_i}
$$

`β_{i,t}` drawn from Layer 1 posterior (vector dynamics implicit).
`M_{a,a'}` from Brazil-adapted Prem et al (2017) contact matrices.
`σ ≈ 1/3.5 days⁻¹`, `γ_h ≈ 1/5 days⁻¹` (Salje, Riswari literature).

### 4.2 Vaccination strategies

| ID | Strategy | Trigger | Coverage |
|---|---|---|---|
| S1 | No vaccination | — | counterfactual baseline |
| S2 | Mass campaign | one-off, 2025 | 50% age 12+ |
| S3 | Age-targeted | one-off, 2025 | 80% age 65+ |
| S4 | Reactive stockpile | when weekly cases > state threshold | 50% in affected health region |
| S5 | Climate-triggered | preceding high-suitability windows | 50% age 18+ |
| S6 | Hybrid | routine in high-FOI municipalities + stockpile elsewhere | mixed |
| S7 | Adaptive optimal *(exploratory)* | strategy chosen per municipality to minimise expected DALYs subject to budget | model-derived |

### 4.3 Outcome metrics

For every (strategy × municipality × posterior draw):
infections, symptomatic cases, hospitalisations (link to existing CHIK
hospitalisation model), deaths, DALYs averted, doses used, number needed
to vaccinate (NNV), ICER vs S1 (link to existing CE pipeline), and
probability of strategy dominance.

*(react: S7 adaptive optimal too ambitious for v1? push to future work?)*

---

## 5. Post-licensure VE Site Selection (uses Layer 1 only)

For each candidate municipality `i` and trial start date `t_0`:

$$
\mathrm{Accrual}_i^{(d)} = \sum_{t = t_0+1}^{t_0+H} \rho_{i,t}^{(d)} \cdot C_{i,t}^{(d)}
$$

where `(d)` indexes posterior draws and `H` is the planned follow-up
horizon (52 or 104 weeks). Site priority is the lower 10% quantile of the
posterior of `Accrual_i`, paired with:

- Design-specific power calculation (test-negative design and
  cluster-randomised options pre-specified);
- Operational feasibility filter (LACEN coverage, primary-care
  infrastructure, ethics-approval lead time);
- Robustness check: re-run ranking under each sensitivity in Section 7.3
  and report Spearman rank stability.

The deliverable is not a single ranked list but a **robust top-tier set**:
municipalities that remain in the top quartile across all sensitivity
scenarios.

*(react: clarify how this connects to actual IXCHIQ Phase 4 protocol?)*

---

## 6. Transferability — Oropouche demonstration

We will refit Layer 1 to the 2023–2024 Oropouche outbreak in Brazil with
disease-specific adaptations:

- Generation interval re-parameterised (Oropouche mean ~4–7 days);
- Vector replaced (*Culicoides paraensis* rather than *Aedes*);
- Reporting fraction informed by the post-2023 testing surge;
- Layer 2 not applied (no Oropouche vaccine currently licensed).

This demonstrates that the reconstruction framework is *disease-agnostic*
and supports the rapid-response value claim for emerging arboviruses.

*(react: should we also pilot dengue or Zika as an additional case
study to broaden the transferability claim?)*

---

## 7. Validation strategy

### 7.1 Hold-out predictive validation

Training period: 2014–2022. Test period: 2023–2024.
Primary metrics:

- Spearman rank correlation between predicted and observed top-50
  municipalities by lab-confirmed cases.
- Top-10 municipality overlap.
- Calibration of 80% predictive intervals (empirical coverage).
- Continuous Ranked Probability Score (CRPS) per municipality.

### 7.2 Cross-implementation consistency

Stan vs INLA on a shared state-level subset. Agreement assessed by
Wasserstein-1 distance between marginal posteriors and by Spearman
correlation of point estimates.

### 7.3 Sensitivity analyses — *reported as main results, not supplementary*

Five pre-specified sensitivities:

1. Generation interval (mean ±50%)
2. Initial seroprevalence (with vs without serology anchoring)
3. Vaccine efficacy uncertainty (Beta priors widened to reflect
   correlate-of-protection licensure)
4. Climate scenario (RCP 4.5, RCP 8.5, observational continuation)
5. Spatial resolution (state, health region, municipality)

For each sensitivity, both site ranking and strategy ranking are
re-evaluated; **robustness of ranking under uncertainty is the primary
deliverable**, distinct from point estimates.

---

## 8. Software, open science, reproducibility

- All code released under MIT license at
  github.com/{owner}/CHIK_benefit_risk_modelling
- Two R packages developed and CRAN-submitted:
  - `chikBR.layer1` — Bayesian reconstruction engine, designed for transfer
    to other arboviruses
  - `chikBR.layer2` — vaccination forward simulator with documented
    interface to Layer 1 outputs
- Workflow reproducible via `targets` pipeline
- All posterior draws archived on Zenodo with DOI
- Pre-registration of validation protocol on OSF before model fitting

---

## 9. Limitations to declare upfront

*Drafted explicitly so reviewers do not have to discover them.*

1. **ρ–C identifiability remains data-limited.** Serology anchoring
   resolves this at *anchored municipalities only*; for others, ρ and C are
   only jointly determined via spatial smoothing assumptions.
2. **12-week+ forecasts of absolute incidence remain hard** for
   chikungunya. The framework's forecast products are explicitly
   *probabilistic / ranking-oriented*, not point-prediction; this is
   stated in deliverables.
3. **Vector compartment is absorbed into FOI.** This is sufficient for
   vaccination strategy comparison but precludes direct evaluation of
   vector control interventions. Such evaluation is out of scope.
4. **Vaccine efficacy assumptions are inherited** from
   correlate-of-protection licensure rather than measured. Sensitivity
   analyses span the plausible range.
5. **Brazil-specific data assets enable this work.** Direct transfer to
   countries without SINAN-equivalent surveillance + serosurvey access
   will require additional development.

*(react: too many? merge with discussion? frame as "design choices"
rather than limitations?)*

---

## 10. Open questions for next iteration

- [ ] Does the title carry the four deliverables clearly?
- [ ] WP structure: 4 WPs (one per deliverable) or 3 WPs (combining
      reconstruction + projection)?
- [ ] Should Layer 2 be split into "vaccination strategy" vs "VE site
      selection" as separate WPs or kept together?
- [ ] Timeline: 3-year project? 5-year?
- [ ] Co-investigators / collaboration anchors to name (FioCruz,
      Cauchemez group, Lowe / Bhatt Imperial?)
- [ ] Specific funder format: MRC / Wellcome / CEPI / Gavi /
      JST? Each has different page limits and structure.
- [ ] Patient & public involvement — required by which funders?
- [ ] Data sharing / IRB considerations specific to SINAN access.

---

## 11. References to cite (working list — verify formatting per funder)

1. Abbott S, *et al.* (2020). Estimating the time-varying reproduction
   number of SARS-CoV-2 using national and subnational case counts.
   *Wellcome Open Research* 5:112.
2. Flaxman S, *et al.* (2020). Estimating the effects of
   non-pharmaceutical interventions on COVID-19 in Europe.
   *Nature* 584:257–261.
3. Watson OJ, *et al.* (2022). Global impact of the first year of
   COVID-19 vaccination: a mathematical modelling study.
   *Lancet Infect Dis* 22:1293–1302.
4. Bhatt S, *et al.* (2013). The global distribution and burden of
   dengue. *Nature* 496:504–507.
5. Lowe R, *et al.* (2013). The development of an early warning system
   for climate-sensitive disease risk with a focus on dengue epidemics
   in Southeast Brazil. *Stat Med* 32:864–883.
6. Lowe R, *et al.* (2016). Evaluating probabilistic dengue risk
   forecasts from a prototype early warning system for Brazil.
   *PLoS Med* 13:e1002000.
7. Lowe R, *et al.* (2021). Combined effects of hydrometeorological
   hazards and urbanisation on dengue risk in Brazil. *Lancet Planet
   Health* 5:e209–e219.
8. Ribeiro Dos Santos G, *et al.* (2025). Global burden of chikungunya
   virus infections and the potential benefit of vaccination campaigns.
   *Nat Med* (online).
9. Mordecai EA, *et al.* (2017). Detecting the impact of temperature on
   transmission of Zika, dengue, and chikungunya using mechanistic
   models. *PLoS Negl Trop Dis* 11:e0005568.
10. Liu-Helmersson J, *et al.* (2014). Vectorial capacity of *Aedes
    aegypti*: effects of temperature and implications for global
    dengue epidemic potential. *PLoS ONE* 9:e89783.
11. Kraemer MUG, *et al.* (2019). Past and future spread of the
    arbovirus vectors *Aedes aegypti* and *Aedes albopictus*. *Nature
    Microbiology* 4:854–863.
12. Prem K, *et al.* (2017). Projecting social contact matrices in 152
    countries using contact surveys and demographic data. *PLoS Comput
    Biol* 13:e1005697.
13. Funk S, *et al.* (2019). Assessing the performance of real-time
    epidemic forecasts: A case study of Ebola in the Western Area
    region of Sierra Leone, 2014–15. *PLoS Comput Biol* 15:e1006785.
14. Brito AF, Cunha MS, Dias L, Donalisio M [serology survey series
    — fill exact citations].

---

## 12. Next 7 days — concrete writing milestones

Day 1 (today): mark up v0 with red pen. Decide WP structure.
Day 2: rewrite Section 1 + 2 to v1 quality.
Day 3: rewrite Section 3 + finalise Layer 1 equations.
Day 4: rewrite Section 4 + finalise vaccination strategy table.
Day 5: rewrite Section 5–6 (VE site + Oropouche transfer).
Day 6: write Section 7 (validation) properly; this is the section that
       discriminates a strong proposal from a weak one.
Day 7: complete pass for grant-format-specific tightening (page limits,
       headers, figure plan).

By end of Day 7: v1 ready to share with one collaborator for honest
feedback. Do not perfect v0 → v1 alone; share at v1 stage.

---

*end of v0 draft*
