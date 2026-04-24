# =============================================================================
# setup_age_props.R
#
# Age-stratified chronic-stage proportions (acute / subac / chr6m / chr12m /
# chr30m). Extends lhs_sample built in setup.R.
#
# Convention
#   - u40 (<40 y): uses the existing Beta parameters already in setup.R
#     (line 689-695: acute, subac, chr6m, chr12m, chr30m). Columns are aliased
#     as *_u40 so the lookup logic is explicit.
#   - o40 (>=40 y, applied exclusively to 65+): new Beta parameters supplied
#     below, sampled with a fresh LHS block.
#
# Age group -> proportion set used
#   1-11   : props_u40                       (all <40)
#   12-17  : props_u40                       (all <40)
#   18-64  : props_1864 = w_u40*u40 + w_o40*o40  (weighted blend)
#   65+    : props_o40                       (all >=40)
#
# Prereq: setup.R must have been sourced so that `lhs_sample` (data.frame)
# exists with columns "acute", "subac", "chr6m", "chr12m", "chr30m",
# and `compute_daly_one()` is in the global env.
# =============================================================================

stopifnot(
  exists("lhs_sample", envir = globalenv()),
  all(c("acute", "subac", "chr6m", "chr12m", "chr30m")
      %in% colnames(lhs_sample)),
  exists("compute_daly_one", envir = globalenv())
)

if (!requireNamespace("lhs", quietly = TRUE)) {
  stop("Package 'lhs' is required (used by setup.R already).")
}

# -----------------------------------------------------------------------------
# 1) Population weights within 18-64 (w_u40 + w_o40 == 1)
#    Two ways to set these:
#      (a) Hard-code from an external census (see IBGE 2022 default below).
#      (b) Compute from an in-project population table via
#          compute_w_u40_from_age_pop().
# -----------------------------------------------------------------------------
# (a) IBGE 2022 (Brazil, national) rough default:
#     18-39: ~92.6 M,  40-64: ~70.5 M  ->  w_u40 = 92.6 / (92.6 + 70.5)
w_u40_ibge2022 <- 92.6 / (92.6 + 70.5)  # ~= 0.568
w_o40_ibge2022 <- 1 - w_u40_ibge2022

# (b) Helper: compute weights from an age-population data.frame.
#     Expects two columns: one numeric/integer age (representative age of the
#     bin, e.g. midpoint or lower edge) and one population count. Rows outside
#     the 18-64 range are ignored.
#
#     Example (using the 5-year bins from brazil_all_draws_ori_v2.R):
#       age_pop_df <- data.frame(
#         age_mid = c(18.5, 22, 27, 32, 37, 42, 47, 52, 57, 62),  # AgeGroup 6:15
#         N       = c(...)                                          # pop in each bin
#       )
#       w <- compute_w_u40_from_age_pop(age_pop_df, "age_mid", "N", cutoff = 40)
#       w_u40 <- w$w_u40; w_o40 <- w$w_o40
compute_w_u40_from_age_pop <- function(age_pop_df,
                                       age_col = "age",
                                       pop_col = "pop",
                                       cutoff  = 40,
                                       age_min = 18,
                                       age_max = 64) {
  stopifnot(
    is.data.frame(age_pop_df),
    age_col %in% names(age_pop_df),
    pop_col %in% names(age_pop_df)
  )
  a <- as.numeric(age_pop_df[[age_col]])
  n <- as.numeric(age_pop_df[[pop_col]])
  keep <- a >= age_min & a <= age_max & is.finite(a) & is.finite(n)
  a <- a[keep]; n <- n[keep]
  if (sum(n) <= 0) stop("Total population in 18-64 range is zero.")
  w_u <- sum(n[a <  cutoff]) / sum(n)
  w_o <- sum(n[a >= cutoff]) / sum(n)
  list(w_u40 = w_u, w_o40 = w_o,
       pop_u40 = sum(n[a <  cutoff]),
       pop_o40 = sum(n[a >= cutoff]),
       cutoff = cutoff)
}

# Active weights used by props_1864 and the downstream wrapper.
# To override: set w_u40 / w_o40 BEFORE sourcing this file, or edit this line.
if (!exists("w_u40", envir = globalenv(), inherits = FALSE)) {
  w_u40 <- w_u40_ibge2022
}
w_o40 <- 1 - w_u40
stopifnot(abs((w_u40 + w_o40) - 1) < 1e-12)

# -----------------------------------------------------------------------------
# 2) o40 (>=40 y) Beta parameters — applied to 65+
#    (User-supplied values)
# -----------------------------------------------------------------------------
o40_beta_params <- list(
  acute  = c(shape1 = 388.2874, shape2 = 742.4988),
  subac  = c(shape1 = 2487.085, shape2 = 6450.815),
  chr6m  = c(shape1 = 5998.872, shape2 = 21654.86),
  chr12m = c(shape1 = 184.6674, shape2 = 1216.058),
  chr30m = c(shape1 = 15.74654, shape2 = 516.3419)
)

# -----------------------------------------------------------------------------
# 3) Add *_u40 and *_o40 columns to lhs_sample
# -----------------------------------------------------------------------------
# *_u40: alias the existing <40 samples (setup.R line 689-695).
lhs_sample$acute_u40  <- lhs_sample$acute
lhs_sample$subac_u40  <- lhs_sample$subac
lhs_sample$chr6m_u40  <- lhs_sample$chr6m
lhs_sample$chr12m_u40 <- lhs_sample$chr12m
lhs_sample$chr30m_u40 <- lhs_sample$chr30m

# *_o40: fresh LHS sampling with user-provided Beta parameters.
set.seed(20260420)
A_o40 <- lhs::randomLHS(n = nrow(lhs_sample), k = 5)

lhs_sample$acute_o40  <- qbeta(
  A_o40[, 1],
  shape1 = o40_beta_params$acute["shape1"],
  shape2 = o40_beta_params$acute["shape2"]
)
lhs_sample$subac_o40  <- qbeta(
  A_o40[, 2],
  shape1 = o40_beta_params$subac["shape1"],
  shape2 = o40_beta_params$subac["shape2"]
)
lhs_sample$chr6m_o40  <- qbeta(
  A_o40[, 3],
  shape1 = o40_beta_params$chr6m["shape1"],
  shape2 = o40_beta_params$chr6m["shape2"]
)
lhs_sample$chr12m_o40 <- qbeta(
  A_o40[, 4],
  shape1 = o40_beta_params$chr12m["shape1"],
  shape2 = o40_beta_params$chr12m["shape2"]
)
lhs_sample$chr30m_o40 <- qbeta(
  A_o40[, 5],
  shape1 = o40_beta_params$chr30m["shape1"],
  shape2 = o40_beta_params$chr30m["shape2"]
)

# -----------------------------------------------------------------------------
# 4) Props lists (one element per draw; length = nrow(lhs_sample))
# -----------------------------------------------------------------------------
props_u40 <- list(
  acute  = lhs_sample$acute_u40,
  subac  = lhs_sample$subac_u40,
  chr6m  = lhs_sample$chr6m_u40,
  chr12m = lhs_sample$chr12m_u40,
  chr30m = lhs_sample$chr30m_u40
)

props_o40 <- list(
  acute  = lhs_sample$acute_o40,
  subac  = lhs_sample$subac_o40,
  chr6m  = lhs_sample$chr6m_o40,
  chr12m = lhs_sample$chr12m_o40,
  chr30m = lhs_sample$chr30m_o40
)

props_1864 <- list(
  acute  = lhs_sample$acute_u40  * w_u40 + lhs_sample$acute_o40  * w_o40,
  subac  = lhs_sample$subac_u40  * w_u40 + lhs_sample$subac_o40  * w_o40,
  chr6m  = lhs_sample$chr6m_u40  * w_u40 + lhs_sample$chr6m_o40  * w_o40,
  chr12m = lhs_sample$chr12m_u40 * w_u40 + lhs_sample$chr12m_o40 * w_o40,
  chr30m = lhs_sample$chr30m_u40 * w_u40 + lhs_sample$chr30m_o40 * w_o40
)

# -----------------------------------------------------------------------------
# 5) Helpers
# -----------------------------------------------------------------------------
# Scalar lookup: given an age category and a single draw_id, return the 5
# proportion values (acute / subac / chr6m / chr12m / chr30m) to use.
get_props_for_agecat_draw <- function(age_cat, draw_id,
                                      w_u40_local = w_u40,
                                      w_o40_local = w_o40) {
  stopifnot(
    length(age_cat) == 1,
    length(draw_id) == 1,
    draw_id >= 1, draw_id <= nrow(lhs_sample)
  )
  pull5 <- function(row_i, suffix) {
    c(
      acute  = lhs_sample[[paste0("acute",  suffix)]][row_i],
      subac  = lhs_sample[[paste0("subac",  suffix)]][row_i],
      chr6m  = lhs_sample[[paste0("chr6m",  suffix)]][row_i],
      chr12m = lhs_sample[[paste0("chr12m", suffix)]][row_i],
      chr30m = lhs_sample[[paste0("chr30m", suffix)]][row_i]
    )
  }
  switch(
    as.character(age_cat),
    "1-11"  = pull5(draw_id, "_u40"),
    "12-17" = pull5(draw_id, "_u40"),
    "65+"   = pull5(draw_id, "_o40"),
    "18-64" = {
      u <- pull5(draw_id, "_u40")
      o <- pull5(draw_id, "_o40")
      u * w_u40_local + o * w_o40_local
    },
    stop("Unknown AgeCat: ", age_cat)
  )
}

# Vector lookup: age-specific proportions for an entire draw vector.
# Returns a data.frame with one row per draw, columns acute/subac/chr6m/
# chr12m/chr30m, already age-resolved.
get_props_df_for_agecat <- function(age_cat,
                                    w_u40_local = w_u40,
                                    w_o40_local = w_o40) {
  switch(
    as.character(age_cat),
    "1-11"  = data.frame(
      acute  = lhs_sample$acute_u40,
      subac  = lhs_sample$subac_u40,
      chr6m  = lhs_sample$chr6m_u40,
      chr12m = lhs_sample$chr12m_u40,
      chr30m = lhs_sample$chr30m_u40
    ),
    "12-17" = data.frame(
      acute  = lhs_sample$acute_u40,
      subac  = lhs_sample$subac_u40,
      chr6m  = lhs_sample$chr6m_u40,
      chr12m = lhs_sample$chr12m_u40,
      chr30m = lhs_sample$chr30m_u40
    ),
    "65+"   = data.frame(
      acute  = lhs_sample$acute_o40,
      subac  = lhs_sample$subac_o40,
      chr6m  = lhs_sample$chr6m_o40,
      chr12m = lhs_sample$chr12m_o40,
      chr30m = lhs_sample$chr30m_o40
    ),
    "18-64" = data.frame(
      acute  = lhs_sample$acute_u40  * w_u40_local + lhs_sample$acute_o40  * w_o40_local,
      subac  = lhs_sample$subac_u40  * w_u40_local + lhs_sample$subac_o40  * w_o40_local,
      chr6m  = lhs_sample$chr6m_u40  * w_u40_local + lhs_sample$chr6m_o40  * w_o40_local,
      chr12m = lhs_sample$chr12m_u40 * w_u40_local + lhs_sample$chr12m_o40 * w_o40_local,
      chr30m = lhs_sample$chr30m_u40 * w_u40_local + lhs_sample$chr30m_o40 * w_o40_local
    ),
    stop("Unknown AgeCat: ", age_cat)
  )
}

# -----------------------------------------------------------------------------
# 6) Drop-in wrapper around compute_daly_one() that injects age-specific props
# -----------------------------------------------------------------------------
# Usage mirrors compute_daly_one(). Pass draw_pars with the usual DW/duration
# parameters AND the 10 age-specific proportion columns:
#   acute_u40, acute_o40, subac_u40, subac_o40, chr6m_u40, chr6m_o40,
#   chr12m_u40, chr12m_o40, chr30m_u40, chr30m_o40.
# If draw_pars already contains unsuffixed "acute"/"subac"/... (e.g., a
# record constructed from the full lhs_sample row), the suffixed values are
# still used for age-aware blending.
compute_daly_one_age_specific <- function(age_group,
                                          deaths_10k       = 0,
                                          hosp_10k         = 0,
                                          nonhosp_symp_10k = 0,
                                          symp_10k         = NULL,
                                          sae_10k          = 0,
                                          deaths_sae_10k   = 0,
                                          draw_pars,
                                          draw_id          = NULL,
                                          w_u40_local = w_u40,
                                          w_o40_local = w_o40) {

  needed <- c("acute_u40", "acute_o40",
              "subac_u40", "subac_o40",
              "chr6m_u40", "chr6m_o40",
              "chr12m_u40","chr12m_o40",
              "chr30m_u40","chr30m_o40")

  missing_cols <- setdiff(needed, names(draw_pars))
  if (length(missing_cols) > 0) {
    # Fallback: pull them from global lhs_sample using draw_id.
    if (is.null(draw_id)) {
      stop(
        "compute_daly_one_age_specific(): draw_pars is missing age-specific ",
        "columns (", paste(missing_cols, collapse = ", "), ") ",
        "and no draw_id was supplied for fallback lookup. Either add these ",
        "columns to draw_pars or pass draw_id = <id>."
      )
    }
    if (!exists("lhs_sample", envir = globalenv())) {
      stop("Global lhs_sample not found; cannot fall back by draw_id.")
    }
    ls_ref <- get("lhs_sample", envir = globalenv())
    if (!all(needed %in% colnames(ls_ref))) {
      stop("Global lhs_sample does not have *_u40/*_o40 columns. ",
           "Did you source('02_Scripts/setup_age_props.R') after setup.R?")
    }
    for (col in missing_cols) {
      draw_pars[[col]] <- ls_ref[[col]][draw_id]
    }
  }

  ac_u  <- draw_pars$acute_u40
  ac_o  <- draw_pars$acute_o40
  sb_u  <- draw_pars$subac_u40
  sb_o  <- draw_pars$subac_o40
  c6_u  <- draw_pars$chr6m_u40
  c6_o  <- draw_pars$chr6m_o40
  c12_u <- draw_pars$chr12m_u40
  c12_o <- draw_pars$chr12m_o40
  c30_u <- draw_pars$chr30m_u40
  c30_o <- draw_pars$chr30m_o40

  blend <- function(u, o) u * w_u40_local + o * w_o40_local

  eff <- switch(
    as.character(age_group),
    "1-11"  = list(acute = ac_u,            subac = sb_u,            chr6m = c6_u,            chr12m = c12_u,             chr30m = c30_u),
    "12-17" = list(acute = ac_u,            subac = sb_u,            chr6m = c6_u,            chr12m = c12_u,             chr30m = c30_u),
    "18-64" = list(acute = blend(ac_u, ac_o), subac = blend(sb_u, sb_o), chr6m = blend(c6_u, c6_o), chr12m = blend(c12_u, c12_o), chr30m = blend(c30_u, c30_o)),
    "65+"   = list(acute = ac_o,            subac = sb_o,            chr6m = c6_o,            chr12m = c12_o,             chr30m = c30_o),
    stop("Unknown age_group: ", age_group)
  )

  draw_pars$acute  <- eff$acute
  draw_pars$subac  <- eff$subac
  draw_pars$chr6m  <- eff$chr6m
  draw_pars$chr12m <- eff$chr12m
  draw_pars$chr30m <- eff$chr30m

  compute_daly_one(
    age_group        = age_group,
    deaths_10k       = deaths_10k,
    hosp_10k         = hosp_10k,
    nonhosp_symp_10k = nonhosp_symp_10k,
    symp_10k         = symp_10k,
    sae_10k          = sae_10k,
    deaths_sae_10k   = deaths_sae_10k,
    draw_pars        = draw_pars
  )
}

# -----------------------------------------------------------------------------
# 7) Sanity summary
# -----------------------------------------------------------------------------
message(
  "[setup_age_props] Added columns to lhs_sample: ",
  paste(c("acute_u40", "subac_u40", "chr6m_u40", "chr12m_u40", "chr30m_u40",
          "acute_o40", "subac_o40", "chr6m_o40", "chr12m_o40", "chr30m_o40"),
        collapse = ", "),
  "\n  w_u40 = ", w_u40, ", w_o40 = ", w_o40,
  "\n  props_u40 / props_1864 / props_o40 available (length ",
  nrow(lhs_sample), " each)",
  "\n  Helpers: get_props_for_agecat_draw(), get_props_df_for_agecat(), ",
  "compute_daly_one_age_specific()"
)

invisible(list(
  o40_beta_params = o40_beta_params,
  w_u40 = w_u40, w_o40 = w_o40,
  o40_summary = data.frame(
    param = c("acute", "subac", "chr6m", "chr12m", "chr30m"),
    shape1 = vapply(o40_beta_params, `[[`, numeric(1), "shape1"),
    shape2 = vapply(o40_beta_params, `[[`, numeric(1), "shape2"),
    mean_o40  = vapply(
      o40_beta_params, function(p) unname(p["shape1"] / sum(p)), numeric(1)
    )
  )
))
