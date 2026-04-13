# =============================================================================
# SECTION 00. Script map (debugging index)
# =============================================================================
# 00) Script map and section index
# 01) Global setup and shared inputs
# 02) Build baseline draw objects (symptomatic burden)
# 03) Shared helper functions (rho scaling)
# 04) Hospitalisation draws
# 05) Fatal draws
# 06) DALY draws
# 07) SAE draws and persistence of core draw objects
# 08) BRR pipeline (draw-level benefit vs risk; no object mixing)
#
# NOTE:
# Functions are intentionally defined before their first use throughout the
# script. For R scripts sourced top-to-bottom, this prevents "object not found"
# issues during execution and makes debugging easier.
# =============================================================================

# =============================================================================
# SECTION 01. Global setup and shared inputs
# =============================================================================

age_groups <- c(mean(0:1),
                mean(1:4),
                mean(5:9),
                mean(10:11),
                mean(12:17),
                mean(18:19),
                mean(20:24),
                mean(25:29),
                mean(30:34),
                mean(35:39),
                mean(40:44),
                mean(45:49),
                mean(50:54),
                mean(55:59),
                mean(60:64),
                mean(65:69),
                mean(70:74),
                mean(75:79),
                mean(80:84),
                mean(85:89)
)

conflicted::conflicts_prefer(dplyr::filter)

# =============================================================================
# SECTION 02. Build baseline draw objects (symptomatic burden)
# =============================================================================
# Step 2.1) Collect posterior rho draws by region
posterior_list <- list(
  ce = posterior_ce,
  bh = posterior_bh,
  pa = posterior_pa,
  pn = posterior_pn,
  rg = posterior_rg,
  pi = posterior_pi,
  ag = posterior_ag,
  tc = posterior_tc,
  mg = posterior_mg,
  se = posterior_se,
  go = posterior_go
)

rho_pool <- purrr::imap_dfr(posterior_list, function(post, region_name){
  tibble(
    Region = region_name,
    rho_id = seq_along(post$rho),
    rho    = as.numeric(post$rho)
  )
})

save(posterior_list, file = "01_Data/posterior_list.RData")

preui_all <- list(
  preui_ce,
  preui_ag,
  preui_bh,
  preui_go,
  preui_mg,
  preui_pa,
  preui_pi,
  preui_pn,
  preui_rg,
  preui_se,
  preui_tc
)

save(preui_all, file = "01_Data/preui_all.RData")


names(preui_all) <- c("Ceará","Alagoas", "Bahia", "Goiás",
                      "Minas Gerais", "Paraíba", "Piauí",
                      "Pernambuco", "Rio Grande do Norte", 
                      "Sergipe", "Tocantins")
setting_key <- c(
  "Ceará"             = "High",
  "Bahia"             = "Low",
  "Paraíba"           = "High",
  "Pernambuco"        = "Moderate",
  "Rio Grande do Norte" = "Low",
  "Piauí"             = "High",
  "Tocantins"         = "Moderate",
  "Alagoas"           = "High",
  "Minas Gerais"      = "Low",
  "Sergipe"           = "Low",
  "Goiás"             = "Low"
)

# Step 2.2) Utility to aggregate draw-level pre/post burden totals
calc_total_impact_draws <- function(pre_list, post_arr) {
  pre_arr <- simplify2array(pre_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

# Step 2.3) Build draw-level symptomatic burden (unscaled)
all_draws_ix <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        draw_totals <- calc_total_impact_draws(pre_list, post_arr)
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})


# =============================================================================
# SECTION 03. Shared helper functions (rho scaling)
# =============================================================================
# Step 3.1) Scale symptomatic burden by rho to recover true burden
set.seed(1)
scale_symp_by_rho_pre <- function(pre_list, rho_vec) {
  # pre_list: length n_draws, each is (age x week) matrix
  # rho_vec: length n_draws
  Map(function(mat, r) mat / r, pre_list, rho_vec)
}

scale_symp_by_rho_post <- function(post_arr, rho_vec) {
  # post_arr: (age x week x draw)
  # rho_vec: length draw
  for (i in seq_along(rho_vec)) post_arr[,,i] <- post_arr[,,i] / rho_vec[i]
  post_arr
}


# Step 3.2) Apply rho scaling to symptomatic burden draws
all_draws_ix_true <- all_draws_ix %>%
  group_by(Region) %>%
  mutate(
    rho = sample(rho_pool$rho[rho_pool$Region == first(Region)],
                 size = n(), replace = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    total_pre_true  = total_pre  / rho,
    total_post_true = total_post / rho
  )

# =============================================================================
# SECTION 04. Hospitalisation draws
# =============================================================================
# Step 4.1) Convert symptomatic draws to hospitalisation draws
make_hosp_draws <- function(symp_list, hosp_rate) {
  lapply(symp_list, function(symp_matrix) {
    sweep(symp_matrix, 1, hosp_rate, `*`)
  })
}

calc_total_hosp_draws <- function(pre_list, post_arr, hosp_rate) {
  # Pre
  pre_hosp_list <- make_hosp_draws(pre_list, hosp_rate)
  pre_arr <- simplify2array(pre_hosp_list)
  if (length(dim(pre_arr)) == 2) pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  
  # Post
  post_list <- lapply(seq_len(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_hosp_list <- make_hosp_draws(post_list, hosp_rate)
  post_arr2 <- simplify2array(post_hosp_list)
  if (length(dim(post_arr2)) == 2) post_arr2 <- array(post_arr2, dim = c(dim(post_arr2), 1))
  
  # Align
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  pre_arr   <- pre_arr[ , , 1:n_draws, drop = FALSE]
  post_arr2 <- post_arr2[ , , 1:n_draws, drop = FALSE]
  
  # Sum
  tibble(
    draw_id    = seq_len(n_draws),
    total_pre  = apply(pre_arr,  3, sum, na.rm = TRUE),
    total_post = apply(post_arr2, 3, sum, na.rm = TRUE)
  )
}

calc_total_hosp_draws_rho <- function(pre_list, post_arr, hosp_rate, rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws, drop = FALSE]
  rho_vec  <- rho_vec[1:n_draws]
  
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  calc_total_hosp_draws(pre_list_true, post_arr_true, hosp_rate = hosp_rate)
}

set.seed(1)

# Step 4.2) Build rho-adjusted hospitalisation draw table
all_draws_hosp_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        n_draws <- min(length(pre_list), dim(post_arr)[3])

        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_hosp_draws_rho(
          pre_list = pre_list,
          post_arr = post_arr,
          hosp_rate = hosp,
          rho_vec = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})

# =============================================================================
# SECTION 05. Fatal draws
# =============================================================================
# Step 5.1) Split symptomatic burden into hospitalised/non-hospitalised fatal burden
make_fatal_hosp_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate) {
  lapply(symp_list, function(symp_matrix) {
    hospitalised     <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised <- symp_matrix - hospitalised
    
    fatal <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    
    list(
      hosp  = hospitalised,
      fatal = fatal
    )
  })
}

# NOTE:
# This function calculates the hospitalisation component using the same
# fatal-model decomposition pipeline. It is intentionally named differently
# from calc_total_hosp_draws(...) above to avoid accidental function override.
calc_total_hosp_from_fatal_model <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_hosp_list <- lapply(pre_conv, `[[`, "hosp")
  pre_arr <- simplify2array(pre_hosp_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq_len(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_hosp_list <- lapply(post_conv, `[[`, "hosp")
  post_arr2 <- simplify2array(post_hosp_list)
  if (length(dim(post_arr2)) == 2) {
    post_arr2 <- array(post_arr2, dim = c(dim(post_arr2), 1))
  }
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  pre_arr   <- pre_arr[ , , 1:n_draws, drop = FALSE]
  post_arr2 <- post_arr2[ , , 1:n_draws, drop = FALSE]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr2, 3, sum, na.rm = TRUE)
  
  tibble(
    draw_id    = seq_len(n_draws),
    total_pre  = total_pre,
    total_post = total_post
  )
}

# Reuse shared rho-scaling helpers defined in SECTION 03.

calc_total_fatal_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_fatal_list <- lapply(pre_conv, `[[`, "fatal")
  pre_arr <- simplify2array(pre_fatal_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_fatal_list <- lapply(post_conv, `[[`, "fatal")
  post_arr <- simplify2array(post_fatal_list)
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

calc_total_fatal_draws_rho <- function(pre_list, post_arr,
                                       hosp_rate, fatal_rate, nh_fatal_rate,
                                       rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  calc_total_fatal_draws(pre_list_true, post_arr_true,
                         hosp_rate, fatal_rate, nh_fatal_rate)
}


set.seed(1)

# Step 5.2) Build rho-adjusted fatal draw table
all_draws_fatal_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # Region-specific rho pool (required)
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_fatal_draws_rho(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          rho_vec       = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})

# =============================================================================
# SECTION 06. DALY draws
# =============================================================================
# Step 6.1) Convert symptomatic burden to DALY components
make_daly_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate,
                            dw_hosp, dw_nonhosp, dw_chronic,
                            dur_acute, dur_subacute, dur_chronic,
                            dw_subacute, subac_prop, chr_prop,
                            le_left_vec) {
  lapply(symp_list, function(symp_matrix) {
    hospitalised       <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised   <- symp_matrix - hospitalised
    fatal              <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    
    yld_acute    <- (hospitalised * dw_hosp     * dur_acute) +
      (non_hospitalised * dw_nonhosp * dur_acute)
    yld_subacute <- (hospitalised * subac_prop * dw_subacute * dur_subacute) +
      (non_hospitalised * chr_prop * dw_subacute * dur_subacute)
    yld_chronic  <- (hospitalised * chr_prop * dw_chronic * dur_chronic) +
      (non_hospitalised * chr_prop * dw_chronic * dur_chronic)
    
    yld_total <- yld_acute + yld_subacute + yld_chronic
    
    # Precomputed life expectancy by age group
    yll <- sweep(fatal, 1, le_left_vec, `*`)
    
    daly_tot <- yld_total + yll
    
    list(
      hosp     = hospitalised,
      fatal    = fatal,
      yld_tot  = yld_total,
      yll      = yll,
      daly_tot = daly_tot
    )
  })
}
le_by_age <- function(age_numeric) {
  if (age_numeric <= 1) return(quantile(le_sample$le_1, 0.5))
  if (age_numeric < 20) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 30) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 40) return(quantile(le_sample$le_3, 0.5))
  if (age_numeric < 50) return(quantile(le_sample$le_4, 0.5))
  if (age_numeric < 60) return(quantile(le_sample$le_5, 0.5))
  if (age_numeric < 70) return(quantile(le_sample$le_6, 0.5))
  if (age_numeric < 80) return(quantile(le_sample$le_7, 0.5))
  return(quantile(le_sample$le_8, 0.5))
}

calc_total_daly_draws <- function(pre_list, post_arr,
                                  hosp_rate, fatal_rate, nh_fatal_rate,
                                  dw_hosp, dw_nonhosp, dw_chronic,
                                  dur_acute, dur_subacute, dur_chronic,
                                  dw_subacute, subac_prop, chr_prop,
                                  le_left_vec) {
  # Pre
  pre_conv <- make_daly_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate,
                              dw_hosp, dw_nonhosp, dw_chronic,
                              dur_acute, dur_subacute, dur_chronic,
                              dw_subacute, subac_prop, chr_prop,
                              le_left_vec)
  pre_daly_list <- lapply(pre_conv, `[[`, "daly_tot")
  pre_arr <- simplify2array(pre_daly_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # Post
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_daly_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate,
                               dw_hosp, dw_nonhosp, dw_chronic,
                               dur_acute, dur_subacute, dur_chronic,
                               dw_subacute, subac_prop, chr_prop,
                               le_left_vec)
  post_daly_list <- lapply(post_conv, `[[`, "daly_tot")
  post_arr <- simplify2array(post_daly_list)
  
  # Align
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

le_left_vec <- sapply(age_groups, le_by_age)

# Reuse shared rho-scaling helpers defined in SECTION 03.

calc_total_daly_draws_rho <- function(pre_list, post_arr,
                                      hosp_rate, fatal_rate, nh_fatal_rate,
                                      dw_hosp, dw_nonhosp, dw_chronic,
                                      dur_acute, dur_subacute, dur_chronic,
                                      dw_subacute, subac_prop, chr_prop,
                                      le_left_vec,
                                      rho_vec) {
  
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  calc_total_daly_draws(
    pre_list = pre_list_true,
    post_arr = post_arr_true,
    hosp_rate     = hosp_rate,
    fatal_rate    = fatal_rate,
    nh_fatal_rate = nh_fatal_rate,
    dw_hosp       = dw_hosp,
    dw_nonhosp    = dw_nonhosp,
    dw_chronic    = dw_chronic,
    dur_acute     = dur_acute,
    dur_subacute  = dur_subacute,
    dur_chronic   = dur_chronic,
    dw_subacute   = dw_subacute,
    subac_prop    = subac_prop,
    chr_prop      = chr_prop,
    le_left_vec   = le_left_vec
  )
}


set.seed(1)

# Step 6.2) Build rho-adjusted DALY draw table
all_draws_daly_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_daly_draws_rho(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          
          dw_hosp       = quantile(lhs_sample_young$dw_hosp, 0.5),
          dw_nonhosp    = quantile(lhs_sample_young$dw_nonhosp, 0.5),
          dw_chronic    = quantile(lhs_sample_young$dw_chronic, 0.5),
          
          dur_acute     = quantile(lhs_sample_young$dur_acute, 0.5),
          dur_subacute  = quantile(lhs_sample_young$dur_subac, 0.5),
          dur_chronic   = quantile(lhs_sample_young$dur_chronic, 0.5),
          
          dw_subacute   = quantile(lhs_sample_young$dw_subac, 0.5),
          subac_prop    = quantile(lhs_sample_young$subac, 0.5),
          chr_prop      = (quantile(lhs_sample_young$chr6m, 0.5) +
                             quantile(lhs_sample_young$chr12m, 0.5) +
                             quantile(lhs_sample_young$chr30m, 0.5)),
          
          le_left_vec   = le_left_vec,
          rho_vec       = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})

# =============================================================================
# SECTION 07. SAE draws and persistence of core draw objects
# =============================================================================
# Step 7.1) SAE = hospitalisation + fatal events
all_draws_sae_true <- all_draws_hosp_true %>%
  dplyr::rename(
    total_pre_hosp  = total_pre,
    total_post_hosp = total_post
  ) %>%
  inner_join(
    all_draws_fatal_true %>%
      dplyr::rename(
        total_pre_fatal  = total_pre,
        total_post_fatal = total_post
      ),
    by = c("draw_id", "Region", "VE", "Coverage", "Scenario")
  ) %>%
  mutate(
    total_pre  = total_pre_hosp  + total_pre_fatal,
    total_post = total_post_hosp + total_post_fatal
  ) %>%
  dplyr::select(draw_id, total_pre, total_post, Region, VE, Coverage, Scenario)



# Step 7.2) Save core draw-level tables
save(all_draws_ix_true, file = "01_Data/all_draws_ix_true.RData")
save(all_draws_hosp_true, file = "01_Data/all_draws_hosp_true.RData")
save(all_draws_daly_true, file = "01_Data/all_draws_daly_true.RData")
save(all_draws_fatal_true, file = "01_Data/all_draws_fatal_true.RData")
save(all_draws_sae_true, file = "01_Data/all_draws_sae_true.RData")


# =============================================================================
# SECTION 08. BRR pipeline (draw-level benefit vs risk; no object mixing)
# Everything is created as *_true and never overwritten.
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 08A. Helper functions and lookup maps
# -----------------------------------------------------------------------------

# -------------------------------
# Step 8.0) Optional safety cleanup for conflicting objects
# -------------------------------
# rm(list = ls(pattern = "^(benefit_draws|risk_draws|joint|draw_level_xy|summary_long|summary_df|brr_)"))
# rm(list = ls(pattern = "^(tot_vacc_map|tot_vacc_map2)$"))

# -------------------------------
# Step 8.1) Helper: scenario -> AgeCat mapping (integer scenarios)
# -------------------------------
map_scenario_agecat_int <- function(scen_int) {
  dplyr::case_when(
    scen_int == 1L ~ "1-11",
    scen_int == 2L ~ "12-17",
    scen_int == 3L ~ "18-64",
    scen_int == 4L ~ "65+",
    TRUE ~ NA_character_
  )
}

# -------------------------------
# Step 8.2) Build tot_vacc_map_true from national vaccine data (cov50 only)
# Output keys: Scenario (int), VE, AgeCat, tot_vacc_grp
# -------------------------------
tot_vacc_map_true <- combined_nnv_df_region_coverage_model %>%
  filter(VC == "cov50") %>%
  mutate(
    AgeCat = case_when(
      AgeGroup %in% 2:4   ~ "1-11",
      AgeGroup == 5       ~ "12-17",
      AgeGroup %in% 6:15  ~ "18-64",
      AgeGroup %in% 16:20 ~ "65+",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(AgeCat)) %>%
  group_by(scenario, region, AgeCat, VE) %>%
  summarise(tot_vacc_grp = sum(tot_vacc, na.rm = TRUE), .groups = "drop") %>%
  # keep only the targeted AgeCat for each scenario
  mutate(
    target = case_when(
      scenario == "Scenario_1" & AgeCat == "1-11"  ~ 1L,
      scenario == "Scenario_2" & AgeCat == "12-17" ~ 1L,
      scenario == "Scenario_3" & AgeCat == "18-64" ~ 1L,
      scenario == "Scenario_4" & AgeCat == "65+"   ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  dplyr::filter(target == 1L) %>%
  transmute(
    scenario = as.integer(gsub("Scenario_", "", scenario)),
    AgeCat, VE, tot_vacc_grp, region
  )%>%
  dplyr::rename(Region   = region,
                Scenario = scenario)

# -------------------------------
# Step 8.3) Build DALY parameter draws and benefit draws
# -------------------------------
make_averted_draws_true <- function(df_true, outcome_name, coverage_keep = "cov50") {
  df_true %>%
    dplyr::filter(Coverage == coverage_keep) %>%
    dplyr::mutate(
      Scenario = as.integer(Scenario),
      outcome  = outcome_name,
      baseline = total_pre,
      post     = total_post,
      averted  = total_pre - total_post
    ) %>%
    dplyr::select(
      draw_id, Region, VE, Coverage, Scenario, outcome,
      baseline, post, averted
    )
}

# -------------------------------
# Step 8.3.a) Reporting/plot helper functions
# -------------------------------
make_pr_gt1_wide <- function(ceac_ob, brr_type_filter) {
  ceac_ob %>%
    filter(
      abs(log10(threshold)) < 1e-12,
      brr_type == brr_type_filter,
      RR_seropos == 0
    ) %>%
    transmute(
      Outcome     = outcome,
      Setting     = setting,
      `Age group` = AgeCat,
      VE_col      = VE_label,
      pr_fmt      = sprintf("%.0f%%", 100 * p_accept)
    ) %>%
    pivot_wider(
      names_from  = VE_col,
      values_from = pr_fmt,
      names_glue  = "{VE_col} Pr(BRR>1)"
    )
}

# Format median and uncertainty interval as "med (lo-hi)"
fmt_ci <- function(med, lo, hi) {
  sprintf("%.2f (%.2f–%.2f)", med, lo, hi)
}


## BRR probability curve and data generation function
plot_brr_ceac_outbreak_ve <- function(ceac_df,
                                      target_outcome = "DALY",
                                      setting_levels = c("Low", "Moderate", "High"),
                                      age_levels     = c("18-64", "65+"),
                                      ve_levels      = NULL,
                                      brr_type_labels = c("brr_base" = "Base (no serostatus)",
                                                          "brr_adj"  = "Adjusted (serostatus)")) {
  df <- ceac_df %>%
    filter(
      outcome    == target_outcome,
      RR_seropos == 0
    ) %>%
    mutate(
      setting   = factor(setting,   levels = setting_levels),
      AgeCat    = factor(AgeCat,    levels = age_levels),
      brr_type  = factor(brr_type,  levels = names(brr_type_labels),
                         labels = brr_type_labels)
    )
  
  if (!is.null(ve_levels)) {
    df <- df %>% mutate(VE_label = factor(VE_label, levels = ve_levels))
  }
  
  ggplot(df, aes(x = threshold, y = p_accept,
                 colour   = AgeCat,
                 linetype = brr_type)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.4) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_linetype_manual(values = c("Base (no serostatus)"  = "solid",
                                     "Adjusted (serostatus)" = "dashed")) +
    facet_grid(setting ~ VE_label) +
    labs(
      x        = "Benefit-risk ratio (BRR)",
      y        = "Probability (BRR > t)",
      title    = paste0("BRR acceptability curve: ", target_outcome),
      colour   = "Age group",
      linetype = "BRR type"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


make_brr_ceac_outbreak <- function(brr_long,
                                   thresholds = 10^seq(-1, 1, by = 0.02),
                                   group_vars = c("setting","VE_label","AgeCat","outcome")) {
  brr_long %>%
    tidyr::crossing(threshold = thresholds) %>%
    group_by(across(all_of(group_vars)), threshold) %>%
    summarise(
      p_accept = mean(brr > threshold),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# SECTION 08B. Core BRR computation tables
# -----------------------------------------------------------------------------



# -------------------------------
# Step 8.6) Build benefit base
# -------------------------------
benefit_base_true <- benefit_draws_true %>%
  dplyr::mutate(
    Scenario = as.integer(Scenario),
    AgeCat   = map_scenario_agecat_int(Scenario),
    risk_band = dplyr::if_else(AgeCat == "65+", "65+", "u65")
  ) %>%
  dplyr::left_join(daly_pars_true, by = "draw_id")



### Step 8.10) Added on 2026-04-04: serostatus-specific risk decomposition -----
age_map <- data.frame(
  age_index = 1:20,
  age_gr = c(
    "<1", "1-4", "5-9", "10-11", "12-17", "18-19", "20-24", "25-29",
    "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
    "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
  ),
  AgeCat = c(
    "1-11", "1-11", "1-11", "1-11",
    "12-17",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "65+", "65+", "65+", "65+", "65+"
  )
)

### Step 8.10.1)
### Baseline harms for 65+ reproduce the original implementation after aligning
### the originally assigned risk draw.
### Under-65 baseline harms intentionally differ because vaccine-attributable
### death risk is fixed at zero in this revised implementation.
calc_q_seromix_for_scenarios_agecat <- function(sim_region_ve_cov, age_map) {
  n_scenarios <- length(sim_region_ve_cov)
  out <- vector("list", n_scenarios)
  
  for (sc in seq_len(n_scenarios)) {
    raw_alloc_array <- sim_region_ve_cov[[sc]]$sim_result$raw_allocation_array
    vacc_to_S_array <- sim_region_ve_cov[[sc]]$sim_result$vacc_to_S_array
    
    stopifnot(all(dim(raw_alloc_array) == dim(vacc_to_S_array)))
    
    n_draws <- dim(raw_alloc_array)[3]
    sc_df <- vector("list", n_draws)
    
    for (d in seq_len(n_draws)) {
      df <- data.frame(
        age_index = seq_len(dim(raw_alloc_array)[1]),
        total_vacc_age = rowSums(raw_alloc_array[, , d, drop = FALSE], na.rm = TRUE),
        seroneg_vacc_age = rowSums(vacc_to_S_array[, , d, drop = FALSE], na.rm = TRUE)
      ) %>%
        dplyr::mutate(
          seropos_vacc_age = pmax(0, total_vacc_age - seroneg_vacc_age)
        ) %>%
        dplyr::left_join(age_map, by = "age_index") %>%
        dplyr::group_by(AgeCat) %>%
        dplyr::summarise(
          total_vacc_age = sum(total_vacc_age, na.rm = TRUE),
          seroneg_vacc_age = sum(seroneg_vacc_age, na.rm = TRUE),
          seropos_vacc_age = sum(seropos_vacc_age, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Scenario = sc,
          draw_id = d,
          q_seroneg_vacc = ifelse(total_vacc_age > 0, seroneg_vacc_age / total_vacc_age, NA_real_),
          q_seropos_vacc = ifelse(total_vacc_age > 0, seropos_vacc_age / total_vacc_age, NA_real_)
        )
      
      sc_df[[d]] <- df
    }
    
    out[[sc]] <- dplyr::bind_rows(sc_df)
  }
  
  dplyr::bind_rows(out)
}

q_all_regions <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      
      q_tmp <- calc_q_seromix_for_scenarios_agecat(
        sim_region_ve_cov = cov_list$scenario_result,
        age_map = age_map
      )
      
      q_tmp %>%
        dplyr::mutate(
          Region = region_name,
          VE = ve_name,
          Coverage = cov_name
        ) %>%
        dplyr::select(
          Region, VE, Coverage, Scenario, draw_id, AgeCat,
          total_vacc_age, seroneg_vacc_age, seropos_vacc_age,
          q_seroneg_vacc, q_seropos_vacc
        )
    })
  })
})


### Step 8.10.2) Generate risk table from LHS sample
risk_draw_df <- bind_rows(
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "1-11",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_u65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]
  ),
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "12-17",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_u65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]
  ),
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "18-64",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_u65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]
  ),
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "65+",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_65"]
  )
)

## Step 8.10.3) Combine risk and regional composition data
rr_vals <- c(1.0, 0.5, 0.1, 0.0)

risk_components_all <- q_all_regions %>%
  tidyr::crossing(RR_seropos = rr_vals) %>%
  dplyr::left_join(
    risk_draw_df,
    by = c("draw_id", "AgeCat")
  )

risk_components_all <- risk_components_all %>%
  dplyr::mutate(
    # conditional risk by serostatus
    p_sae_vacc_seroneg   = p_sae_vacc_base,
    p_death_vacc_seroneg = p_death_vacc_base,
    
    p_sae_vacc_seropos   = p_sae_vacc_base * RR_seropos,
    p_death_vacc_seropos = p_death_vacc_base * RR_seropos,
    
    # contribution to vaccinated cohort
    p_sae_vacc_seroneg_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seroneg_vacc * p_sae_vacc_seroneg,
      0
    ),
    p_death_vacc_seroneg_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seroneg_vacc * p_death_vacc_seroneg,
      0
    ),
    
    p_sae_vacc_seropos_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seropos_vacc * p_sae_vacc_seropos,
      0
    ),
    p_death_vacc_seropos_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seropos_vacc * p_death_vacc_seropos,
      0
    ),
    
    # total adjusted
    p_sae_vacc_adj = p_sae_vacc_seroneg_contrib + p_sae_vacc_seropos_contrib,
    p_death_vacc_adj = p_death_vacc_seroneg_contrib + p_death_vacc_seropos_contrib
  )

# Step 8.10.4) Estimate SAE and death per 10,000 vaccinated
risk_components_all <- risk_components_all %>%
  dplyr::mutate(
    # baseline
    sae_10k_base   = 1e4 * p_sae_vacc_base + 1e4 * p_death_vacc_base,
    death_10k_base = 1e4 * p_death_vacc_base,
    
    # seronegative contribution
    sae_10k_seroneg   = 1e4 * p_sae_vacc_seroneg_contrib + 1e4 * p_death_vacc_seroneg_contrib,
    death_10k_seroneg = 1e4 * p_death_vacc_seroneg_contrib,
    
    # seropositive contribution
    sae_10k_seropos   = 1e4 * p_sae_vacc_seropos_contrib + 1e4 * p_death_vacc_seropos_contrib,
    death_10k_seropos = 1e4 * p_death_vacc_seropos_contrib,
    
    # adjusted total
    sae_10k_adj   = sae_10k_seroneg + sae_10k_seropos,
    death_10k_adj = death_10k_seroneg + death_10k_seropos
  )

# Step 8.10.5) Join serostatus risk decomposition to benefit data
risk_join_all <- risk_components_all %>%
  dplyr::select(
    Region, VE, Coverage, Scenario, draw_id, AgeCat, RR_seropos,
    total_vacc_age, q_seroneg_vacc, q_seropos_vacc,
    p_sae_vacc_base, p_death_vacc_base,
    p_sae_vacc_seroneg, p_death_vacc_seroneg,
    p_sae_vacc_seropos, p_death_vacc_seropos,
    p_sae_vacc_seroneg_contrib, p_death_vacc_seroneg_contrib,
    p_sae_vacc_seropos_contrib, p_death_vacc_seropos_contrib,
    p_sae_vacc_adj, p_death_vacc_adj,
    sae_10k_base, death_10k_base,
    sae_10k_seroneg, death_10k_seroneg,
    sae_10k_seropos, death_10k_seropos,
    sae_10k_adj, death_10k_adj
  )

# Step 8.10.6) Build draw-level burden table
benefit_draw_df <- benefit_base_true %>%
  dplyr::left_join(
    tot_vacc_map_true,
    by = c("Region", "Scenario", "VE", "AgeCat")
  ) %>%
  dplyr::mutate(
    Scenario = as.integer(Scenario),
    baseline_10k = (baseline / tot_vacc_grp) * 1e4,
    post_10k     = (post / tot_vacc_grp) * 1e4,
    averted_10k  = (averted / tot_vacc_grp) * 1e4,
    pct_reduction = dplyr::if_else(
      baseline > 0,
      100 * averted / baseline,
      NA_real_
    )
  )

benefit_draw_df <- tidyr::crossing(
  benefit_draw_df,
  RR_seropos = rr_vals
)

# join benefit + risk decomposition
draw_level_xy_serostatus <- benefit_draw_df %>%
  dplyr::left_join(
    risk_join_all,
    by = c("Region", "VE", "Coverage", "Scenario", "draw_id", "AgeCat", "RR_seropos")
  )

# Step 8.10.7) Compute vaccine-caused DALY components
draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    daly_10k_base = if (outcome == "DALY") {
      compute_daly_one(
        age_group      = AgeCat,
        sae_10k        = sae_10k_base,
        deaths_sae_10k = death_10k_base,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    daly_10k_seroneg = if (outcome == "DALY") {
      compute_daly_one(
        age_group      = AgeCat,
        sae_10k        = sae_10k_seroneg,
        deaths_sae_10k = death_10k_seroneg,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    daly_10k_seropos = if (outcome == "DALY") {
      compute_daly_one(
        age_group      = AgeCat,
        sae_10k        = sae_10k_seropos,
        deaths_sae_10k = death_10k_seropos,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    daly_10k_adj = if (outcome == "DALY") {
      daly_10k_seroneg + daly_10k_seropos
    } else {
      NA_real_
    }
  ) %>%
  dplyr::ungroup()

# Step 8.10.8) Add reporting labels
draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
  dplyr::mutate(
    VE_label = factor(
      VE,
      levels = c("VE0", "VE98.9"),
      labels = c("Disease blocking only", "Disease and infection blocking")
    ),
    setting = unname(setting_key[Region])
  )

draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
  dplyr::mutate(
    x_10k_base = dplyr::case_when(
      outcome == "SAE"   ~ sae_10k_base,
      outcome == "Death" ~ death_10k_base,
      outcome == "DALY"  ~ daly_10k_base,
      TRUE ~ NA_real_
    ),
    x_10k_adj = dplyr::case_when(
      outcome == "SAE"   ~ sae_10k_adj,
      outcome == "Death" ~ death_10k_adj,
      outcome == "DALY"  ~ daly_10k_adj,
      TRUE ~ NA_real_
    ),
    brr_base = dplyr::if_else(
      is.na(x_10k_base) | x_10k_base == 0,
      NA_real_,
      averted_10k / x_10k_base
    ),
    brr_adj = dplyr::if_else(
      is.na(x_10k_adj) | x_10k_adj == 0,
      NA_real_,
      averted_10k / x_10k_adj
    )
  )

## Step 8.10.9) Summary across all age groups (including 1-17)
brr_draw_summary_true <- draw_level_xy_serostatus %>% 
  dplyr::mutate(
    brr_base = ifelse(is.infinite(brr_base), NA, brr_base),
    brr_adj  = ifelse(is.infinite(brr_adj),  NA, brr_adj),
    setting  = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  dplyr::group_by(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting) %>%  
  dplyr::summarise(
    # BRR
    brr_base_med = quantile(brr_base, 0.50,  na.rm = TRUE),
    brr_base_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
    brr_base_hi  = quantile(brr_base, 0.975, na.rm = TRUE),
    
    brr_adj_med  = quantile(brr_adj,  0.50,  na.rm = TRUE),
    brr_adj_lo   = quantile(brr_adj,  0.025, na.rm = TRUE),
    brr_adj_hi   = quantile(brr_adj,  0.975, na.rm = TRUE),
    
    # Averted burden (outcome-specific; distinguished by outcome column)
    av_med = quantile(averted_10k, 0.50,  na.rm = TRUE),
    av_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
    av_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),
    
    # Caused — base
    ca_sae_base_med   = quantile(sae_10k_base,   0.50,  na.rm = TRUE),
    ca_sae_base_lo    = quantile(sae_10k_base,   0.025, na.rm = TRUE),
    ca_sae_base_hi    = quantile(sae_10k_base,   0.975, na.rm = TRUE),
    
    ca_death_base_med = quantile(death_10k_base, 0.50,  na.rm = TRUE),
    ca_death_base_lo  = quantile(death_10k_base, 0.025, na.rm = TRUE),
    ca_death_base_hi  = quantile(death_10k_base, 0.975, na.rm = TRUE),
    
    ca_daly_base_med  = quantile(daly_10k_base,  0.50,  na.rm = TRUE),
    ca_daly_base_lo   = quantile(daly_10k_base,  0.025, na.rm = TRUE),
    ca_daly_base_hi   = quantile(daly_10k_base,  0.975, na.rm = TRUE),
    
    # Caused — adj (serostatus-adjusted)
    ca_sae_adj_med    = quantile(sae_10k_adj,    0.50,  na.rm = TRUE),
    ca_sae_adj_lo     = quantile(sae_10k_adj,    0.025, na.rm = TRUE),
    ca_sae_adj_hi     = quantile(sae_10k_adj,    0.975, na.rm = TRUE),
    
    ca_death_adj_med  = quantile(death_10k_adj,  0.50,  na.rm = TRUE),
    ca_death_adj_lo   = quantile(death_10k_adj,  0.025, na.rm = TRUE),
    ca_death_adj_hi   = quantile(death_10k_adj,  0.975, na.rm = TRUE),
    
    ca_daly_adj_med   = quantile(daly_10k_adj,   0.50,  na.rm = TRUE),
    ca_daly_adj_lo    = quantile(daly_10k_adj,   0.025, na.rm = TRUE),
    ca_daly_adj_hi    = quantile(daly_10k_adj,   0.975, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    scenario  = Scenario,
    age_group = AgeCat
  ) %>%
  dplyr::arrange(outcome, Scenario, setting, AgeCat, VE_label, RR_seropos)

## Step 8.10.10) Filter to >=18 years and RR_seropos = 0
brr_draw_summary_true_filtered <- brr_draw_summary_true %>%
  dplyr::filter(
    AgeCat %in% c("18-64", "65+"),
    RR_seropos == 0.0
  ) %>%
  dplyr::select(-scenario, -age_group)  


## Step 8.10.11) Build wide-format reporting tables
# Step 8.10.12) Merge caused-burden columns into wide table
ca_summary <- brr_draw_summary_true_filtered %>%
  filter(outcome == "DALY") %>%  
  select(Scenario, AgeCat, VE_label, RR_seropos, setting,
         starts_with("ca_"))

brr_table_long_serostatus <- brr_draw_summary_true_filtered %>%
  select(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting,
         brr_base_med, brr_base_lo, brr_base_hi,
         brr_adj_med,  brr_adj_lo,  brr_adj_hi,
         av_med, av_lo, av_hi) %>%
  left_join(ca_summary,
            by = c("Scenario", "AgeCat", "VE_label", "RR_seropos", "setting")) %>%
  mutate(
    brr_base_formatted = sprintf("%.2f (%.2f–%.2f)", brr_base_med, brr_base_lo, brr_base_hi),
    brr_adj_formatted  = sprintf("%.2f (%.2f–%.2f)", brr_adj_med,  brr_adj_lo,  brr_adj_hi),
    av_formatted       = sprintf("%.2f (%.2f–%.2f)", av_med,       av_lo,       av_hi),
    ca_sae_base_formatted   = sprintf("%.2f (%.2f–%.2f)", ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
    ca_death_base_formatted = sprintf("%.2f (%.2f–%.2f)", ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
    ca_daly_base_formatted  = sprintf("%.2f (%.2f–%.2f)", ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
    ca_sae_adj_formatted    = sprintf("%.2f (%.2f–%.2f)", ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
    ca_death_adj_formatted  = sprintf("%.2f (%.2f–%.2f)", ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
    ca_daly_adj_formatted   = sprintf("%.2f (%.2f–%.2f)", ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
  )


brr_long_serostatus <- draw_level_xy_serostatus %>%
  mutate(
    outcome = recode(outcome, "sae" = "SAE", "death" = "Death", "daly" = "DALY"),
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    AgeCat  = factor(AgeCat,  levels = c("18-64", "65+"))
  ) %>%
  filter(
    is.finite(brr_base), brr_base > 0,
    is.finite(brr_adj),  brr_adj  > 0
  ) %>%
  transmute(
    Region, setting, VE_label, AgeCat, outcome, RR_seropos,
    brr_base,
    brr_adj
  )

brr_long_serostatus_long <- brr_long_serostatus %>%
  filter(!is.na(AgeCat)) %>%
  pivot_longer(
    cols      = c(brr_base, brr_adj),
    names_to  = "brr_type",
    values_to = "brr"
  ) %>%
  filter(is.finite(brr), brr > 0)

brr_range <- brr_long_serostatus_long %>%
  filter(RR_seropos == 0) %>%
  summarise(
    min_brr = min(brr, na.rm = TRUE),
    max_brr = max(brr, na.rm = TRUE)
  )

brr_min <- brr_range$min_brr
brr_max <- brr_range$max_brr

thresholds_auto <- 10^seq(
  floor(log10(brr_min)),   
  ceiling(log10(brr_max)), 
  by = 0.02
)
ceac_ob <- make_brr_ceac_outbreak(
  brr_long_serostatus_long,
  group_vars = c("setting", "VE_label", "AgeCat", "outcome", "RR_seropos", "brr_type")
)

pr_gt1_wide_setting_base <- make_pr_gt1_wide(ceac_ob, "brr_base")
pr_gt1_wide_setting_adj  <- make_pr_gt1_wide(ceac_ob, "brr_adj")

ca_summary <- brr_draw_summary_true_filtered %>%
  filter(RR_seropos == 0, outcome == "DALY") %>%
  select(Scenario, AgeCat, VE_label, setting,
         ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi,
         ca_death_base_med, ca_death_base_lo, ca_death_base_hi,
         ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi,
         ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi,
         ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi,
         ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)

brr_draw_for_wide <- brr_draw_summary_true_filtered %>%
  filter(RR_seropos == 0) %>%
  select(outcome, Scenario, AgeCat, VE_label, setting,
         av_med,       av_lo,       av_hi,
         brr_base_med, brr_base_lo, brr_base_hi,
         brr_adj_med,  brr_adj_lo,  brr_adj_hi) %>%
  left_join(ca_summary, by = c("Scenario", "AgeCat", "VE_label", "setting"))

brr_table_wide_serostatus <- brr_draw_for_wide %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    
    av_formatted            = fmt_ci(av_med,       av_lo,       av_hi),
    brr_base_formatted      = fmt_ci(brr_base_med, brr_base_lo, brr_base_hi),
    brr_adj_formatted       = fmt_ci(brr_adj_med,  brr_adj_lo,  brr_adj_hi),
    ca_sae_base_formatted   = fmt_ci(ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
    ca_death_base_formatted = fmt_ci(ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
    ca_daly_base_formatted  = fmt_ci(ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
    ca_sae_adj_formatted    = fmt_ci(ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
    ca_death_adj_formatted  = fmt_ci(ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
    ca_daly_adj_formatted   = fmt_ci(ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
  ) %>%
  select(outcome, Scenario, AgeCat, setting, VE_label,
         av_formatted, brr_base_formatted, brr_adj_formatted,
         ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
         ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted) %>%
  pivot_wider(
    names_from  = VE_label,
    values_from = c(av_formatted, brr_base_formatted, brr_adj_formatted,
                    ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
                    ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted),
    names_glue  = "{VE_label}_{.value}"
  ) %>%
  arrange(outcome, Scenario, setting, AgeCat) %>%
  rename(
    Outcome     = outcome,
    `Age group` = AgeCat,
    Setting     = setting
  )

brr_table_wide_serostatus <- brr_table_wide_serostatus %>%
  dplyr::rename(
    # Disease blocking only
    `DB (Averted)`           = `Disease blocking only_av_formatted`,
    `DB (BRR base)`          = `Disease blocking only_brr_base_formatted`,
    `DB (BRR adj)`           = `Disease blocking only_brr_adj_formatted`,
    `DB (SAE caused base)`   = `Disease blocking only_ca_sae_base_formatted`,
    `DB (SAE caused adj)`    = `Disease blocking only_ca_sae_adj_formatted`,
    `DB (Death caused base)` = `Disease blocking only_ca_death_base_formatted`,
    `DB (Death caused adj)`  = `Disease blocking only_ca_death_adj_formatted`,
    `DB (DALY caused base)`  = `Disease blocking only_ca_daly_base_formatted`,
    `DB (DALY caused adj)`   = `Disease blocking only_ca_daly_adj_formatted`,
    # Disease and infection blocking
    `DIB (Averted)`           = `Disease and infection blocking_av_formatted`,
    `DIB (BRR base)`          = `Disease and infection blocking_brr_base_formatted`,
    `DIB (BRR adj)`           = `Disease and infection blocking_brr_adj_formatted`,
    `DIB (SAE caused base)`   = `Disease and infection blocking_ca_sae_base_formatted`,
    `DIB (SAE caused adj)`    = `Disease and infection blocking_ca_sae_adj_formatted`,
    `DIB (Death caused base)` = `Disease and infection blocking_ca_death_base_formatted`,
    `DIB (Death caused adj)`  = `Disease and infection blocking_ca_death_adj_formatted`,
    `DIB (DALY caused base)`  = `Disease and infection blocking_ca_daly_base_formatted`,
    `DIB (DALY caused adj)`   = `Disease and infection blocking_ca_daly_adj_formatted`
  )

pr_gt1_wide_serostatus_base <- pr_gt1_wide_setting_base %>%
  dplyr::rename_with(
    ~ gsub("Pr(BRR>1)", "Pr(BRR_base>1)", .x, fixed = TRUE),
    .cols = -c(Outcome, Setting, `Age group`)
  )

pr_gt1_wide_serostatus_adj <- pr_gt1_wide_setting_adj %>%
  dplyr::rename_with(
    ~ gsub("Pr(BRR>1)", "Pr(BRR_adj>1)", .x, fixed = TRUE),
    .cols = -c(Outcome, Setting, `Age group`)
  )

brr_table_wide_serostatus2 <- brr_table_wide_serostatus %>%
  left_join(pr_gt1_wide_serostatus_base, by = c("Outcome", "Setting", "Age group")) %>%
  left_join(pr_gt1_wide_serostatus_adj,  by = c("Outcome", "Setting", "Age group"))

brr_table_wide_serostatus2 <- brr_table_wide_serostatus2 %>%
  relocate(`DB (Averted)`,                                  .after = `Age group`) %>%
  relocate(`DB (BRR base)`,                                 .after = `DB (Averted)`) %>%
  relocate(`DB (BRR adj)`,                                  .after = `DB (BRR base)`) %>%
  relocate(`Disease blocking only Pr(BRR_base>1)`,          .after = `DB (BRR adj)`) %>%
  relocate(`Disease blocking only Pr(BRR_adj>1)`,           .after = `Disease blocking only Pr(BRR_base>1)`) %>%
  relocate(`DB (SAE caused base)`,                          .after = `Disease blocking only Pr(BRR_adj>1)`) %>%
  relocate(`DB (SAE caused adj)`,                           .after = `DB (SAE caused base)`) %>%
  relocate(`DB (Death caused base)`,                        .after = `DB (SAE caused adj)`) %>%
  relocate(`DB (Death caused adj)`,                         .after = `DB (Death caused base)`) %>%
  relocate(`DB (DALY caused base)`,                         .after = `DB (Death caused adj)`) %>%
  relocate(`DB (DALY caused adj)`,                          .after = `DB (DALY caused base)`) %>%
  relocate(`DIB (Averted)`,                                 .after = `DB (DALY caused adj)`) %>%
  relocate(`DIB (BRR base)`,                                .after = `DIB (Averted)`) %>%
  relocate(`DIB (BRR adj)`,                                 .after = `DIB (BRR base)`) %>%
  relocate(`Disease and infection blocking Pr(BRR_base>1)`, .after = `DIB (BRR adj)`) %>%
  relocate(`Disease and infection blocking Pr(BRR_adj>1)`,  .after = `Disease and infection blocking Pr(BRR_base>1)`) %>%
  relocate(`DIB (SAE caused base)`,                         .after = `Disease and infection blocking Pr(BRR_adj>1)`) %>%
  relocate(`DIB (SAE caused adj)`,                          .after = `DIB (SAE caused base)`) %>%
  relocate(`DIB (Death caused base)`,                       .after = `DIB (SAE caused adj)`) %>%
  relocate(`DIB (Death caused adj)`,                        .after = `DIB (Death caused base)`) %>%
  relocate(`DIB (DALY caused base)`,                        .after = `DIB (Death caused adj)`) %>%
  relocate(`DIB (DALY caused adj)`,                         .after = `DIB (DALY caused base)`) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ tidyr::replace_na(as.character(.x), "beneficial")
    )
  )

# -----------------------------------------------------------------------------
# SECTION 08C. Reporting table assembly and document export
# -----------------------------------------------------------------------------

# Step 8.10.13) Build publication-ready BRR table and export to Word
# 1) Extract and rename columns for "Disease blocking only" (base + adjusted BRR/risk)
part1 <- brr_table_wide_serostatus2 %>%
  dplyr::mutate(
    Risk_base = dplyr::case_when(
      Outcome == "SAE"   ~ `DB (SAE caused base)`,
      Outcome == "Death" ~ `DB (Death caused base)`,
      Outcome == "DALY"  ~ `DB (DALY caused base)`,
      TRUE ~ NA_character_
    ),
    Risk_adj = dplyr::case_when(
      Outcome == "SAE"   ~ `DB (SAE caused adj)`,
      Outcome == "Death" ~ `DB (Death caused adj)`,
      Outcome == "DALY"  ~ `DB (DALY caused adj)`,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `DB (Averted)`,
    Risk_base,
    Risk_adj,
    BRR_base = `DB (BRR base)`,
    BRR_adj  = `DB (BRR adj)`,
    prob_base = `Disease blocking only Pr(BRR_base>1)`,
    prob_adj  = `Disease blocking only Pr(BRR_adj>1)`
  ) %>%
  dplyr::mutate(mechanism = "Disease blocking only")

# 2) Extract and rename columns for "Disease and infection blocking"
part2 <- brr_table_wide_serostatus2 %>%
  dplyr::mutate(
    Risk_base = dplyr::case_when(
      Outcome == "SAE"   ~ `DIB (SAE caused base)`,
      Outcome == "Death" ~ `DIB (Death caused base)`,
      Outcome == "DALY"  ~ `DIB (DALY caused base)`,
      TRUE ~ NA_character_
    ),
    Risk_adj = dplyr::case_when(
      Outcome == "SAE"   ~ `DIB (SAE caused adj)`,
      Outcome == "Death" ~ `DIB (Death caused adj)`,
      Outcome == "DALY"  ~ `DIB (DALY caused adj)`,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `DIB (Averted)`,
    Risk_base,
    Risk_adj,
    BRR_base = `DIB (BRR base)`,
    BRR_adj  = `DIB (BRR adj)`,
    prob_base = `Disease and infection blocking Pr(BRR_base>1)`,
    prob_adj  = `Disease and infection blocking Pr(BRR_adj>1)`
  ) %>%
  dplyr::mutate(mechanism = "Disease and infection blocking")

# 3) Stack both mechanisms into one long-format reporting table
brr_table_final_long <- dplyr::bind_rows(part1, part2) %>%
  dplyr::mutate(
    Setting = factor(Setting, levels = c("High", "Moderate", "Low"))
  ) %>%
  dplyr::arrange(Outcome, Setting, `Age group`) %>%
  dplyr::select(
    Outcome, Setting, `Age group`, mechanism,
    Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj, prob_base, prob_adj
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj),
      ~ {
        x <- gsub(" \\(", "\n(", .x, fixed = FALSE)
        dplyr::case_when(
          x %in% c(
            "NA (NA-NA)", "NA (NA–NA)",
            "NA\n(NA-NA)", "NA\n(NA–NA)",
            "NA (NA NA)", "NA\n(NA NA)"
          ) ~ "beneficial",
          TRUE ~ x
        )
      }
    )
  )

ft_brr <- flextable::flextable(brr_table_final_long) %>%
  flextable::set_header_labels(
    Outcome     = "Outcome",
    Setting     = "Setting",
    `Age group` = "Age group",
    mechanism   = "Vaccine protection\nmechanism",
    Benefit     = "Benefit:\nOutcomes averted\n(per 10,000)",
    Risk_base   = "Risk (base):\nOutcomes attributable\n(per 10,000)",
    Risk_adj    = "Risk (adjusted):\nOutcomes attributable\n(per 10,000)",
    BRR_base    = "BRR (base):\n(Prevented per 1 caused)",
    BRR_adj     = "BRR (adjusted):\n(Prevented per 1 caused)",
    prob_base   = "Probability\n(BRR_base > 1)\n(%)",
    prob_adj    = "Probability\n(BRR_adj > 1)\n(%)"
  ) %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:4, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Outcome", "Setting", "Age group")) %>%
  flextable::valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

ft_brr

doc <- officer::read_docx() %>%
  officer::body_add_par(
    "Benefit-Risk Ratio (BRR) by Outcome, Age group, and VE",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft_brr)

print(doc, target = "06_Results/BRR_table_ori_setting.docx")


# -----------------------------------------------------------------------------
# SECTION 08D. BRR probability curve plots
# -----------------------------------------------------------------------------

## BRR probability curves 

p_daly_ceac_outbreak <- plot_brr_ceac_outbreak_ve(ceac_ob, "DALY") +
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag          = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )


plot_sae   <- plot_brr_ceac_outbreak_ve(ceac_ob, "SAE")+
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag          = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

plot_death <- plot_brr_ceac_outbreak_ve(ceac_ob, "Death")+
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag          = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

