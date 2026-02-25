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
# make all draws for symptomatic cases
# draw base impact
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

### before scaling
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


# 1. symptomatic cases
## scaling rho
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

# 2. Hospitalisation 
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

# 3. Fatal
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

calc_total_hosp_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
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

all_draws_fatal_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool (필수)
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

# 4. DALY
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
    
    # precomputed 
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

scale_symp_by_rho_pre <- function(pre_list, rho_vec) {
  # pre_list: list of (age x week) matrices
  Map(function(mat, r) mat / r, pre_list, rho_vec)
}

scale_symp_by_rho_post <- function(post_arr, rho_vec) {
  # post_arr: (age x week x draw)
  for (i in seq_along(rho_vec)) post_arr[,,i] <- post_arr[,,i] / rho_vec[i]
  post_arr
}

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

# 5. SAE
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



save(all_draws_ix_true, file = "01_Data/all_draws_ix_true.RData")
save(all_draws_hosp_true, file = "01_Data/all_draws_hosp_true.RData")
save(all_draws_daly_true, file = "01_Data/all_draws_daly_true.RData")
save(all_draws_fatal_true, file = "01_Data/all_draws_fatal_true.RData")
save(all_draws_sae_true, file = "01_Data/all_draws_sae_true.RData")


# add setting key variable
# ============================================================
# BRR pipeline (NO object mixing) — draw-level benefit vs risk
# Everything is created as *_true and never overwritten.
# ============================================================

# -------------------------------
# 0) Safety: clean conflicting objects (optional)
# -------------------------------
# rm(list = ls(pattern = "^(benefit_draws|risk_draws|joint|draw_level_xy|summary_long|summary_df|brr_)"))
# rm(list = ls(pattern = "^(tot_vacc_map|tot_vacc_map2)$"))

# -------------------------------
# 1) Helper: scenario -> AgeCat mapping (integer scenarios)
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
# 2) Build tot_vacc_map_true from national vaccine data (cov50 only)
#    Output keys: scenario(int), VE, AgeCat, tot_vacc_grp
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
# 3) Build benefit_draws_true from *_true outcomes (cov50 only)
#    Input expected: all_draws_daly_true, all_draws_hosp_true, all_draws_fatal_true
#    Each must contain: draw_id, Region, VE, Coverage, Scenario, setting, total_pre, total_post
# -------------------------------
make_averted_draws_true <- function(df_true, outcome_name, coverage_keep = "cov50") {
  df_true %>%
    filter(Coverage == coverage_keep) %>%
    mutate(
      Scenario = as.integer(Scenario),
      outcome  = outcome_name,
      averted  = total_pre - total_post
    ) %>%
    dplyr::select(draw_id, Region, VE, Coverage, Scenario, outcome, averted)
}

benefit_draws_true <- bind_rows(
  make_averted_draws_true(all_draws_daly_true,  "DALY"),
  make_averted_draws_true(all_draws_sae_true,   "SAE"),
  make_averted_draws_true(all_draws_fatal_true, "Death")
)

# -------------------------------
# 4) Risk draws (independent): one table with everything we need
#    Input expected: lhs_sample (your LHS draws)
# -------------------------------
risk_draws_true <- lhs_sample %>%
  mutate(risk_id = row_number()) %>%
  transmute(
    risk_id,
    
    # SAE / Death probabilities
    p_sae_u65   = p_sae_vacc_u65,
    p_sae_65    = p_sae_vacc_65,
    p_death_u65 = p_death_vacc_u65,
    p_death_65  = p_death_vacc_65,
    
    # DALY parameters used by compute_daly_one()
    le_lost_1_11, le_lost_12_17, le_lost_18_64, le_lost_65,
    dw_hosp, dw_nonhosp, dw_subac, dw_chronic,
    dur_acute, dur_nonhosp, dur_subac, dur_6m, dur_12m, dur_30m,
    acute, subac, chr6m, chr12m, chr30m
  )

# -------------------------------
# 5) Attach independent risk to each benefit row (bootstrap-style)
#    IMPORTANT: set seed ONCE here for reproducibility
# -------------------------------
set.seed(1)
risk_idx_true <- sample(risk_draws_true$risk_id, size = nrow(benefit_draws_true), replace = TRUE)

joint_true <- benefit_draws_true %>%
  mutate(risk_id = risk_idx_true) %>%
  left_join(risk_draws_true, by = "risk_id")

# -------------------------------
# 6) Build draw_level_xy_true (single pass, no repeats)
#    - derive AgeCat from Scenario
#    - join tot_vacc
#    - compute y_10k and x_10k for SAE/Death
#    - compute x_10k for DALY via compute_daly_one()
# -------------------------------
draw_level_xy_true <- joint_true %>%
  mutate(
    Scenario = as.integer(Scenario),
    AgeCat   = map_scenario_agecat_int(Scenario)
  ) %>%
  left_join(
    tot_vacc_map_true,
    by = c("Region", "Scenario", "VE", "AgeCat")
  ) %>%
  # hard checks: join must succeed
  { 
    if (any(is.na(.$tot_vacc_grp))) {
      bad <- . %>% filter(is.na(tot_vacc_grp)) %>% distinct(Scenario, AgeCat, VE) %>% head(20)
      stop("tot_vacc_map_true join failed for some keys. Examples:\n", paste(capture.output(print(bad)), collapse="\n"))
    }
    .
  } %>%
  mutate(
    # benefit per 10k vaccinated
    y_10k = (averted / tot_vacc_grp) * 1e4,
    
    # choose probabilities by AgeCat
    p_sae   = if_else(AgeCat == "65+", p_sae_65,   p_sae_u65),
    p_death = if_else(AgeCat == "65+", p_death_65, p_death_u65),
    
    # compute input units for compute_daly_one(): per 10k
    sae_10k        = p_sae   * 1e4,
    deaths_sae_10k = p_death * 1e4
  ) %>%
  rowwise() %>%
  mutate(
    # DALY risk per 10k vaccinated, from SAE & vaccine-death per 10k
    x_daly_10k = if (outcome == "DALY") {
      compute_daly_one(
        age_group      = AgeCat,
        sae_10k        = sae_10k,
        deaths_sae_10k = deaths_sae_10k,
        draw_pars = as.list(list(
          le_lost_1_11 = le_lost_1_11, le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64, le_lost_65 = le_lost_65,
          dw_hosp = dw_hosp, dw_nonhosp = dw_nonhosp, dw_subac = dw_subac, dw_chronic = dw_chronic,
          dur_acute = dur_acute, dur_nonhosp = dur_nonhosp, dur_subac = dur_subac,
          dur_6m = dur_6m, dur_12m = dur_12m, dur_30m = dur_30m,
          acute = acute, subac = subac, chr6m = chr6m, chr12m = chr12m, chr30m = chr30m
        ))
      )$daly_sae
    } else NA_real_,
    
    # x-axis (excess adverse outcomes) per 10k vaccinated
    x_10k = case_when(
      outcome == "SAE"   ~ sae_10k,
      outcome == "Death" ~ deaths_sae_10k,
      outcome == "DALY"  ~ x_daly_10k,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# -------------------------------
# 7) VE label (plot-friendly)
# -------------------------------
draw_level_xy_true <- draw_level_xy_true %>%
  mutate(
    VE_label = factor(VE,
                      levels = c("VE0", "VE98.9"),
                      labels = c("Disease blocking only", "Disease and infection blocking")
    )
  )%>%
  mutate(setting = unname(setting_key[Region]))
# -------------------------------
# 8) Summary for plotting (x/y quantiles)
#    Output matches your create_br_plot() expectation:
#    columns: outcome, Scenario, AgeCat, VE_label, x_lo/x_med/x_hi, y_lo/y_med/y_hi
# -------------------------------
summary_long_true <- draw_level_xy_true %>%
  group_by(outcome, Scenario, AgeCat, VE_label) %>%
  summarise(
    x_lo  = quantile(x_10k, 0.025, na.rm = TRUE),
    x_med = quantile(x_10k, 0.50,  na.rm = TRUE),
    x_hi  = quantile(x_10k, 0.975, na.rm = TRUE),
    y_lo  = quantile(y_10k, 0.025, na.rm = TRUE),
    y_med = quantile(y_10k, 0.50,  na.rm = TRUE),
    y_hi  = quantile(y_10k, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------
# 9) BRR summary + table (med + 95% UI)
# -------------------------------
brr_draw_summary_true <- draw_level_xy_true %>%
  mutate(brr = y_10k / pmax(x_10k, 1e-12)) %>%
  group_by(outcome, Scenario, AgeCat, VE_label) %>%
  summarise(
    brr_med = quantile(brr, 0.50,  na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario  = Scenario,
    age_group = AgeCat
  )

brr_table_long_true <- brr_draw_summary_true %>%
  mutate(
    brr_formatted = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi),
    VE_col = as.character(VE_label)
  ) %>%
  dplyr::select(outcome, scenario, age_group, VE_col, brr_formatted)

brr_table_wide_true <- brr_table_long_true %>%
  pivot_wider(names_from = VE_col, values_from = brr_formatted) %>%
  arrange(outcome, scenario, age_group) %>%
  dplyr::select(-scenario) %>%
  dplyr::rename(`Outcome` = outcome, `Age group` = age_group)

idx_outcome_true <- table(brr_table_wide_true$`Outcome`)

kable(
  brr_table_wide_true,
  format  = "html",
  caption = "Benefit–Risk Ratio (BRR) by Outcome, Scenario, Age Group, and VE",
  align   = "l"
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size  = 12
  ) %>%
  pack_rows(index = idx_outcome_true) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)

# -------------------------------
# 10) by setting
# -------------------------------
summary_long_setting <- draw_level_xy_true %>%
  mutate(
    setting = factor(setting, levels = c("Low","Moderate","High"))
  )%>%
  group_by(outcome, Scenario, AgeCat, VE_label, setting) %>%
  summarise(
    x_lo  = quantile(x_10k, 0.025, na.rm = TRUE),
    x_med = quantile(x_10k, 0.50,  na.rm = TRUE),
    x_hi  = quantile(x_10k, 0.975, na.rm = TRUE),
    y_lo  = quantile(y_10k, 0.025, na.rm = TRUE),
    y_med = quantile(y_10k, 0.50,  na.rm = TRUE),
    y_hi  = quantile(y_10k, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, Scenario, setting, AgeCat, VE_label)

# -------------------------------
# 9) BRR summary + table (med + 95% UI)
# -------------------------------
brr_draw_summary_setting <- draw_level_xy_true %>%
  mutate(brr = y_10k / pmax(x_10k, 1e-12)) %>%
  group_by(outcome, Scenario, AgeCat, VE_label, setting) %>%
  summarise(
    brr_med = quantile(brr, 0.50,  na.rm = TRUE),
    brr_lo  = quantile(brr, 0.025, na.rm = TRUE),
    brr_hi  = quantile(brr, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario  = Scenario,
    age_group = AgeCat
  )

brr_table_long_setting <- brr_draw_summary_setting %>%
  mutate(
    brr_formatted = sprintf("%.2f [%.2f–%.2f]", brr_med, brr_lo, brr_hi),
    VE_col = as.character(VE_label)
  ) %>%
  dplyr::select(outcome, scenario, setting, age_group, VE_col, brr_formatted)

brr_table_wide_setting <- brr_table_long_setting %>%
  mutate(
    setting = factor(setting, levels = c("Low","Moderate","High"))
  ) %>%
  pivot_wider(names_from = VE_col, values_from = brr_formatted) %>%
  arrange(outcome, scenario, setting, age_group) %>%
  dplyr::select(-scenario) %>%
  dplyr::rename(`Outcome` = outcome, `Age group` = age_group, `Setting` = setting)

idx_outcome_true <- table(brr_table_wide_setting$`Outcome`)

kable(
  brr_table_wide_setting,
  format  = "html",
  caption = "Benefit–Risk Ratio (BRR) by Outcome, Scenario, Age Group, and VE",
  align   = "l"
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size  = 12
  ) %>%
  pack_rows(index = idx_outcome_true) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE)

# ------------------------------------------------------------------------------
# Final Visualization for assessment plot
# ------------------------------------------------------------------------------
log_min <- -2
log_max <- 2
log_range <- seq(log_min, log_max, by = 1)
brr_labels <- c("0.01", "0.1", "1", "10", "100")

create_br_plot <- function(data, target_outcome, log_min = -1, log_max = 3,
                           grid_n = 150, pad = 1.10,
                           show_prop = TRUE) {
  
  plot_data <- data %>%
    filter(outcome == target_outcome)
  
  panel_limits <- plot_data %>%
    group_by(setting, AgeCat) %>%
    summarise(
      x_max_raw = max(x_hi, na.rm = TRUE),
      y_max_raw = max(y_hi, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_max_raw = ifelse(is.finite(x_max_raw), x_max_raw, NA_real_),
      y_max_raw = ifelse(is.finite(y_max_raw), y_max_raw, NA_real_),
      x_min = 1e-6,
      y_min = 0,
      x_max = x_max_raw * pad,
      y_max = y_max_raw * pad
    ) %>%
    filter(!is.na(x_max), !is.na(y_max))
  
  bg_grid_specific <- panel_limits %>%
    rowwise() %>%
    do({
      panel <- .
      x_seq <- seq(panel$x_min, panel$x_max, length.out = grid_n)
      y_seq <- seq(panel$y_min, panel$y_max, length.out = grid_n)
      
      grid_df <- expand.grid(x = x_seq, y = y_seq)
      
      grid_df$brr <- grid_df$y / grid_df$x
      grid_df$log10_brr <- log10(grid_df$brr)
      
      grid_df$log10_brr[!is.finite(grid_df$log10_brr)] <- NA_real_
      grid_df$log10_brr <- pmax(pmin(grid_df$log10_brr, log_max), log_min)
      
      grid_df$is_fav <- !is.na(grid_df$log10_brr) & grid_df$log10_brr > 0
      
      grid_df$setting <- panel$setting
      grid_df$AgeCat  <- panel$AgeCat
      grid_df
    }) %>%
    ungroup() %>%
    filter(!is.na(log10_brr))
  
  panel_prop <- bg_grid_specific %>%
    group_by(setting, AgeCat) %>%
    summarise(
      prop_fav = mean(is_fav, na.rm = TRUE),
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = ifelse(prop_fav < 0.005, "BRR>1: <1%",
                     sprintf("BRR>1: %.0f%%", 100 * prop_fav)),
      x_lab = x_min + 0.02 * (x_max - x_min),
      y_lab = y_max - 0.02 * (y_max - y_min)
    )
  
  p <- ggplot() +
    geom_raster(
      data = bg_grid_specific,
      aes(x = x, y = y, fill = log10_brr),
      alpha = 0.88,
      interpolate = TRUE
    ) +
    scale_fill_gradient2(
      name = "Benefit–risk ratio",
      low = "#ca0020", mid = "#f7f7f7", high = "#0571b0",
      midpoint = 0,
      limits = c(log_min, log_max),
      oob = scales::squish,
      breaks = log_range,
      labels = brr_labels,
      na.value = "white"
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      alpha = 0.55,
      linewidth = 0.9,
      colour = "grey35"
    ) +
    geom_errorbar(
      data = plot_data,
      aes(x = x_med, ymin = y_lo, ymax = y_hi, color = outcome),
      width = 0, linewidth = 0.3
    ) +
    geom_errorbarh(
      data = plot_data,
      aes(y = y_med, xmin = x_lo, xmax = x_hi, color = outcome),
      height = 0, linewidth = 0.3
    ) +
    geom_point(
      data = plot_data,
      aes(x = x_med, y = y_med, shape = VE_label, color = outcome),
      size = 2.9,
      fill = "white",
      stroke = 1.1,
      alpha = 0.95
    )
  
  # if (show_prop) {
  #   p <- p + geom_label(
  #     data = panel_prop,
  #     aes(x = x_lab, y = y_lab, label = label),
  #     inherit.aes = FALSE,
  #     hjust = 0, vjust = 1,
  #     size = 3,
  #     fill = "white",
  #     alpha = 0.85,
  #     colour = "black",
  #     label.size = 0.25,
  #     label.padding = unit(0.10, "lines")
  #   )
  # }
  
  p +
    facet_wrap(~ setting + AgeCat, scales = "free", ncol = 4) +
    scale_color_manual(
      name = "Outcome",
      values = c("SAE" = "#1B7F1B", "Death" = "#B8860B", "DALY" = "#A23B72")
    ) +
    scale_shape_manual(
      name = "Vaccine protection mechanism",
      values = c("Disease blocking only" = 21,
                 "Disease and infection blocking" = 24)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = paste("Benefit–risk Assessment:", target_outcome),
      x = "Vaccine attributable adverse outcome (per 10,000 vaccinated individuals)",
      y = "Outcomes averted by vaccination (per 10,000 vaccinated individuals)"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.6),
      strip.background = element_rect(fill = "gray95", linewidth = 0.6),
      axis.title = element_text(size = 11),
      axis.text  = element_text(size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      panel.spacing = unit(0.35, "lines")
    )
}

# by setting
plot_sae   <- create_br_plot(summary_long_setting, "SAE")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "D") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  ) 

plot_death <- create_br_plot(summary_long_setting, "Death")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "C") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  ) 

plot_daly  <- create_br_plot(summary_long_setting, "DALY")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  ) 

# national
#plot_sae   <- create_br_plot(summary_long_true, "SAE")
#plot_death <- create_br_plot(summary_long_true, "Death")
#plot_daly  <- create_br_plot(summary_long_true, "DALY")

plot_daly
plot_sae
plot_death


ggsave("06_Results/brr_daly_ori_setting.pdf", plot = plot_daly, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brr_death_ori_setting.pdf", plot = plot_death, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brr_sae_ori_setting.pdf", plot = plot_sae, width = 10, height = 8, device = cairo_pdf)

# ------------------------------------------------------------------------------
# Final Visualization for BRRAC
# ------------------------------------------------------------------------------

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

plot_brr_ceac_outbreak_ve <- function(ceac_df,
                                      target_outcome = "DALY",
                                      setting_levels = c("Low","Moderate","High"),
                                      age_levels = c("1-11","12-17","18-64","65+"),
                                      ve_levels = NULL) {
  
  df <- ceac_df %>%
    filter(outcome == target_outcome) %>%
    mutate(
      setting = factor(setting, levels = setting_levels),
      AgeCat  = factor(AgeCat,  levels = age_levels)
    )
  
  if (!is.null(ve_levels)) {
    df <- df %>% mutate(VE_label = factor(VE_label, levels = ve_levels))
  }
  
  ggplot(df, aes(x = threshold, y = p_accept, colour = AgeCat)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.6) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    facet_grid(setting ~ VE_label) +
    labs(
      x = "BRR threshold (t)",
      y = "Probability(BRR > t)",
      title = paste0("BRR acceptability curve: ", target_outcome),
      colour = "Age group"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

brr_long_ob <- draw_level_xy_true %>%
  mutate(
    brr = y_10k / pmax(x_10k, 1e-12),
    outcome = recode(outcome, "sae"="SAE","death"="Death","daly"="DALY"),
    setting = factor(setting, levels = c("Low","Moderate","High")),
    AgeCat  = factor(AgeCat, levels = c("1-11","12-17","18-64","65+"))
  ) %>%
  filter(is.finite(brr), brr > 0) %>%
  transmute(Region, setting, VE_label, AgeCat, outcome, brr)

brr_max <- quantile(brr_long_ob$brr, 0.999, na.rm = TRUE)  # 99.9%
brr_min <- quantile(brr_long_ob$brr, 0.001, na.rm = TRUE)  # 0.1%

lo_exp <- floor(log10(max(0.01, brr_min)))
hi_exp <- ceiling(log10(min(1e3, brr_max)))

thresholds_auto <- 10^seq(lo_exp, hi_exp, by = 0.02)
thresholds_auto

ceac_ob <- make_brr_ceac_outbreak(brr_long_ob,
                                  thresholds = thresholds_auto)


p_daly_ceac_outbreak <- plot_brr_ceac_outbreak_ve(ceac_ob, "DALY")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

p_death_ceac_outbreak <- plot_brr_ceac_outbreak_ve(ceac_ob, "Death")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "F") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

p_sae_ceac_outbreak <- plot_brr_ceac_outbreak_ve(ceac_ob, "SAE")+ 
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "G") +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1) 
  )

ggsave("06_Results/brrac_daly_ori.pdf", plot = p_daly_ceac_outbreak, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brrac_death_ori.pdf", plot = p_death_ceac_outbreak, width = 10, height = 8, device = cairo_pdf)
ggsave("06_Results/brrac_sae_ori.pdf", plot = p_sae_ceac_outbreak, width = 10, height = 8, device = cairo_pdf)


write.csv(ceac_ob, file = "06_Results/ceac_ob.csv")


ceac_t1 <- ceac_ob %>%
  filter(abs(log10(threshold)) < 1e-12) %>%   # == threshold=1
  dplyr::select(outcome, setting, AgeCat, p_accept, threshold, VE_label)

ceac_t1_wide <- ceac_t1 %>%
  mutate(p_fmt = sprintf("%.0f%%", 100 * p_accept)) %>%
  dplyr::select(outcome, setting, AgeCat, threshold, VE_label, p_fmt) %>%
  tidyr::pivot_wider(names_from = VE_label, values_from = p_fmt)

write_xlsx(ceac_t1_wide, "06_Results/ceac_travel.xlsx")
write_xlsx(ceac_df, "06_Results/ceac_travel_all.xlsx")
