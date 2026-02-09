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

names(preui_all) <- c("Ceará","Alagoas", "Bahia", "Goiás",
                      "Minas Gerais", "Paraíba", "Piauí",
                      "Pernambuco", "Rio Grande do Norte", 
                      "Sergipe", "Tocantins")
setting_key <- c(
  "Ceará"             = "10%+",
  "Bahia"             = "<1%",
  "Paraíba"           = "1-10%",
  "Pernambuco"        = "1-10%",
  "Rio Grande do Norte" = "1-10%",
  "Piauí"             = "10%+",
  "Tocantins"         = "1-10%",
  "Alagoas"           = "10%+",
  "Minas Gerais"      = "<1%",
  "Sergipe"           = "1-10%",
  "Goiás"             = "<1%"
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
calc_total_hosp_draws_rho <- function(pre_list, post_arr, hosp_rate, rho_vec) {
  
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling (true symptomatic)
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  # 3) hosp totals
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

calc_total_fatal_draws_rho <- function(pre_list, post_arr,
                                       hosp_rate, fatal_rate, nh_fatal_rate,
                                       rho_vec) {
  # 1) draw 길이 맞추기
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  # 3) 기존 함수 그대로 사용
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
  
  # 1) draw 길이 맞추기
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  # 3) 기존 DALY 함수 그대로
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
        
        # 해당 조합에서 가능한 draw 수
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        # rho는 독립적으로 리샘플링 (500 rho pool -> 1000 draw도 가능)
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
  select(draw_id, total_pre, total_post, Region, VE, Coverage, Scenario)


# add setting key variable
all_draws_ix_true <- all_draws_ix_true %>%
  mutate(setting = unname(setting_key[Region]))

all_draws_hosp_true <- all_draws_hosp_true %>%
  mutate(setting = unname(setting_key[Region]))

all_draws_fatal_true <- all_draws_fatal_true %>%
  mutate(setting = unname(setting_key[Region]))

all_draws_daly_true <- all_draws_daly_true %>%
  mutate(setting = unname(setting_key[Region]))
