
## Quantile function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.5, 0.975, 0.995), names = TRUE)
}


## Smooth trend model to evaluate difference between two time series
STM_diff <- function(data1, data2, gene_idx, condition_name, ps,
                     data_start, data_end){
  
  out_name <- paste0(condition_name, "_", names(data1)[gene_idx], "rep")
  
  # Load stan model
  model <- cmdstan_model("stan_model/SSM_diff_SunShade.stan")
  
  # Prepare data_list
  data_list <- list(
    N_time = length(unique(data1$Time)),
    N1 = as.numeric(table(data1$Time)),
    N2 = as.numeric(table(data2$Time)),
    sumN1 = sum(as.numeric(table(data1$Time))),
    sumN2 = sum(as.numeric(table(data2$Time))),
    Y1 = data1[,gene_idx],
    Y2 = data2[,gene_idx]
  )
  
  # SSM model
  fit <- model$sample(
    data = data_list,
    init = function() { list(mu1 = as.numeric(tapply(data1[,gene_idx], data1$Time, mean))) },
    seed = 10,
    iter_warmup = 1000,
    iter_sampling = 3000,
    thin = 3,
    chains = 4,
    parallel_chains = 4,
    max_treedepth = 15,
    adapt_delta = 0.99,
    refresh = 1000,
    show_messages = F,
    sig_figs = 4,
    output_dir = "/tmp",
    output_basename = out_name
  )
  
  # 99% Bayesian credible intervals
  outcsv_name <- list.files("/tmp")
  outcsv_name <- outcsv_name[grep(out_name, outcsv_name)]
  tmp_csv_mu1 <- NULL
  tmp_csv_mu2 <- NULL
  tmp_csv_diff <- NULL
  tmp_csv_s_mu1 <- NULL
  tmp_csv_s_diff <- NULL
  tmp_csv_s_Y <- NULL
  for(j in 1:length(outcsv_name)){
    tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[j])))
    tmp_csv_mu1 <- rbind(tmp_csv_mu1, tmp_csv[,str_starts(names(tmp_csv), "mu1")])
    tmp_csv_mu2 <- rbind(tmp_csv_mu2, tmp_csv[,str_starts(names(tmp_csv), "mu2")])
    tmp_csv_diff <- rbind(tmp_csv_diff, tmp_csv[,str_starts(names(tmp_csv), "diff\\.")])
    tmp_csv_s_mu1 <- c(tmp_csv_s_mu1, tmp_csv[,str_starts(names(tmp_csv), "s_mu1")])
    tmp_csv_s_diff <- c(tmp_csv_s_diff, tmp_csv[,str_starts(names(tmp_csv), "s_diff")])
    tmp_csv_s_Y <- c(tmp_csv_s_Y, tmp_csv[,str_starts(names(tmp_csv), "s_Y")])
  }
  
  # Calculate quantiles
  df_mu1 <- as.data.frame(t(apply(tmp_csv_mu1, 2, quantile99)))[data_start:data_end,]
  df_mu2 <- as.data.frame(t(apply(tmp_csv_mu2, 2, quantile99)))[data_start:data_end,]
  df_diff <- as.data.frame(t(apply(tmp_csv_diff, 2, quantile99)))[data_start:data_end,]
  df_s <- t(data.frame(s_mu1 = quantile99(tmp_csv_s_mu1),
                       s_diff = quantile99(tmp_csv_s_diff),
                       s_Y = quantile99(tmp_csv_s_Y)))
  df_s <- cbind(data.frame(s_name = row.names(df_s)), df_s)
  
  # Save output
  colnames(df_mu1) <- paste0("mu1_", colnames(df_mu1))
  colnames(df_mu2) <- paste0("mu2_", colnames(df_mu2))
  colnames(df_diff) <- paste0("diff_", colnames(df_diff))
  df <- cbind(data.frame(time = unique(data1$Time)[data_start:data_end]),
              data.frame(data1 = tapply(data1[,gene_idx], data1$Time, mean)[data_start:data_end],
                         data2 = tapply(data2[,gene_idx], data2$Time, mean)[data_start:data_end]),
              df_mu1, df_mu2, df_diff)
  fwrite(df, file = paste0(out, "STM_diff_", out_name, ".csv"))
  fwrite(df_s, file = paste0(out, "STM_diff_", out_name, "_sd.csv"))
  
  # Diagnosis of MCMC
  
  # Check of Rhat
  bayesplot::color_scheme_set("viridisC")
  bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = ps+2, base_family = "sans"))
  g <- bayesplot::mcmc_rhat(bayesplot::rhat(fit))
  ggsave(paste0(out, "STM_diff_rhat_", out_name, ".pdf"),
         g, height = ps*20, width = ps*20, units = "mm")
  max_rhat <- names(which.max(bayesplot::rhat(fit)))
  
  # Confirmation of convergence
  g <- bayesplot::mcmc_combo(
    fit$draws(),
    combo = c("dens_overlay", "trace"),
    widths = c(1, 1),
    pars = c("mu1[13]", "mu1[25]", "diff[13]", "diff[25]", "s_Y"),
    gg_theme = theme_classic(base_size = ps+2)
  )
  ggsave(paste0(out, "STM_diff_combo_", out_name, ".pdf"),
         g, height = ps*20, width = ps*20, units = "mm")
  
  # Remove temporary files
  file.remove(paste0("/tmp/", outcsv_name))
  
  return(fit)
}

