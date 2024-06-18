
# Load packages
library(tidyverse)
library(data.table)
library(cmdstanr)
set_cmdstan_path("~/cmdstan/")

# Create output directory
out <- "03_LLMX_MarSep_BLRP/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}


### Preparation of data

# Load data of March and September
data <- read.csv("data/march_data.csv")
data_MarSun <- subset(data, Condition=="Sun")
data_MarShade <- subset(data, Condition=="Shade")
data <- read.csv("data/september_data.csv")
data_SepSun <- subset(data, Condition=="Sun")
data_SepShade <- subset(data, Condition=="Shade")
environment_data_march <- read.csv("data/environment_data_corrected_march.csv")
environment_data_march$Condition = factor(environment_data_march$Condition, levels=c("Sun", "Shade"))
environment_data_march <- environment_data_march %>%
  arrange(Condition, Time)
environment_data_september <- read.csv("data/environment_data_corrected_september.csv")
environment_data_september$Condition = factor(environment_data_september$Condition, levels=c("Sun", "Shade"))
environment_data_september <- environment_data_september %>%
  arrange(Condition, Time)
env_MarSun <- environment_data_march[environment_data_march$Condition=="Sun",]
env_MarShade <- environment_data_march[environment_data_march$Condition=="Shade",]
env_SepSun <- environment_data_september[environment_data_september$Condition=="Sun",]
env_SepShade <- environment_data_september[environment_data_september$Condition=="Shade",]

env_MarSun2 <- NULL; env_MarShade2 <- NULL
env_SepSun2 <- NULL; env_SepShade2 <- NULL
for(i in 1:length(unique(data_MarSun$Time))){
  env_MarSun2 <- rbind(env_MarSun2, env_MarSun[env_MarSun$Time == ts(unique(data_MarSun$Time))[i],])
  env_MarShade2 <- rbind(env_MarShade2, env_MarShade[env_MarShade$Time == ts(unique(data_MarShade$Time))[i],])
  env_SepSun2 <- rbind(env_SepSun2, env_SepSun[env_SepSun$Time == ts(unique(data_SepSun$Time))[i],])
  env_SepShade2 <- rbind(env_SepShade2, env_SepShade[env_SepShade$Time == ts(unique(data_SepShade$Time))[i],])
}

env_MarSun2 <- rbind(env_MarSun2, env_MarSun[nrow(env_MarSun),])
env_MarShade2 <- rbind(env_MarShade2, env_MarShade[nrow(env_MarShade),])
env_SepSun2 <- rbind(env_SepSun2, env_SepSun[nrow(env_SepSun),])
env_SepShade2 <- rbind(env_SepShade2, env_SepShade[nrow(env_SepShade),])
env_SepSun2 <- rbind(env_SepSun[1,], env_SepSun2)

# SIG5
SIG5_MarSun = data_MarSun %>%
  group_by(Time) %>%
  summarize(Mean_SIG5 = mean(SIG5))
SIG5_MarShade = data_MarShade %>%
  group_by(Time) %>%
  summarize(Mean_SIG5 = mean(SIG5))
SIG5_SepSun = data_SepSun %>%
  group_by(Time) %>%
  summarize(Mean_SIG5 = mean(SIG5))
SIG5_SepShade = data_SepShade %>%
  group_by(Time) %>%
  summarize(Mean_SIG5 = mean(SIG5))

# Lag of environmental effect
Lag_env <- expand.grid(list(
  Lag_temp=c(0),
  Lag_light=c(0),
  Lag_SIG5=c(0, 1, 2, 3, 4, 5)))
  

## Load stan model
model <- cmdstan_model("stan_model/LLMX_MarSep_BLRP.stan")


## Function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.5, 0.975, 0.995), names = TRUE)
}


### SSM for AhgBLRP
for(i in 1){
  
  max_Lag <- max(c(Lag_env$Lag_temp[i], Lag_env$Lag_light[i], Lag_env$Lag_SIG5[i]))
  N_MarSun = as.numeric(table(data_MarSun$Time))
  N_MarShade = as.numeric(table(data_MarShade$Time))
  N_SepSun = as.numeric(table(data_SepSun$Time))
  N_SepShade = as.numeric(table(data_SepShade$Time))
  
  data_list <- list(
    N_time = length(unique(data_MarSun$Time)),
    Lag_temp = Lag_env$Lag_temp[i],
    Lag_light = Lag_env$Lag_light[i],
    Lag_SIG5 = Lag_env$Lag_SIG5[i],
    max_Lag = max_Lag,
    N_MarSun = N_MarSun,
    N_MarShade = N_MarShade,
    N_SepSun = N_SepSun,
    N_SepShade = N_SepShade,
    sumN_MarSun = sum(N_MarSun),
    sumN_MarShade = sum(N_MarShade),
    sumN_SepSun = sum(N_SepSun),
    sumN_SepShade = sum(N_SepShade),
    temp_MarSun = env_MarSun2$Temperature,
    temp_MarShade = env_MarShade2$Temperature,
    temp_SepSun = env_SepSun2$Temperature,
    temp_SepShade = env_SepShade2$Temperature,
    light_MarSun = env_MarSun2$Irradiance,
    light_MarShade = env_MarShade2$Irradiance,
    light_SepSun = env_SepSun2$Irradiance,
    light_SepShade = env_SepShade2$Irradiance,
    SIG5_MarSun = SIG5_MarSun$Mean_SIG5,
    SIG5_MarShade = SIG5_MarShade$Mean_SIG5,
    SIG5_SepSun = SIG5_SepSun$Mean_SIG5,
    SIG5_SepShade = SIG5_SepShade$Mean_SIG5,
    Y_MarSun = data_MarSun$BLRP,
    Y_MarShade = data_MarShade$BLRP,
    Y_SepSun = data_SepSun$BLRP,
    Y_SepShade = data_SepShade$BLRP
    )
  
  out_name <- paste0("LLMX_MarSep_BLRP_", "Lagtemp", data_list$Lag_temp, "_Laglight", data_list$Lag_light, "_LagSIG5", data_list$Lag_SIG5)
  ps <- 7
  
  # SSM model
  fit <- model$sample(
    data = data_list,
    init = function() { list(alpha_MarSun = as.numeric(tapply(data_MarSun$BLRP, data_MarSun$Time, mean, na.rm=T)),
                             alpha_MarShade = as.numeric(tapply(data_MarShade$BLRP, data_MarShade$Time, mean, na.rm=T)),
                             alpha_SepSun = as.numeric(tapply(data_SepSun$BLRP, data_SepSun$Time, mean, na.rm=T)),
                             alpha_SepShade = as.numeric(tapply(data_SepShade$BLRP, data_SepShade$Time, mean, na.rm=T))) },
    seed = 10,
    iter_warmup = 3000,
    iter_sampling = 1000,
    #thin = 3,
    chains = 4,
    parallel_chains = 4,
    max_treedepth = 15,
    adapt_delta = 0.99,
    refresh = 1000,
    show_messages = F,
    sig_figs = 4,
    output_dir = "/tmp",
    output_basename = paste0("BLRP_LLMX", i)
  )
  
  # 99% Bayesian credible intervals
  outcsv_name <- list.files("/tmp")
  outcsv_name <- outcsv_name[grep(paste0("BLRP_LLMX", i), outcsv_name)]
  tmp_csv_mu <- NULL
  tmp_csv_b_temp <- NULL
  tmp_csv_b_light <- NULL
  tmp_csv_b_SIG5 <- NULL
  tmp_csv_s_mu <- NULL
  tmp_csv_s_Y <- NULL
  tmp_csv_alpha_MarSun <- NULL
  tmp_csv_alpha_MarShade <- NULL
  tmp_csv_alpha_SepSun <- NULL
  tmp_csv_alpha_SepShade <- NULL
  tmp_csv_log_lik_MarSun <- NULL
  tmp_csv_log_lik_MarShade <- NULL
  tmp_csv_log_lik_SepSun <- NULL
  tmp_csv_log_lik_SepShade <- NULL
  for(j in 1:length(outcsv_name)){
    tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[j])))
    tmp_csv_mu <- rbind(tmp_csv_mu, tmp_csv[,str_starts(names(tmp_csv), "mu\\.")])
    tmp_csv_b_temp <- c(tmp_csv_b_temp, tmp_csv[,str_starts(names(tmp_csv), "b_temp")])
    tmp_csv_b_light <- c(tmp_csv_b_light, tmp_csv[,str_starts(names(tmp_csv), "b_light")])
    tmp_csv_b_SIG5 <- c(tmp_csv_b_SIG5, tmp_csv[,str_starts(names(tmp_csv), "b_SIG5")])
    tmp_csv_s_mu <- c(tmp_csv_s_mu, tmp_csv[,str_starts(names(tmp_csv), "s_mu")])
    tmp_csv_s_Y <- c(tmp_csv_s_Y, tmp_csv[,str_starts(names(tmp_csv), "s_Y")])
    tmp_csv_alpha_MarSun <- rbind(tmp_csv_alpha_MarSun, tmp_csv[,str_starts(names(tmp_csv), "alpha_MarSun")])
    tmp_csv_alpha_MarShade <- rbind(tmp_csv_alpha_MarShade, tmp_csv[,str_starts(names(tmp_csv), "alpha_MarShade")])
    tmp_csv_alpha_SepSun <- rbind(tmp_csv_alpha_SepSun, tmp_csv[,str_starts(names(tmp_csv), "alpha_SepSun")])
    tmp_csv_alpha_SepShade <- rbind(tmp_csv_alpha_SepShade, tmp_csv[,str_starts(names(tmp_csv), "alpha_SepShade")])
    tmp_csv_log_lik_MarSun <- rbind(tmp_csv_log_lik_MarSun, tmp_csv[,str_starts(names(tmp_csv), "log_lik_MarSun")])
    tmp_csv_log_lik_MarShade <- rbind(tmp_csv_log_lik_MarShade, tmp_csv[,str_starts(names(tmp_csv), "log_lik_MarShade")])
    tmp_csv_log_lik_SepSun <- rbind(tmp_csv_log_lik_SepSun, tmp_csv[,str_starts(names(tmp_csv), "log_lik_SepSun")])
    tmp_csv_log_lik_SepShade <- rbind(tmp_csv_log_lik_SepShade, tmp_csv[,str_starts(names(tmp_csv), "log_lik_SepShade")])
  }
  
  df_mu <- as.data.frame(round(t(apply(tmp_csv_mu, 2, quantile99)), digits = 4))
  df_b_temp <- as.data.frame(round(t(quantile99(tmp_csv_b_temp)), digits = 4))
  row.names(df_b_temp) <- "b_temp"
  df_b_light <- as.data.frame(round(t(quantile99(tmp_csv_b_light)), digits = 4))
  row.names(df_b_light) <- "b_light"
  df_b_SIG5 <- as.data.frame(round(t(quantile99(tmp_csv_b_SIG5)), digits = 4))
  row.names(df_b_SIG5) <- "b_SIG5"
  df_s_mu <- as.data.frame(round(t(quantile99(tmp_csv_s_mu)), digits = 4))
  row.names(df_s_mu) <- "s_mu"
  df_s_Y <- as.data.frame(round(t(quantile99(tmp_csv_s_Y)), digits = 4))
  row.names(df_s_Y) <- "s_Y"
  df_alpha_MarSun <- as.data.frame(round(t(apply(tmp_csv_alpha_MarSun, 2, quantile99)), digits = 4))
  df_alpha_MarShade <- as.data.frame(round(t(apply(tmp_csv_alpha_MarShade, 2, quantile99)), digits = 4))
  df_alpha_SepSun <- as.data.frame(round(t(apply(tmp_csv_alpha_SepSun, 2, quantile99)), digits = 4))
  df_alpha_SepShade <- as.data.frame(round(t(apply(tmp_csv_alpha_SepShade, 2, quantile99)), digits = 4))
  
  df_all <- rbind(df_mu, df_b_temp, df_b_light, df_b_SIG5, df_s_mu, df_s_Y, df_alpha_MarSun, df_alpha_MarShade, df_alpha_SepSun, df_alpha_SepShade)
  df_all <- df_all %>%
    mutate(par = row.names(df_all)) %>%
    relocate(par)
  fwrite(df_all, file = paste0(out, out_name, ".csv"))
  
  # MCMC samples
  tmp_csv_all <- cbind(tmp_csv_mu, tmp_csv_b_temp, tmp_csv_b_light, tmp_csv_b_SIG5, tmp_csv_s_mu, tmp_csv_s_Y, tmp_csv_alpha_MarSun, tmp_csv_alpha_MarShade, tmp_csv_alpha_SepSun, tmp_csv_alpha_SepShade)
  names(tmp_csv_all) <- str_replace(names(tmp_csv_all), "tmp_csv_", "")
  fwrite(tmp_csv_all, file = paste0(out, out_name, "_MCMC.csv"))
  
  # log likelihood
  df_log_lik <- as.data.frame(cbind(tmp_csv_log_lik_MarSun, tmp_csv_log_lik_MarShade, 
                                    tmp_csv_log_lik_SepSun, tmp_csv_log_lik_SepShade))
  fwrite(df_log_lik, file = paste0(out, out_name, "_loglik.csv"))

  
  ## Diagnosis of MCMC
  
  # Check of Rhat
  # Remove NAs
  rhat_fit <- bayesplot::rhat(fit)
  rhat_fit2 <- rhat_fit[!str_detect(names(rhat_fit), paste0("\\[", data_list$max_Lag,"\\]"))]
  rhat_fit3 <- rhat_fit2[!is.na(rhat_fit2)]
  
  bayesplot::color_scheme_set("viridisC")
  bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = ps+2, base_family = "sans"))
  g <- bayesplot::mcmc_rhat(rhat_fit3)
  ggsave(paste0(out, out_name, "_rhat.pdf"),
         g, height = ps*20, width = ps*20, units = "mm")
  max_rhat <- names(which.max(rhat_fit3))
  
  # Get NA position
  fit_draws <- fit$draws()
  na_pos <- NULL
  for(l in 1:dim(fit_draws)[3]){
    if(sum(is.na(fit_draws[,,l])) != 0){
      na_pos <- c(na_pos, l)
    }
  }
  
  if(is.null(na_pos)){
    fit_draws2 <- fit_draws
  }else{
    fit_draws2 <- fit_draws[,,-na_pos]
  }
  
  # Confirmation of convergence
  g <- bayesplot::mcmc_combo(
    fit_draws2,
    combo = c("dens_overlay", "trace"),
    widths = c(1, 1),
    pars = c(paste0("mu[", data_list$max_Lag+1,"]"), paste0("mu[", data_list$N_time, "]"), "b_temp", "b_light", "b_SIG5"),
    gg_theme = theme_classic(base_size = ps+2)
  )
  ggsave(paste0(out, out_name, "_combo.pdf"),
         g, height = ps*20, width = ps*20, units = "mm")
  
  # Remove temporary files
  file.remove(paste0("/tmp/", outcsv_name))

}

