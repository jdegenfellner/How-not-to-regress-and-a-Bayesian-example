# Bayes 2

# change working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(pacman)
p_load(NHANES, tidyverse,  dagitty, ggdag, 
       car, rethinking, data.table, performance)
?NHANES

colnames(NHANES)
dim(NHANES) # 10000  86

sample_size_for_analysis <- 50

# Q: Total effect of Physical Activity on Systolic Blood Pressure---------------

# Distribution of blood pressure in the whole sample (adults and children)------
hist(NHANES$BPSysAve, breaks = 30, 
     main = "Distribution of Systolic Blood Pressure in NHANES (N=10000)", 
     xlab = "Systolic Blood Pressure (BPSysAve) mmHg", col = "lightblue")
summary(NHANES$BPSysAve)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 76.0   106.0   116.0   118.2   127.0   226.0    1449 
sd(NHANES$BPSysAve, na.rm = TRUE)
# 17.24817

# Define DAG--------------------------------------------------------------------
dag <- dagitty('dag {
  PhysActive -> BPSysAve
  Age -> PhysActive
  Age -> BPSysAve
  PhysActive -> BMI
  BMI -> BPSysAve
  Gender -> PhysActive
  Gender -> BMI
  Gender -> BPSysAve
}')

dagitty::coordinates(dag) <- list(
  x = c(PhysActive = 0, BPSysAve = 2, 
        Age = 1, BMI = 1, Gender = 0.5),
  y = c(PhysActive = 1, BPSysAve = 1, 
        Age = 2, BMI = 1.5, Gender = 2)
)

ggdag(dag) + 
  theme_minimal() + 
  geom_dag_point(size = 20, color = "black") +  
  geom_dag_text(size = 2.5, color = "white") +
  ggtitle("Hypothesized (Causal) Relationships") + 
  theme(plot.title = element_text(hjust = 0.5))

# Which covariates to put into the regression model?----------------------------

adjustmentSets(dag, exposure = "PhysActive", 
                    outcome = "BPSysAve")
# { Age, Gender }

impliedConditionalIndependencies(dag)

# Age _||_ BMI | Gndr, PhyA
# Age _||_ Gndr

# ---



# Check implied conditional independence:---------------------------------------

# Age _||_ BMI | Gnder, PhyA
mod <-lm(Age ~ BMI + Gender + PhysActive, 
           data = NHANES %>% dplyr::filter(Age >=20))
summary(mod)
coef(mod) # 0.04996778
confint(mod)
# BMI           -0.008479335  0.1084149

df_age <- NHANES %>% dplyr::filter(Age >=20)
df_age %>%
  ggplot(aes(x = BPSysAve)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, fill = "lightblue", alpha = 0.6) +  
  geom_density(color = "blue", linewidth = 1) +  
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(df_age$BPSysAve, na.rm = TRUE), 
                sd = sd(df_age$BPSysAve, na.rm = TRUE)), 
    color = "red", linewidth = 1, linetype = "dashed"
  ) +  # Theoretical normal curve
  labs(
    x = "Systolic Blood Pressure (BPSysAve)", 
    y = "Density", 
    title = "Distribution of Systolic Blood Pressure with Normal Curve"
  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))
# -> good enough for now?

# "Outlier": One could start with an analysis by excluding the large value:-----
# while being fully aware that this is a correct value!
df_age <- df_age %>% dplyr::filter(BPSysAve < 180) %>%
  sample_n(sample_size_for_analysis)


# PPC--------------------
# Prior predictive checks
set.seed(123)
n_sims <- dim(df_age)[1] 
beta_0_vec <- rnorm(n_sims, 125, 15)
beta_1_vec <- rnorm(n_sims, 0, 17/2) # ~ half a SD of BPSysAve in the whole population (one would look to a range of non-pharmacological blood pressure lowering interventions)? -> if DAG correct, this measures the causal effect (in an ideal world)
beta_2_vec <- rnorm(n_sims, 0, 5) # per year, not soo much change; literature ...
beta_3_vec <- rnorm(n_sims, 0, 10)
sigma_vec <- runif(n_sims, 0, 15) 
BPSysAve_sim <- rnorm(n_sims, beta_0_vec + beta_1_vec + beta_2_vec + beta_3_vec, sigma_vec)

dens <- density(BPSysAve_sim, na.rm = TRUE)
dens_obs <- density(df_age$BPSysAve, na.rm = TRUE)
hist(BPSysAve_sim, 
     #breaks = 30, 
     main = "Prior Predictive Distribution of BPSysAve", 
     xlab = "Simulated BPSysAve", col = "lightblue",
     ylim = c(0, max(range(dens$y)[2], range(dens_obs$y)[2])),
     freq = FALSE)
lines(dens, lwd = 3)
polygon(dens, 
        col = adjustcolor("darkblue", alpha.f = 0.2), 
        border = NA)
summary(BPSysAve_sim)
lines(dens_obs, lwd = 3, col = "red")
# -> not so bad a priori.
# One could argue that there are hardly people with <80 in the population.
sum(NHANES$BPSysAve < 80, na.rm = TRUE)/length(NHANES$BPSysAve)*100 


# SEED/filter/prepare-------------
set.seed(125)

df <- NHANES %>%
  dplyr::filter(Age >= 20) %>%
  dplyr::sample_n(50)

df <- df %>%
  dplyr::filter(BPSysAve < 180)

df <- df %>%
  dplyr::select(BPSysAve, PhysActive, Age, Gender, BMI) %>%
  dplyr::mutate(
    PhysActive = dplyr::case_when(
      PhysActive == "Yes" ~ 2L,
      PhysActive == "No"  ~ 1L,
      TRUE ~ NA_integer_
    ),
    PhysActive = as.integer(PhysActive)
  ) %>%
  tidyr::drop_na()

Age_mean <- mean(df$Age, na.rm = TRUE)

dim(df)
Age_mean


# Fit model--------------
df <- df %>% 
  dplyr::select(BPSysAve, PhysActive, Age, Gender, BMI)
dim(df) # 48 5

dim(df) #
Age_mean
m_NHANES <- quap(
  alist(
    BPSysAve ~ dnorm(mu, sigma), # might be necessary to adapt shape of likelihood and homoscedascity-assumption
    
    # likelihood
    mu <- beta_0 + beta_1[PhysActive] + beta_2 * (Age - Age_mean) + beta_3[Gender],
    
    beta_0 ~ dnorm(125, 15),
    
    # Prior 1:
    #beta_1[PhysActive] ~ dnorm(0, 10),  # in both phyActive Yes/No, normal
    
    # OR
    
    # Prior 2:
    beta_1[1] ~ dnorm(0, 10),
    beta_1[2] ~ dnorm(-20, 10), # = 10*MCID (2 mmHg)
    
    # OR
    
    # Prior 3:
    #beta_1[1] ~ dnorm(0, 5),
    #beta_1[2] ~ dnorm(-20, 5), # = 10*MCID (2 mmHg)
    
    beta_2 ~ dnorm(0, 17/2),
    beta_3[Gender] ~ dnorm(0, 10),  
    sigma ~ dunif(0, 15) 
  ),
  data = df
)
precis(m_NHANES, depth = 2)

post <- extract.samples(m_NHANES)
post$diff_physactive <- post$beta_1[, 2] - post$beta_1[, 1] # 2=Yes,1=No

dens(post$diff_physactive, 
     main = "Posterior Distribution of the Effect of Physical Activity (Yes/No) on BPSysAve", 
     xlab = "Difference in BPSysAve", 
     col = "lightblue",
     lwd = 4,
     show.HPDI = 0.89,
     show.zero = TRUE
     )
post_mean <- mean(post$diff_physactive)
abline(v = post_mean, col = "red", lwd = 3)
text(post_mean, par("usr")[4]*0.9,
     labels = paste0("Mean = ", round(post_mean, 3)),
     col = "red", pos = 4)
rope <- c(-2, 2)
points(rope, c(0, 0), 
       pch = 16, 
       col = "darkgreen", 
       cex = 1.6)
segments(x0 = rope[1], y0 = 0,
         x1 = rope[2], y1 = 0,
         col = "darkgreen",
         lwd = 2)
text(rope, c(0, 0),
     labels = c("-2", "+2"),
     col = "darkgreen",
     pos = 3)

precis(post, depth = 2)
# -> using a ROPE of +/-2mmHG -> we accept the null effect.


# -> Note that this value is -0.78 (-1.3460111  -0.1354783)
#    in the large sample using an improved model!


#________________------------------
# VARIABILITY-------------------------------------------------------------------

plot_prior_diff <- function(mu1, sd1,
                            mu2, sd2,
                            label = "Prior",
                            n = 10000,
                            path = "./NHANES_Results/0_Priors/",
                            width = 7, height = 5, dpi = 300,
                            make_plot = TRUE, save_plot = TRUE) {
  
  if (save_plot) dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  # Prior-Samples ziehen
  beta_1_1 <- rnorm(n, mu1, sd1)
  beta_1_2 <- rnorm(n, mu2, sd2)
  
  diff_prior <- beta_1_2 - beta_1_1
  
  draw_plot <- function() {
    hist(diff_prior, breaks = 40, freq = FALSE,
         main = paste0(label, ": Prior for Effect of PhysActive on BPSysAve"),
         xlab = expression(beta[Yes] - beta[No]),
         col = "lightblue", border = "white",
         cex.main = 0.85)
    
    dens <- density(diff_prior, na.rm = TRUE)
    lines(dens, lwd = 3)
    
    abline(v = mean(diff_prior), col = "red", lwd = 3)
    text(mean(diff_prior), par("usr")[4]*0.9,
         labels = paste0("Mean = ", round(mean(diff_prior),2)),
         col = "red", pos = 4)
  }
  
  if (make_plot) draw_plot()
  
  if (save_plot) {
    file_name <- paste0(path, gsub(" ", "_", label), ".png")
    
    png(file_name, width = width, height = height, units = "in", res = dpi)
    draw_plot()
    dev.off()
    
    cat("Saved:", file_name, "\n")
  }
  
  invisible(diff_prior)
}

fit_quap_prior1 <- function(df, Age_mean) {
  quap(
    alist(
      BPSysAve ~ dnorm(mu, sigma),
      mu <- beta_0 + beta_1[PhysActive] + beta_2 * (Age - Age_mean) + beta_3[Gender],
      beta_0 ~ dnorm(125, 15),
      
      # Prior 1:
      beta_1[PhysActive] ~ dnorm(0, 10),
      
      beta_2 ~ dnorm(0, 17/2),
      beta_3[Gender] ~ dnorm(0, 10),
      sigma ~ dunif(0, 15)
    ),
    data = df
  )
}

fit_quap_prior2 <- function(df, Age_mean) {
  quap(
    alist(
      BPSysAve ~ dnorm(mu, sigma),
      mu <- beta_0 + beta_1[PhysActive] + beta_2 * (Age - Age_mean) + beta_3[Gender],
      beta_0 ~ dnorm(125, 15),
      
      # Prior 2:
      beta_1[1] ~ dnorm(0, 10),
      beta_1[2] ~ dnorm(-20, 10),
      
      beta_2 ~ dnorm(0, 17/2),
      beta_3[Gender] ~ dnorm(0, 10),
      sigma ~ dunif(0, 15)
    ),
    data = df
  )
}

fit_quap_prior3 <- function(df, Age_mean) {
  quap(
    alist(
      BPSysAve ~ dnorm(mu, sigma),
      mu <- beta_0 + beta_1[PhysActive] + beta_2 * (Age - Age_mean) + beta_3[Gender],
      beta_0 ~ dnorm(125, 15),
      
      # Prior 3:
      beta_1[1] ~ dnorm(0, 5),
      beta_1[2] ~ dnorm(-20, 5),
      
      beta_2 ~ dnorm(0, 17/2),
      beta_3[Gender] ~ dnorm(0, 10),
      sigma ~ dunif(0, 15)
    ),
    data = df
  )
}

fit_quap_prior4 <- function(df, Age_mean) {
  quap(
    alist(
      BPSysAve ~ dnorm(mu, sigma),
      mu <- beta_0 + beta_1[PhysActive] + beta_2 * (Age - Age_mean) + beta_3[Gender],
      beta_0 ~ dnorm(125, 15),
      
      # Prior 4:
      beta_1[1] ~ dnorm(0, 5),
      beta_1[2] ~ dnorm(-5, 5),
      
      beta_2 ~ dnorm(0, 17/2),
      beta_3[Gender] ~ dnorm(0, 10),
      sigma ~ dunif(0, 15)
    ),
    data = df
  )
}

fit_quap_prior5 <- function(df, Age_mean) {
  quap(
    alist(
      BPSysAve ~ dnorm(mu, sigma),
      mu <- beta_0 + beta_1[PhysActive] + beta_2 * (Age - Age_mean) + beta_3[Gender],
      beta_0 ~ dnorm(125, 15),
      
      # Prior 3:
      beta_1[1] ~ dnorm(0, 0.5),
      beta_1[2] ~ dnorm(-1, 0.5),
      
      beta_2 ~ dnorm(0, 17/2),
      beta_3[Gender] ~ dnorm(0, 10),
      sigma ~ dunif(0, 15)
    ),
    data = df
  )
}

#set to parent directory
setwd(dirname(getwd()))
getwd()

# PRIOR 1  (unspecific)
plot_prior_diff(mu1 = 0,  sd1 = 10,
                mu2 = 0,  sd2 = 10,
                label = "Prior_1_Weak_Symmetric")

# PRIOR 2  (expects strong negative effect)
plot_prior_diff(mu1 = 0,   sd1 = 10,
                mu2 = -20, sd2 = 10,
                label = "Prior_2_Strong_Negative")

# PRIOR 3  (same mean but tighter belief)
plot_prior_diff(mu1 = 0,   sd1 = 5,
                mu2 = -20, sd2 = 5,
                label = "Prior_3_Strong_Tight")

# PRIOR 4  (moderate negative belief)
plot_prior_diff(mu1 = 0,  sd1 = 5,
                mu2 = -5, sd2 = 5,
                label = "Prior_4_Moderate")

# PRIOR 5  (very informative, close to expected real effect)
plot_prior_diff(mu1 = 0,  sd1 = 0.5,
                mu2 = -1, sd2 = 0.5,
                label = "Prior_5_Very_Informative")

# Main function: sample + fit + extract diff
run_nhanes_physactive <- function(seed,
                                  prior = c(1, 2, 3, 4, 5),
                                  n = 50,
                                  age_min = 20,
                                  drop_bp_over = 180,
                                  rope = c(-2, 2),
                                  hpdi_prob = 0.89,
                                  verbose = TRUE,
                                  make_plot = TRUE,
                                  save_plot = TRUE,
                                  width = 7,
                                  height = 5,
                                  dpi = 300) {
  
  prior <- as.integer(prior[1])
  if (!prior %in% c(1, 2, 3, 4, 5)) stop("prior must be 1, 2, 3, 4 or 5")
  
  set.seed(seed)
  
  df0 <- NHANES %>%
    dplyr::filter(Age >= age_min) %>%
    dplyr::sample_n(n)
  
  if (!is.null(drop_bp_over)) {
    df0 <- df0 %>% dplyr::filter(BPSysAve < drop_bp_over)
  }
  
  df0 <- df0 %>%
    dplyr::select(BPSysAve, PhysActive, Age, Gender, BMI) %>%
    dplyr::mutate(
      PhysActive = dplyr::case_when(
        PhysActive == "Yes" ~ 2L,
        PhysActive == "No"  ~ 1L,
        TRUE ~ NA_integer_
      )
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      Gender = factor(Gender, levels = c("female", "male")),
      PhysActive = factor(PhysActive, levels = c(1, 2))
    )
  
  Age_mean <- mean(df0$Age, na.rm = TRUE)
  df0 <- as.data.frame(df0)
  
  # ---- error catcher: if quap (or prior function) fails -> return NULL + no save 
  m <- tryCatch(
    {
      switch(
        as.character(prior),
        "1" = fit_quap_prior1(df0, Age_mean),
        "2" = fit_quap_prior2(df0, Age_mean),
        "3" = fit_quap_prior3(df0, Age_mean),
        "4" = fit_quap_prior4(df0, Age_mean),
        "5" = fit_quap_prior5(df0, Age_mean)
      )
    },
    error = function(e) {
      if (verbose) {
        message("quap failed for seed=", seed, " prior=", prior)
        message("Reason: ", e$message)
      }
      return(NULL)
    }
  )
  
  if (is.null(m)) return(NULL)
  
  post <- extract.samples(m)
  diff_physactive <- post$beta_1[, 2] - post$beta_1[, 1]
  
  draw_plot <- function() {
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    dens(
      diff_physactive,
      main = paste0(
        "Posterior: Effect of Physical Activity (Yes/No) on BPSysAve\n",
        "seed=", seed, ", prior=", prior, ", n=", nrow(df0)
      ),
      xlab = "Difference in BPSysAve mmHg (89% CredibleInt)",
      col = "lightblue",
      lwd = 4,
      show.HPDI = hpdi_prob,
      show.zero = TRUE
    )
    
    post_mean <- mean(diff_physactive, na.rm = TRUE)
    abline(v = post_mean, col = "red", lwd = 3)
    
    text(post_mean, par("usr")[4] * 0.9,
         labels = paste0("Mean = ", round(post_mean, 3)),
         col = "red", pos = 4)
    
    rope_sorted <- sort(as.numeric(rope))
    
    segments(rope_sorted[1], 0, rope_sorted[2], 0,
             col = "darkgreen", lwd = 2)
    
    points(rope_sorted, c(0, 0),
           pch = 16, col = "darkgreen", cex = 1.6)
    
    y_text <- par("usr")[3] + 0.04 * diff(par("usr")[3:4])
    
    text(rope_sorted, c(y_text, y_text),
         labels = c(paste0(rope_sorted[1]), paste0("+", rope_sorted[2])),
         col = "darkgreen", pos = 3)
  }
  
  if (make_plot) draw_plot()
  
  file_name <- NULL
  if (save_plot) {
    dir.create("./NHANES_Results", showWarnings = FALSE)
    
    file_name <- paste0(
      "./NHANES_Results/",
      "Posterior_seed", seed,
      "_prior", prior,
      "_n", nrow(df0),
      ".png"
    )
    
    png(file_name, width = width, height = height, units = "in", res = dpi)
    draw_plot()
    dev.off()
  }
  
  out <- list(
    seed = seed,
    prior = prior,
    n_used = nrow(df0),
    Age_mean = Age_mean,
    model = m,
    diff_physactive = diff_physactive,
    file_name = file_name
  )
  
  if (verbose) {
    if (!is.null(file_name)) cat("\nSaved plot to:", file_name, "\n")
    print(precis(m, depth = 2))
  }
  
  return(out)
}


#RUN------

for(seed_for_plot in 21:40){
  run_nhanes_physactive(seed = seed_for_plot, prior = 1, n = 50, # rather unsure.
                        verbose = TRUE, make_plot = TRUE)
  run_nhanes_physactive(seed = seed_for_plot, prior = 2, n = 50, # not so sure, strong negative effect
                        verbose = TRUE, make_plot = TRUE)
  run_nhanes_physactive(seed = seed_for_plot, prior = 3, n = 50, # rather sure, strong negative effect
                        verbose = TRUE, make_plot = TRUE)
  run_nhanes_physactive(seed = seed_for_plot, prior = 4, n = 50, # rather sure and negative effect
                        verbose = TRUE, make_plot = TRUE)  
  run_nhanes_physactive(seed = seed_for_plot, prior = 5, n = 50, # rather sure, effect size similar to large sample
                        verbose = TRUE, make_plot = TRUE)  
}


summarise_run <- function(res, hpdi_prob = 0.89, rope = c(-2, 2)) {
  if (is.null(res)) return(NULL)
  
  d <- res$diff_physactive
  hp <- HPDI(d, prob = hpdi_prob)
  
  data.frame(
    seed = res$seed,
    prior = res$prior,
    n_used = res$n_used,
    mean = mean(d),
    sd = sd(d),
    hpdi_low = hp[1],
    hpdi_high = hp[2],
    rope_low = rope[1],
    rope_high = rope[2],
    rope_prop = mean(d >= min(rope) & d <= max(rope)),
    stringsAsFactors = FALSE
  )
}

results <- list()

for (seed_for_plot in 21:40) {
  for (p in 1:5) {
    results[[length(results) + 1]] <- run_nhanes_physactive(
      seed = seed_for_plot, prior = p, n = 50,
      verbose = FALSE, make_plot = FALSE, save_plot = FALSE
    )
  }
}

tab <- dplyr::bind_rows(lapply(results, summarise_run, hpdi_prob = 0.89, rope = c(-2, 2))) %>%
  dplyr::arrange(prior, seed)

tab

knitr::kable(
  tab,
  digits = 3,
  caption = "Posterior summaries for Δ = beta_Yes − beta_No (BPSysAve)"
)

# plot:
tab_plot <- tab %>%
  mutate(prior = factor(prior))

ggplot(tab_plot, aes(x = prior, y = mean)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Posterior effect size (Δ = Yes − No) stratified by prior",
    x = "Prior",
    y = "Posterior mean of Δ [mmHg]"
  ) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = c(-2, 2), 
             color = "darkgreen",
             linewidth = 2)



# Compare with lm()------------
coef_vec_physAct <- numeric(15)
# HERE--------
for(i in 1:15) {
  set.seed(100 + i)
  df <- NHANES %>%
    dplyr::filter(Age >= 20) %>%
    dplyr::sample_n(50) %>%
    dplyr::filter(BPSysAve < 180)
  df$Age_centered <- df$Age - Age_mean
  lm_fit <- lm(BPSysAve ~ PhysActive + Age_centered + Gender,
               data = df)
  #check_model(lm_fit) 
  #check_model(lm_fit, check = "pp_check") # maybe improve a little...
  #qqPlot(lm_fit) # ok
  
  #plot(residuals(lm_fit) ~ fitted(lm_fit)) # not so bad
  #abline(h = 0, col = "red", lwd = 2)
  
  #summary(lm_fit)
  coef_vec_physAct[i] <- coef(lm_fit)[2] # 2.9543297
  #confint(lm_fit, level = 0.89) # -4.6437203  10.5523796
}
  
hist(coef_vec_physAct, breaks = 30, 
     main = "Distribution of Estimated Effect of PhysActive on BPSysAve (lm)", 
     xlab = "Estimated Coefficient for PhysActive (Yes vs. No)", 
     col = "lightblue")
boxplot(coef_vec_physAct, 
        main = "Boxplot of Estimated Effect of PhysActive on BPSysAve (lm)", 
        ylab = "Estimated Coefficient for PhysActive (Yes vs. No)", 
        col = "lightblue")


# _What if we add BMI?-------------
# BMI is on a pipe between PhysActive and BPSysAve 
lm_fit_bmi <- lm(BPSysAve ~ PhysActive + Age_centered + Gender + BMI,
                 data = df)
check_model(lm_fit_bmi)
check_model(lm_fit_bmi, check = "pp_check") # maybe improve a little
qqPlot(lm_fit_bmi) # seems worse in the middle
plot(residuals(lm_fit_bmi) ~ fitted(lm_fit_bmi)) # not so bad
abline(h = 0, col = "red", lwd = 2)
summary(lm_fit_bmi)
confint(lm_fit_bmi, level = 0.89) # -3.16194636  11.9645231
coef(lm_fit_bmi) # 4.4012884 vs. 2.9543297 (+50%) before without BMI
