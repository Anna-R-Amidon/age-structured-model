library(bbmle)


age_data <- read.csv("WATS6460_Homework3_AgeStructuredData.csv",
                     header = TRUE, row.names = 1, check.names = FALSE)

harv_data <- read.csv("WATS6460_Homework3_HarvestRates.csv",
                      header = TRUE, row.names = 1, check.names = FALSE)


n_years <- ncol(age_data)   
n_start <- as.numeric(age_data$`Year 1`)                # Initial population vector 
harv_mat <- as.matrix(harv_data)                        # Harvest matrix

alpha <- 0.1                                            
beta <- 0.0001





# 1. Age-Structured Model
ages.mod <- function(alpha, beta, n_years, n_start, harvest_data){
  
  f1 <- 0
  f2 <- 10
  f3 <- 20
  f4 <- 50
  
  s1 <- 0.3
  s2 <- 0.5
  s3 <- 0.7
  s4 <- 0.8
  
  nage <- length(n_start)
  Nt <- matrix(0, nrow = nage, ncol = n_years)
  Nt[,1] <- n_start
  
  s0_vec <- rep(NA, n_years)
  harvest <- as.matrix(harvest_data)
  
  for(i in 2:n_years){
    
    
    s0 <- alpha * exp(-beta * sum(Nt[,i-1]))       # Density-dependent 
    s0_vec[i] <- s0
    
    leslie1 <- matrix(c(
      s0*f1, s0*f2, s0*f3, s0*f4,
      s1,    0,     0,     0,
      0,     s2,    0,     0,
      0,     0,     s3,    s4
    ), nrow = 4, byrow = TRUE)

    Nt_temp <- leslie1 %*% Nt[,i-1]
  
    Nt[,i] <- Nt_temp * (1 - harvest[,i-1])          # Apply harvest
    
    Nt[,i][Nt[,i] < 0] <- 0                          # Prevent negative values
  }
  
  return(list(Nt = Nt, s0 = s0_vec))
}

out <- ages.mod(alpha, beta, n_year = n_years, n_start = n_start, harvest_data = harv_mat)
str(out)




# 2. Negative Log-Likelihood
nll_fun <- function(log_alpha, log_beta, log_sigma, age_data, harv_data){
  
  alpha <- exp(log_alpha)
  beta  <- exp(log_beta)
  sigma <- exp(log_sigma)
  
  n_years <- ncol(age_data)
  n_start <- as.numeric(age_data$`Year 1`)
  harv_mat <- as.matrix(harv_data)
  
  out <- ages.mod(alpha = alpha,
                  beta = beta,
                  n_years = n_years,
                  n_start = n_start,
                  harvest_data = harv_mat)
  
  Npred <- out$Nt
  
  epsilon <- 1e-6                                              # Prevent log(0)
  Npred[Npred < epsilon] <- epsilon
  
  if(any(!is.finite(Npred))) return(1e12)
  
  Nobs <- as.matrix(age_data)
  

  nll <- -sum(dlnorm(Nobs,                                  # Log-normal likelihood
                     meanlog = log(Npred),
                     sdlog = sigma,
                     log = TRUE))
  
  return(nll)
}

nll_fun(
  log_alpha = log(0.1),
  log_beta  = log(1e-4),
  log_sigma = log(0.5),
  age_data  = age_data,
  harv_data = harv_data
)





# 3. Fit Model 
fit_sann <- mle2(
  nll_fun,
  start = list(log_alpha = log(0.5),
               log_beta  = log(0.0001),
               log_sigma = log(1)),           # large variance will catch more data points
  data = list(age_data = age_data,
              harv_data = harv_data),
  method = "SANN",
  control = list(maxit = 10000)
)

summary(fit_sann)




# 4. Parameter Estimates
coef_est <- coef(fit_sann)

alpha_hat <- exp(coef_est["log_alpha"])
beta_hat  <- exp(coef_est["log_beta"])
sigma_hat <- exp(coef_est["log_sigma"])

alpha_hat
beta_hat





# 5. 95% Profile Likelihood CIs

prof <- profile(fit_sann)
CI_log <- confint(prof)  
CI <- exp(CI_log)
CI





# 6. Plot Model Fit
pred <- ages.mod(alpha = alpha_hat,
                 beta  = beta_hat,
                 n_years = n_years,
                 n_start = n_start,
                 harvest_data = harv_mat)

Npred <- pred$Nt
Nobs  <- as.matrix(age_data)






# 7. 95% prediction intervals for log-normal error
lwrN <- qlnorm(0.025, meanlog = log(Npred), sdlog = sigma_hat)
uprN <- qlnorm(0.975, meanlog = log(Npred), sdlog = sigma_hat)


png("sablefish_model_fit.png", width = 800, height = 600)

# Plot model fit with uncertainty
matplot(t(Npred),
        type = "n",       
        xlab = "Year",
        ylab = "Abundance",
        ylim = c(0, max(uprN)*1.1))  

colors <- c("sienna2", "salmon4","paleturquoise3","peru")
for(age in 1:4){
  polygon(
    x = c(1:n_years, rev(1:n_years)),
    y = c(lwrN[age,], rev(uprN[age,])),
    col = adjustcolor(colors[age], 0.25),
    border = NA
  )
}

# predicted lines
matlines(t(Npred),
         lwd = 2,
         lty = 1,
         col = colors)

# observed points
matpoints(t(Nobs),
          pch = 16,
          col = colors)

legend("topright",
       legend = paste("Age", 1:4),
       col = colors,
       lwd = 2,
       pch = 16)

dev.off()


# Plot density-dependent age-0 survival

# Create a range of total population sizes
Nt_seq <- seq(0, max(colSums(Nobs)) * 1.1, length.out = 200)

# Density-dependent survival curve
s0_curve <- alpha_hat * exp(-beta_hat * Nt_seq)


png("density_dependent_survival.png", width = 700, height = 550)

plot(Nt_seq, s0_curve,
     type = "l",
     lwd = 3,
     col = "black",
     xlab = "Total population size (Nt)",
     ylab = "Age-0 survival probability (s0)",
     main = "Density-dependent age-0 survival")

abline(h = alpha_hat, lty = 2)

text(x = max(Nt_seq) * 0.7,
     y = alpha_hat * 0.95,
     labels = expression(paste("Baseline survival (", alpha, ")")),
     pos = 3)

dev.off()

