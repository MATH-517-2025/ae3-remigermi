
# Simulation study for plug-in h_AMISE using blockwise degree-4 OLS fits



set.seed(12345)



library(ggplot2)
library(dplyr)
library(tidyr)


# Regression function
m_fun <- function(x) {
  sin((x/3 + 0.1)^(-1))
}


# Estimate h_AMISE for a single simulated dataset
estimate_h <- function(n, N, alpha, beta, sigma2 = 1.0) {
  # we simulate the random variable
  X <- rbeta(n, alpha, beta)
  Y <- m_fun(X) + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  # we separate in N blocks
  block_idx <- rep(1:N, each = ceiling(n/N))[1:n]
  
  # --- Fit polynomial regression of degree 4 ---
  fitted_vals <- numeric(n)
  second_deriv <- numeric(n)
  for (j in seq_len(N)) {
    idx <- which(block_idx == j)
    Xj <- X[idx]; Yj <- Y[idx]
    fit <- lm(Yj ~ Xj + I(Xj^2) + I( Xj^3) + I( Xj^4))
    fitted_vals[idx] <- fitted.values(fit)
    #second derivative
    coefs <- coef(fit)
    second_deriv[idx] <- 2*coefs[2] + 6*coefs[3]*Xj + 12*coefs[4]*(Xj^2)
  }
  
  
  # theta22_hat
  theta22_hat <- mean(second_deriv^2)
  # sigma^2_hat: denominator n - 5N (fallback to n-1 if nonpositive)
  denom <- n - 5 * N
  if (denom <= 0) denom <- n - 1
  sigma2_hat <- sum((Y - fitted_vals)^2, na.rm = TRUE) / denom
  # RSS(N)
  RSS_hat <- denom*sigma2_hat
  
  # h_AMISE estimation :
  supp_len <- 1.0 # Beta on [0,1]
  if (theta22_hat <= 1e-12) { # if theta is too small, h_hat is not defined
    h_hat <- NA_real_
  } else {
    h_hat <- n^(-1/5) * ((35 * sigma2_hat * supp_len / theta22_hat)^(1/5))
  }
  list(h = h_hat, theta22 = theta22_hat, sigma2 = sigma2_hat, RSS = RSS_hat )
}

# Grid and simulation settings (you can adjust)
ns <- c(200,500, 1000, 2000, 5000)
betas_params <- list(c(2,2), c(0.5,0.5),c(1,3), c(3,1), c(5,1), c(1,5), c(5,5))


# Storage
results <- list()
selected_N <- data.frame()
# Main simulation loop
for (bb in betas_params) {
  alpha <- bb[1]; beta <- bb[2]
  for (n in ns) {
    Nmax <- max(min(floor(n / 20), 5), 1)
    N_best <-Nmax
    rss_Nmax <-estimate_h(n, Nmax, alpha, beta, sigma2 = 1.0)$RSS
    Cp_min <- rss_Nmax / (rss_Nmax / (n - 5 * Nmax)) - (n - 10 * Nmax)
    
    for (Nblocks in seq(Nmax)) {
      est <- estimate_h(n, Nblocks, alpha, beta, sigma2 = 1.0)
      
     # choose N_best who minimize Cp
      Cp_N <- est$RSS / (rss_Nmax / (n - 5 * Nmax)) - (n - 10 * Nblocks)
      if(Cp_N < Cp_min ){
        N_best <-Nblocks
        Cp_min <- Cp_N
      }
      
      results[[length(results) + 1]] <- data.frame(
        alpha = alpha, beta = beta, n = n, N = Nblocks,
        h = est$h, theta = est$theta22, sigma = est$sigma2
      )
     
    }
    selected_N <- bind_rows(selected_N, data.frame(alpha = alpha, beta = beta, n = n, N_sel = N_best))
  }
}


df_res <- bind_rows(results)

df_res2 <- left_join(df_res, selected_N, by = c("alpha", "beta", "n"))



# Plot 1: h_mean vs N for Beta(2,2)


chosen <- df_res %>% filter(alpha == 2, beta == 2)
p1 <- ggplot(chosen, aes(x = factor(N), y = h, group = factor(n), color = factor(n))) +
  geom_line(aes(group = factor(n)), lwd = 1) + geom_point(size = 2) +
  labs(x = "Number of blocks N", y = "Estimated h_AMISE",
       title = "h_AMISE vs N (Beta(2,2))", color = "n") +
  theme_minimal()
ggsave("plot_h_vs_Nblock.png", p1, width = 8, height = 6, dpi = 150)
plot(p1)


# Plot 2: h_mean vs n (using selected N by Cp)
p2data <- df_res2 %>% group_by(alpha, beta, n) %>%
  filter(N == N_sel) %>% ungroup()
p2 <- ggplot(p2data, aes(x = n, y = h, color = factor(paste0("Beta(", alpha, ",", beta, ")")))) +
  geom_line(aes(group = interaction(alpha, beta)), size = 1) + geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "Sample size n (log scale)", y = "Estimated h_AMISE (mean)",
       title = "h_AMISE vs n (N selected by Cp-like rule)", color = "Beta dist") +
  theme_minimal()
ggsave("plot_h_vs_n.png", p2, width = 8, height = 6, dpi = 150)
plot(p2)


# Plot 3: effect of Beta shape on h (n fixed = 1000)
n_fix <- 1000
p3data <- df_res2 %>% filter(n == n_fix) %>%
  group_by(alpha, beta) %>% filter(N == N_sel) %>% ungroup()
p3 <- ggplot(p3data,aes(x = factor(sprintf("(%.1f, %.1f), N=%d", alpha, beta, N)), y = h)) +
  geom_col(fill = "steelblue") + labs(x = "Beta shape (and selected N)", y = "Estimated h_AMISE",
                                      title = sprintf("Effect of Beta shape on h_AMISE (n=%d)", n_fix)) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("plot_h_vs_beta.png", p3, width = 9, height = 5, dpi = 150)
plot(p3)

# Print brief summary
print("Simulation finished. Files produced:")
print(list.files(pattern = "plot_h_.*png"))
print("Saved tables: simulation_results_summary.csv, selected_N_by_Cp.csv")
