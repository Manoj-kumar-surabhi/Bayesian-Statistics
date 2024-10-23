##problem, 1.1
mu_0 <- 5
sigma_0_squared <- 4
kappa_0 <- 1
v_0 <- 2

school1 <- scan("school1.dat")
school2 <- scan("school2.dat")
school3 <- scan("school3.dat")

posterior_sampling <- function(data, mu_0, sigma_0_squared, kappa_0, v_0, num_samples = 10000) {
  n <- length(data)
  sample_mean <- mean(data)
  sample_var <- var(data)
  
  # Posterior parameters
  kappa_n <- kappa_0 + n
  mu_n <- (kappa_0 * mu_0 + n * sample_mean) / kappa_n
  v_n <- v_0 + n
  sigma_n_squared <- (v_0 * sigma_0_squared + (n - 1) * sample_var + 
                        (kappa_0 * n * (sample_mean - mu_0)^2) / kappa_n) / v_n
  
# Monte Carlo sampling
  sigma_squared_samples <- 1 / rgamma(num_samples, v_n / 2, v_n * sigma_n_squared / 2)
  theta_samples <- rnorm(num_samples, mu_n, sqrt(sigma_squared_samples / kappa_n))
  
  return(list(theta_samples = theta_samples, sigma_squared_samples = sqrt(sigma_squared_samples)))
}

set.seed(123)  
posterior_school1 <- posterior_sampling(school1, mu_0, sigma_0_squared, kappa_0, v_0)
posterior_school2 <- posterior_sampling(school2, mu_0, sigma_0_squared, kappa_0, v_0)
posterior_school3 <- posterior_sampling(school3, mu_0, sigma_0_squared, kappa_0, v_0)

compute_summary <- function(samples) {
  mean_val <- mean(samples)
  ci <- quantile(samples, probs = c(0.025, 0.975))
  return(list(mean = mean_val, ci = ci))
}

# Summarize results
school1_theta_summary <- compute_summary(posterior_school1$theta_samples)
school1_sigma_summary <- compute_summary(posterior_school1$sigma_squared_samples)

school2_theta_summary <- compute_summary(posterior_school2$theta_samples)
school2_sigma_summary <- compute_summary(posterior_school2$sigma_squared_samples)

school3_theta_summary <- compute_summary(posterior_school3$theta_samples)
school3_sigma_summary <- compute_summary(posterior_school3$sigma_squared_samples)

# Print results
print("School 1: Theta Posterior Mean and 95% CI")
print(school1_theta_summary)

print("School 1: Sigma Posterior Mean and 95% CI")
print(school1_sigma_summary)

print("School 2: Theta Posterior Mean and 95% CI")
print(school2_theta_summary)

print("School 2: Sigma Posterior Mean and 95% CI")
print(school2_sigma_summary)

print("School 3: Theta Posterior Mean and 95% CI")
print(school3_theta_summary)

print("School 3: Sigma Posterior Mean and 95% CI")
print(school3_sigma_summary)

###1.2
theta_samples <- cbind(posterior_school1$theta_samples, 
                       posterior_school2$theta_samples, 
                       posterior_school3$theta_samples)

compute_prob <- function(theta_samples, i, j, k) {
  return(mean(theta_samples[, i] < theta_samples[, j] & theta_samples[, j] < theta_samples[, k]))
}

prob_123 <- compute_prob(theta_samples, 1, 2, 3)
prob_132 <- compute_prob(theta_samples, 1, 3, 2)
prob_213 <- compute_prob(theta_samples, 2, 1, 3)
prob_231 <- compute_prob(theta_samples, 2, 3, 1)
prob_312 <- compute_prob(theta_samples, 3, 1, 2)
prob_321 <- compute_prob(theta_samples, 3, 2, 1)

cat("Posterior probabilities for theta_i < theta_j < theta_k:\n")
cat("P(theta_1 < theta_2 < theta_3) =", prob_123, "\n")
cat("P(theta_1 < theta_3 < theta_2) =", prob_132, "\n")
cat("P(theta_2 < theta_1 < theta_3) =", prob_213, "\n")
cat("P(theta_2 < theta_3 < theta_1) =", prob_231, "\n")
cat("P(theta_3 < theta_1 < theta_2) =", prob_312, "\n")
cat("P(theta_3 < theta_2 < theta_1) =", prob_321, "\n")

##1.3

generate_posterior_predictive <- function(theta_samples, sigma_samples, num_samples) {
  return(rnorm(num_samples, mean = theta_samples, sd = sigma_samples))
}

# Number of samples
num_samples <- length(posterior_school1$theta_samples)

# Generate posterior predictive samples for each school
y_tilde_school1 <- generate_posterior_predictive(posterior_school1$theta_samples, posterior_school1$sigma_squared_samples, num_samples)
y_tilde_school2 <- generate_posterior_predictive(posterior_school2$theta_samples, posterior_school2$sigma_squared_samples, num_samples)
y_tilde_school3 <- generate_posterior_predictive(posterior_school3$theta_samples, posterior_school3$sigma_squared_samples, num_samples)

# Combine the posterior predictive samples into a matrix
y_tilde_samples <- cbind(y_tilde_school1, y_tilde_school2, y_tilde_school3)

# Function to compute probability for each permutation
compute_prob_y_tilde <- function(y_tilde_samples, i, j, k) {
  return(mean(y_tilde_samples[, i] < y_tilde_samples[, j] & y_tilde_samples[, j] < y_tilde_samples[, k]))
}

# Compute the posterior probability for all six permutations
prob_y_123 <- compute_prob_y_tilde(y_tilde_samples, 1, 2, 3)
prob_y_132 <- compute_prob_y_tilde(y_tilde_samples, 1, 3, 2)
prob_y_213 <- compute_prob_y_tilde(y_tilde_samples, 2, 1, 3)
prob_y_231 <- compute_prob_y_tilde(y_tilde_samples, 2, 3, 1)
prob_y_312 <- compute_prob_y_tilde(y_tilde_samples, 3, 1, 2)
prob_y_321 <- compute_prob_y_tilde(y_tilde_samples, 3, 2, 1)

# Print the results
cat("Posterior probabilities for Y_i < Y_j < Y_k:\n")
cat("P(Y_1 < Y_2 < Y_3) =", prob_y_123, "\n")
cat("P(Y_1 < Y_3 < Y_2) =", prob_y_132, "\n")
cat("P(Y_2 < Y_1 < Y_3) =", prob_y_213, "\n")
cat("P(Y_2 < Y_3 < Y_1) =", prob_y_231, "\n")
cat("P(Y_3 < Y_1 < Y_2) =", prob_y_312, "\n")
cat("P(Y_3 < Y_2 < Y_1) =", prob_y_321, "\n")


#1.4:
# Posterior probability that theta_1 > theta_2 and theta_3
prob_theta_1_greater <- mean(posterior_school1$theta_samples > posterior_school2$theta_samples & 
                               posterior_school1$theta_samples > posterior_school3$theta_samples)

# Posterior probability that Y_tilde_1 > Y_tilde_2 and Y_tilde_3
prob_y_tilde_1_greater <- mean(y_tilde_school1 > y_tilde_school2 & y_tilde_school1 > y_tilde_school3)

cat("P(theta_1 > theta_2 and theta_3) =", prob_theta_1_greater, "\n")
cat("P(Y_tilde_1 > Y_tilde_2 and Y_tilde_3) =", prob_y_tilde_1_greater, "\n")



###2.1
mu <- 0
sigma <- 0.25
n <- 10
y <- 6

theta_grid <- seq(-2, 2, length.out = 1000)

likelihood <- function(theta, y, n) {
  p <- exp(theta) / (1 + exp(theta)) # logistic function
  return((p^y) * ((1 - p)^(n - y)))
}

prior <- function(theta, mu, sigma) {
  return(dnorm(theta, mean = mu, sd = sigma))
}

unnormalized_posterior <- likelihood(theta_grid, y, n) * prior(theta_grid, mu, sigma)

posterior <- unnormalized_posterior / sum(unnormalized_posterior)

pr_theta_greater_than_0 <- sum(posterior[theta_grid > 0])

cdf <- cumsum(posterior)
lower_bound <- theta_grid[which(cdf >= 0.025)[1]]
upper_bound <- theta_grid[which(cdf >= 0.975)[1]]

cat("Pr(theta > 0 | y) =", pr_theta_greater_than_0, "\n")
cat("95% credible interval for theta: [", lower_bound, ",", upper_bound, "]\n")

