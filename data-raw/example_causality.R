# Script to generate example causality data

set.seed(42)
n <- 200

# Generate X as a random walk
cause_x <- cumsum(rnorm(n, sd = 1))

# Generate Y that depends on lagged X (X Granger-causes Y)
# Y_t = 0.7 * X_{t-1} + epsilon_t
effect_y <- c(0, 0.7 * cause_x[1:(n-1)]) + rnorm(n, sd = 0.5)

# Create data frame
example_causality <- data.frame(
  time = 1:n,
  cause_x = cause_x,
  effect_y = effect_y
)

# Save to data directory
usethis::use_data(example_causality, overwrite = TRUE)
