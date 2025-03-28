# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Define helper functions

# Custom sign function: returns 1 for >= 0, -1 for < 0
custom_sign <- function(x) {
  sapply(x, function(xi) ifelse(xi >= 0, 1, -1))
}

# Synchronous update rule for Hopfield network: new state = sign(W %*% s)
hopfield_update <- function(s, W) {
  custom_sign(W %*% s)
}

# Convergence function: iterate until the state is stable (or max_iter reached)
converge_state <- function(s, W, max_iter = 100) {
  for (i in 1:max_iter) {
    new_s <- hopfield_update(s, W)
    if (all(new_s == s)) break
    s <- new_s
  }
  return(s)
}

# Compute basin size by testing convergence from sampled states
# Enumerate all states
compute_basin_size <- function(W, target_pattern, N) {
  # First check if the target is a stable fixed point
  stable_target <- converge_state(target_pattern, W)
  if (!all(stable_target == target_pattern)) {
    return(0)
  }
  # Generate all possible states
  states <- as.matrix(expand.grid(replicate(N, c(-1, 1), simplify = FALSE)))
  basin_count <- 0
  for (i in 1:nrow(states)) {
    s_init <- as.numeric(states[i, ])
    final_state <- converge_state(s_init, W)
    if (all(final_state == target_pattern)) {
      basin_count <- basin_count + 1
    }
  }
  return(basin_count)
}

# Function to run one batch of simulations using a multi-memory initialization
run_simulation_batch_multi <- function(num_runs, T2, N, p, p_pres) {
  # Pre-allocate matrix: rows = runs, columns = T2 target presentation steps
  batch_results <- matrix(NA, nrow = num_runs, ncol = T2)
  # Determine progress update intervals for inner and outer loops
  inner_update_interval <- ceiling(T2 / 10)
  outer_update_interval <- ceiling(num_runs / 10)
  for (run in 1:num_runs) {
    # Pre-store p random memories
    memory_pres <- p_pres*matrix(sample(c(-1, 1), p * N, replace = TRUE), 
                                 nrow = p, ncol = N, byrow = TRUE)
    
    # Choose a target stimulus
    stimulus_target <- sample(c(-1, 1), N, replace = TRUE)

    # Pre-allocate vector to store basin sizes for this run
    basin_sizes_run <- numeric(T2)
    
    # For each target presentation count j (Phase 2), build the complete memory matrix
    for (j in 1:T2) {
      # Create memory_phase2: j copies of the target stimulus
      memory_phase2 <- matrix(rep(stimulus_target, j), nrow = j, byrow = TRUE)
      
      # Combine the pre-stored memories with the target presentations.
      x <- rbind(memory_pres, memory_phase2)
      
      # Compute weight matrix using the Hebbian rule (W = t(x) %*% x) and set the diagonal to zero.
      W <- t(x) %*% x
      diag(W) <- 0
      
      # Compute the basin size for the target stimulus.
      basin_sizes_run[j] <- compute_basin_size(W, stimulus_target, N)
      
      # Print inner progress update.
      if (j %% inner_update_interval == 0) {
        cat("Run", run, ":", round((j / T2) * 100), "% target presentations processed.\n")
      }
    }
    batch_results[run, ] <- basin_sizes_run
    
    # Print outer progress update.
    if (run %% outer_update_interval == 0) {
      cat("Batch progress:", round((run / num_runs) * 100), "% runs completed (run", run, "of", num_runs, ").\n")
    }
  }
  
  # Convert matrix to data frame (wide format) and add RunID
  df <- as.data.frame(batch_results)
  df$RunID <- 1:num_runs
  return(df)
}

# Run the simulation with multi-memory initialization
set.seed(42)  # set seed once

# Set simulation parameters
N <- 10 # number of neurons per pattern
T2 <- 1000 # number of target stimulus presentations (Phase 2 steps)
p <- 50 # number of pre-stored random memories
p_pres=10 # weighting factor for pre-stored memories
num_runs <- 30  # number of independent simulation runs

results_multi <- run_simulation_batch_multi(num_runs = num_runs, T2 = T2, N = N, p = p)

# Convert results to long format for analysis
df_long <- pivot_longer(results_multi, 
                        cols = -RunID, 
                        names_to = "PresentationStep", 
                        values_to = "BasinSize")

# The column names for presentation steps may be V1, V2, ... so convert them to numeric
df_long$PresentationStep <- as.numeric(gsub("V", "", df_long$PresentationStep))

# Compute summary statistics: mean and standard error at each presentation step
summary_df <- df_long %>%
  group_by(PresentationStep) %>%
  summarise(MeanBasin = mean(BasinSize, na.rm = TRUE),
            SE = sd(BasinSize, na.rm = TRUE) / sqrt(n()))

# Plot the average growth of the target attractor basin
p1 <- ggplot(summary_df, aes(x = PresentationStep, y = MeanBasin)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Average Growth of Target Attractor Basin",
       x = "Number of Target Presentations",
       y = "Mean Basin Size (Number of States)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

print(p1) 

# Extract the basin size columns (assuming they are named V1, V2, ..., V_T2)
T2 <- ncol(results_multi) - 1  # subtracting the RunID column
basin_matrix <- as.matrix(results_multi[, 1:T2])

# For each run (row), compute differences (jump sizes) between consecutive presentations
# This will create a matrix with one fewer column than the basin_matrix
jump_matrix <- t(apply(basin_matrix, 1, diff))

# Combine all jump sizes into a single vector
all_jump_sizes <- as.vector(jump_matrix)

# Extract positive jumps
positive_jumps <- all_jump_sizes[all_jump_sizes > 0]

