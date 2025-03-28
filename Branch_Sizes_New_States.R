############################################
# 1. Libraries
############################################

library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

############################################
# 2. Parameters
############################################
N <- 10                # number of neurons per pattern
T2 <- 1000             # number of target presentations (Phase 2)
p <- 50                # number of pre-stored memories
p_pres <- 10           # pretraining weight factor (each memory multiplied by p_pres)
num_runs <- 30          # number of independent simulation runs
jump_threshold <- 0    # jump detection threshold

############################################
# 3. Hopfield Helper Functions
############################################

# Custom sign: returns 1 if x>=0, -1 otherwise
custom_sign <- function(x) {
  ifelse(x >= 0, 1, -1)
}

# One synchronous update
hopfield_update <- function(s, W) {
  custom_sign(W %*% s)
}

# Convergence: iterate until state is stable (or max_iter reached)
converge_state <- function(s, W, max_iter = 100) {
  for (i in 1:max_iter) {
    new_s <- hopfield_update(s, W)
    if (all(new_s == s)) break
    s <- new_s
  }
  return(s)
}

# Compute basin size by enumerating all states
compute_basin_size <- function(W, target_pattern, N) {
  stable_target <- converge_state(target_pattern, W)
  if (!all(stable_target == target_pattern)) return(0)
  states <- as.matrix(expand.grid(rep(list(c(-1,1)), N)))
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

# Enumerate all states
enumerate_states <- function(N) {
  as.matrix(expand.grid(rep(list(c(-1,1)), N)))
}

# Label each state with the attractor it converges to under W
label_all_attractors <- function(W, N) {
  states <- enumerate_states(N)
  state_labels <- apply(states, 1, paste, collapse = "")
  attractor_label <- character(length(state_labels))
  for (i in seq_along(state_labels)) {
    s_init <- as.numeric(states[i, ])
    final_s <- converge_state(s_init, W)
    attractor_label[i] <- paste(final_s, collapse = "")
  }
  names(attractor_label) <- state_labels
  return(attractor_label)
}

# Parser: Convert a state label (e.g., "1-1-1-1") back to a numeric vector
# (Assumes labels were produced with no delimiter and using the rule that every number ends with "1")
parse_attractor_string <- function(s) {
  n <- nchar(s)
  i <- 1
  result <- numeric(0)
  while (i <= n) {
    ch <- substr(s, i, i)
    if (ch == "1") {
      result <- c(result, 1)
      i <- i + 1
    } else if (ch == "-") {
      if (i + 1 <= n) {
        next_ch <- substr(s, i + 1, i + 1)
        if (next_ch == "1") {
          result <- c(result, -1)
          i <- i + 2
        } else {
          stop("Unexpected character following '-' in attractor string.")
        }
      } else {
        stop("Attractor string ends with '-' without a following digit.")
      }
    } else {
      stop(paste("Unexpected character in attractor string:", ch))
    }
  }
  return(result)
}

# Get immediate parent of a state using W
get_parent_label <- function(state_label, W) {
  s <- parse_attractor_string(state_label)
  parent <- hopfield_update(s, W)
  parent_label <- paste(parent, collapse = "")
  return(parent_label)
}

############################################
# 4. Graph Construction and Efficient Branch Size Computation for New States
############################################

# We build a graph only on the new state labels
# new_states: vector of state labels that are new (i.e. in post basin but not in pre basin)
# Returns a vector of branch sizes (one value per branch/component)
compute_branch_sizes_new <- function(W, new_states) {
  # Build edges only among new states
  edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
  for (v in new_states) {
    parent <- get_parent_label(v, W)
    # If the parent's label is also in new_states, add an edge (v -> parent)
    if (parent %in% new_states) {
      edges <- rbind(edges, data.frame(from = v, to = parent, stringsAsFactors = FALSE))
    }
  }
  # Create a graph using only the new states
  if(nrow(edges) > 0) {
    g_new <- graph_from_data_frame(edges, directed = TRUE, vertices = new_states)
  } else {
    # If there are no edges, then each state is isolated
    g_new <- make_empty_graph(n = length(new_states), directed = TRUE)
    V(g_new)$name <- new_states
  }
  # Reverse the graph so that "descendants" represent the new states draining to a branch head
  g_new_rev <- reverse_edges(g_new, eids = E(g_new))
  # Compute connected components
  comps <- components(g_new_rev)
  # For each connected component, the branch size is the number of nodes in that component
  branch_sizes <- comps$csize
  # Return as a vector
  return(branch_sizes)
}

############################################
# 5. Simulation: Run Batch of Simulations with Multi-Memory Initialization
############################################

run_simulation_batch_multi <- function(num_runs, T2, N, p, p_pres) {
  all_basin_sizes <- matrix(NA, nrow = num_runs, ncol = T2)
  all_weight_matrices <- vector("list", num_runs)
  target_patterns <- vector("list", num_runs)
  
  inner_update_interval <- ceiling(T2 / 10)
  outer_update_interval <- ceiling(num_runs / 10)
  
  for (run in 1:num_runs) {
    memory_pres <- p_pres * matrix(sample(c(-1, 1), p * N, replace = TRUE),
                                   nrow = p, ncol = N, byrow = TRUE)
    repeat {
      stimulus_target <- sample(c(-1, 1), N, replace = TRUE)
      if (!any(apply(memory_pres, 1, function(row) all(row == stimulus_target)))) break
    }
    target_patterns[[run]] <- stimulus_target
    
    basin_sizes_run <- numeric(T2)
    weight_list_run <- vector("list", T2)
    
    for (j in 1:T2) {
      memory_phase2 <- matrix(rep(stimulus_target, j), nrow = j, byrow = TRUE)
      x <- rbind(memory_pres, memory_phase2)
      W <- t(x) %*% x
      diag(W) <- 0
      weight_list_run[[j]] <- W
      basin_sizes_run[j] <- compute_basin_size(W, stimulus_target, N)
      
      if (j %% inner_update_interval == 0) {
        cat("Run", run, ":", round((j / T2) * 100), "% target presentations processed.\n")
      }
    }
    all_basin_sizes[run, ] <- basin_sizes_run
    all_weight_matrices[[run]] <- weight_list_run
    
    if (run %% outer_update_interval == 0) {
      cat("Batch progress:", round((run / num_runs) * 100),
          "% runs completed (run", run, "of", num_runs, ").\n")
    }
  }
  
  df_basin <- as.data.frame(all_basin_sizes)
  df_basin$RunID <- 1:num_runs
  return(list(basin_df = df_basin, weight_matrices = all_weight_matrices, target_patterns = target_patterns))
}

############################################
# 6. Jump and Branch Analysis
############################################

# For a single simulation run, detect jump events and compute branch sizes.
analyze_jumps <- function(W_list, target_pattern, N, jump_threshold) {
  basin_sizes <- sapply(W_list, function(W) compute_basin_size(W, target_pattern, N))
  jump_events <- which(diff(basin_sizes) > jump_threshold)
  cat("Detected jump events at steps:", jump_events, "\n")
  
  branch_sizes <- c()
  for (j in jump_events) {
    W_pre <- W_list[[j]]
    W_post <- W_list[[j + 1]]
    
    basin_target_pre <- get_basin_membership(W_pre, target_pattern, N)
    basin_target_post <- get_basin_membership(W_post, target_pattern, N)
    
    new_states <- setdiff(basin_target_post, basin_target_pre)
    cat("Jump at step", j, ": Newly captured states =", length(new_states), "\n")
    
    if (length(new_states) > 0) {
      # Compute branch sizes only for new states.
      branch_sizes_new <- compute_branch_sizes_new(W_post, new_states)
      cat("Branch sizes at jump", j, ":", branch_sizes_new, "\n")
      branch_sizes <- c(branch_sizes, branch_sizes_new)
    }
  }
  return(branch_sizes)
}

############################################
# 7. Main Script: Run Simulations and Analyze Branch Sizes at Jumps
############################################

set.seed(42)
batch2 <- run_simulation_batch_multi(num_runs, T2, N, p, p_pres)

batch_next <- run_simulation_batch_multi(num_runs, T2, N, p, p_pres)

batch1=sim_result

basin_df <- sim_result$basin_df
W_matrices_all <- sim_result$weight_matrices
target_patterns <- sim_result$target_patterns

basin_df <- c(batch1$basin_df, batch2$basin_df)
W_matrices_all <- c(batch1$weight_matrices, batch2$weight_matrices)
target_patterns <- c(batch1$target_patterns, batch2$target_patterns)

basin_df <- c(basin_df, batch_next$basin_df)
W_matrices_all <- c(W_matrices_all, batch_next$weight_matrices)
target_patterns <- c(target_patterns, batch_next$target_patterns)

# For each run, detect jumps and collect branch sizes
all_branch_sizes <- c()
for (run in 1:num_runs) {
  cat("Analyzing run", run, "\n")
  W_list <- W_matrices_all[[run]]
  target_pattern_run <- target_patterns[[run]]
  branch_sizes_run <- analyze_jumps(W_list, target_pattern_run, N, jump_threshold)
  all_branch_sizes <- c(all_branch_sizes, branch_sizes_run)
}

cat("Total branches measured:", length(all_branch_sizes), "\n")

# Use Distribution_Analysis script to analyze distribution
