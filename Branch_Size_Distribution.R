############################################
# 1. Libraries
############################################

library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

############################################
# 2. Hopfield Helper Functions
############################################

# Returns 1 if x >= 0 and -1 if x < 0
custom_sign <- function(x) {
  ifelse(x >= 0, 1, -1)
}

# One synchronous update
hopfield_update <- function(s, W) {
  custom_sign(W %*% s)
}

# Iteratively update until convergence (or max_iter reached)
converge_state <- function(s, W, max_iter = 100) {
  for (i in seq_len(max_iter)) {
    new_s <- hopfield_update(s, W)
    if (all(new_s == s)) break
    s <- new_s
  }
  return(s)
}

# Compute the basin size for a given target pattern by enumerating all states
compute_basin_size <- function(W, target_pattern, N) {
  stable_target <- converge_state(target_pattern, W)
  if (!all(stable_target == target_pattern)) return(0)
  
  states <- as.matrix(expand.grid(rep(list(c(-1,1)), N)))
  basin_count <- 0
  for (i in seq_len(nrow(states))) {
    s_init <- as.numeric(states[i, ])
    final_state <- converge_state(s_init, W)
    if (all(final_state == target_pattern)) {
      basin_count <- basin_count + 1
    }
  }
  return(basin_count)
}

# Enumerate all possible states (each row is one state)
enumerate_states <- function(N) {
  as.matrix(expand.grid(rep(list(c(-1,1)), N)))
}

# Label each state with the attractor it converges to under W
# We convert the numeric attractor vector to a string using a comma delimiter
label_all_attractors <- function(W, N) {
  states <- enumerate_states(N)
  state_labels <- apply(states, 1, paste, collapse = ",")  # use comma as delimiter
  attractor_label <- character(length(state_labels))
  for (i in seq_along(state_labels)) {
    s_init <- as.numeric(states[i, ])
    final_s <- converge_state(s_init, W)
    attractor_label[i] <- paste(final_s, collapse = ",")
  }
  names(attractor_label) <- state_labels
  return(attractor_label)
}

# Conversion function: from a comma-separated attractor string back to a numeric vector
convert_attractor_string_to_numeric <- function(attr_str) {
  as.numeric(strsplit(attr_str, ",")[[1]])
}

############################################
# 3. Training a Hopfield Network
############################################

# Train a Hopfield network using repeated random stimuli
# Each network is trained on num_memories random patterns, each presented rep_per_memory times
train_hopfield_network <- function(N, num_memories, rep_per_memory) {
  # Generate num_memories random patterns of length N
  patterns <- matrix(sample(c(-1,1), num_memories * N, replace = TRUE), nrow = num_memories, ncol = N)
  # Repeat each pattern rep_per_memory times
  repeated_patterns <- do.call(rbind, lapply(1:num_memories, function(i) patterns[i, , drop=FALSE][rep(1, rep_per_memory), ]))
  # Compute weight matrix using Hebbian learning
  W <- t(repeated_patterns) %*% repeated_patterns
  diag(W) <- 0
  return(W)
}

############################################
# 4. Graph Construction and Branch Size Calculation
############################################

# Build the full state transition graph
build_full_graph <- function(W, N) {
  states <- enumerate_states(N)
  state_labels <- apply(states, 1, paste, collapse = "")
  edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(states))) {
    s <- as.numeric(states[i, ])
    next_s <- hopfield_update(s, W)
    edges <- rbind(edges, data.frame(
      from = state_labels[i],
      to = paste(next_s, collapse = ""),
      stringsAsFactors = FALSE
    ))
  }
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = state_labels)
  return(g)
}

# Compute branch sizes (descendant counts) for every state
compute_branch_sizes <- function(W, N) {
  # Build the full transition graph
  g <- build_full_graph(W, N)
  # Reverse the graph so that for any vertex v, the reachable set (mode="out") is the set
  # of states that drain into v
  g_rev <- reverse_edges(g, eids = E(g))
  
  descendant_counts <- sapply(V(g_rev)$name, function(vname) {
    reachable <- subcomponent(g_rev, vname, mode = "out")
    count <- length(reachable) - 1  # subtract one to not count the vertex itself
    if (count == 0) 1 else count
  })
  
  return(descendant_counts)
}


############################################
# 5. Overall Analysis Script
############################################

# Parameters for the network and training:
N <- 10                 # number of neurons per pattern 
num_memories <- 50      # number of random training patterns
rep_per_memory <- 10    # number of presentations per pattern

# How many random networks to generate:
num_networks <- 1000

all_branch_sizes <- c()

for (net in 1:num_networks) {
  cat("Training network", net, "of", num_networks, "\n")
  # Train a network
  W <- train_hopfield_network(N, num_memories, rep_per_memory)
  # Compute branch sizes (descendant counts) across all states in the reversed graph
  branch_sizes <- compute_branch_sizes(W, N)
  all_branch_sizes <- c(all_branch_sizes, branch_sizes)
}

# Convert branch sizes to a data frame for analysis
df_branch <- data.frame(branch_size = all_branch_sizes)

# Use Distribution_Analysis script to analyze distribution
