setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(doParallel)

# ==============================================================================
# Filenames
# ==============================================================================

seeds_file <- "seeds.csv"
models_file <- "models.rds"
results_file <- "results.csv"

batch_prefix <- "batch-"

# ==============================================================================
# Seeds and Models
# ==============================================================================

# Read seeds from file
read_csv(seeds_file, col_types = "ic") |> 
    arrange(seed) |> 
    write_csv(seeds_file)
seeds <- read_csv(seeds_file, col_types = "ic") |> 
    pull(seed)

# Compute models for any seeds not found in saved models file
# Tries to read from "models.csv" in the current directory, if possible, 
# and writes to that file if there are models haven't yet been generated.
models <- map(seeds, 
              \(s) generate_simulation_models(s, models_file = models_file))

# ==============================================================================
# Sanity Checks
# ==============================================================================

# Is there variety of model dimensionalities? 
models |> 
    map_chr(\(m) str_flatten(map_dbl(m$dimnames, length), collapse = "x")) |> 
    str_flatten(collapse = " ")
# Is there variety (but no extremes) of probability of group 0? 
models |> 
    map_dbl(\(m) m$group_dist[1]) 
# Is there variety of distributions of ability conditioned on group?
models |> 
    map(\(m) m$ability_dist)
# How long did Markov move computations take?
models |> 
    map(\(m) m$moves_computation_time)
# What are the KL divergences of all of the joint distributions like?
models |> 
    map(\(m) m$config) |> 
    bind_rows()

# ==============================================================================
# Sample Sizes and Number of Simulations
# ==============================================================================

# These are the desired number of sample sizes and the desired number of
# simulations per configuration

sample_sizes <- seq(25, 2500, by = 25)
goal_num_sims <- 500

# ==============================================================================
# Simulation Batch Setup
# ==============================================================================

if (file.exists(results_file)) {
    # Read in existing data to update required numbers of simulation
    counts <- read_csv(results_file, show_col_types = FALSE) |> 
        group_by(seed, type, norm, sample_size) |> 
        summarize(num_sims = goal_num_sims - n(), .groups = "drop")
    
    # Visual check for batches with more than desired number of sims
    counts |> 
        filter(num_sims < 0) |> 
        print()
} else {
    # Initialize empty tibble
    counts <- tibble(seed = numeric(), 
                     type = character(), 
                     norm = numeric(), 
                     sample_size = numeric(), 
                     num_sims = numeric())
}

# Set up batches of simulations
batches <- models |>
    map(\(m) m$config) |> 
    bind_rows() |> 
    expand_grid(sample_size = sample_sizes) |> 
    mutate(num_sims = goal_num_sims) |> 
    rows_update(counts, by = c("seed", "type", "norm", "sample_size")) |> 
    filter(num_sims > 0) |> 
    mutate(id = seq_along(seed)) |>
    relocate(id)

# Run one batch of simulations. The configurations of the batch are specified 
# by c, which is one row of the batches tibble as it is set up above. 
# If prefix is a not NULL, the results are written to a temporary csv file 
# whose name is prefix followed by the batch id followed by ".csv". 
# It also returns the results as a tibble. 
run_batch <- function(c, prefix = NULL) {
    moves <- models[[match(c$seed, seeds)]]$moves
    
    lambda <- function(x) {
        t <- simulate_table(c$dist[[1]], c$sample_size)
        bind_cols(seed = c$seed, 
                  type = c$type,
                  norm = c$norm,
                  sample_size = c$sample_size,
                  table_task(t, moves))
    }
    
    results <- map(1:(c$num_sims), lambda) |> 
        bind_rows()
    if (!is_null(prefix)) {
        filename <- paste(prefix, c$id, ".csv", sep="")
        write_csv(results, filename)
    }
    results
}

# ==============================================================================
# Parallelized Simulations
# ==============================================================================

num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)
registerDoParallel(cl = cl)
clusterEvalQ(cl, library(tidyverse))

start_timer <- Sys.time()
foreach(i = 1:nrow(batches), 
        .combine = \(x) NULL, 
        .inorder = FALSE) %dopar% {
    run_batch(slice(batches, i), batch_prefix)
}
stopCluster(cl = cl)
simulation_computation_time <- Sys.time() - start_timer

# ==============================================================================
# Clean Up
# ==============================================================================

# Read all batch files
batch_files <- str_c(batch_prefix, batches$id, ".csv") |> 
    keep(file.exists)
results <- batch_files |> 
    map(\(f) read_csv(f, show_col_types = FALSE)) |> 
    bind_rows()

# Consolidate with existing results
if (file.exists(results_file)) {
    results <- read_csv(results_file, show_col_types = FALSE) |> 
        bind_rows(results)
}

# Write into results
write_csv(results, results_file)

# Remove batch files
file.remove(batch_files)

# Display computation time
simulation_computation_time