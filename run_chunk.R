# run_chunk.R

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

source("sim_setup.R")

# Fixed simulation parameters
n <- 200
p <- 15
eps_vec <- c(0.17, 0.2, 0.25, 0.33, 0.5, 0.67, 0.75, 0.8, 0.83)
sigma <- sqrt(2)

d_list <- list(
  c(100,80,60,40,20),
  c(150,120,120),
  c(1, 0.9, 0.8),
  c(100,90,80),
  c(80,80,60)
)

# Map Slurm task ID -> trial + counter
reps_per_d <- 200

trial <- ((task_id - 1) %/% reps_per_d) + 1
counter <- ((task_id - 1) %% reps_per_d) + 1

d <- d_list[[trial]]
seed <- task_id

res <- run_one_replicate(
  seed = seed,
  counter = task_id,
  d = d,
  n = n,
  p = p,
  sigma = sigma,
  eps_vec = eps_vec
)

dir.create("results", showWarnings = FALSE)
outfile <- sprintf("results/sim_%04d.csv", task_id)
readr::write_csv(res, outfile)
