##################################################################
## Gasser-Muller estimator with Dirichlet kernel on the simplex ##
##################################################################

## Written by Frederic Ouimet (November 2024)

require("cubature")             # for integrals
require("doFuture")             # for parallel execution with foreach
require("fs")                   # for filesystem operations (dependency)
require("future.batchtools")    # for batchtools integration with future
require("ggplot2")              # for plotting
require("LaplacesDemon")        # for the Dirichlet distribution
require("parallel")             # for parallel execution of calculations
require("tidyverse")            # for data manipulation and visualization
require("writexl")              # to export output to Excel files

##############################
## Parallelization on cores ##
##############################

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "cubature",
  "doFuture",
  "fs",
  "future.batchtools",
  "ggplot2",
  "LaplacesDemon",
  "parallel",
  "tidyverse",
  "writexl"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "BB",
  "JJ",
  "KK",
  "LSCV",
  "LSCV_MC",
  "MCsim",
  "MM",
  "RR",
  "adaptIntegrate",
  "b_grid",
  "b_opt",
  "b_opt_MC",
  "b_opt_MC_grid",
  "b_opt_LOOCV",
  "b_opt_grid",
  "b_opt_value_MC_grid",
  "d",
  "elapsed_time",
  "hat_m",
  "integrand",
  "j",
  "LOOCV",
  "m",
  "mesh",
  "method",
  "mm",
  "n_to_k",
  "num_cores",
  "path",
  "results",
  "results_grid",
  "RSS_LOOCV",
  "s_grid",
  "start_time",
  "tol1",
  "tol2",
  "unif_sample",
  "xx",
  "y"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- detectCores() - 1
  cl <- makeCluster(num_cores) # Create the cluster
  
  # Export the list of libraries to the worker nodes
  clusterExport(cl, varlist = "libraries_to_load")
  
  # Load necessary libraries on each cluster node
  invisible(clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))
  
  # Export all necessary objects, functions, and parameters to the worker nodes
  clusterExport(cl, varlist = vars_to_export)
  
  return(cl) # Return the cluster object
}

# Initialize all variables in the list as NULL except vars_to_export and setup_parallel_cluster

invisible(
  lapply(
    vars_to_export[!(vars_to_export %in% c("vars_to_export", "setup_parallel_cluster"))],
    function(x) assign(x, NULL, envir = .GlobalEnv)
  )
)

################
# Set the path #
################

# path <- file.path(
#   "C://Users//fred1//Desktop//Github_GasserMuller",
#   fsep = .Platform$file.sep
# )
path <- getwd()
# setwd(path)

##############
# Parameters #
##############

d <- 2 # dimension of simplex
MCsim <- 10 ^ 3 # number of uniforms sampled for integral MC estimates

cores_per_node <- 39 # number of cores for each node in the super-computer
BB <- seq(0.01, 1, length.out = cores_per_node) # bandwidths for LSCV_MC graphs

MM <- list("GM", "LL", "NW") # list of Dirichlet kernel methods
KK <- c(7, 10, 14) # indices for the mesh
JJ <- 1:3 # target regression function indices
RR <- 1:100 # replication indices

tol1 <- 1e-3
tol2 <- 1e-1

##############################
## Parallelization on nodes ##
##############################

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "180G",
  walltime = "24:00:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)

###########################
## Mesh of design points ##
###########################

# Mesh of the 2-dim simplex

mesh <- function(k) {
  w <- (k - 1 / sqrt(2)) / (k - 1)
  res <- list()
  for (i in 1:k) {
    for (j in i:k) {
      res <- append(res, list((c(w * (i - 1) + 0.5, w * (k - j) + 0.5) / (k + 1))))
    }
  }
  return(res)
}

# Number of points on each side of the mesh

n_to_k <- function(n) {
  (-1 + sqrt(1 + 8 * n)) / 2
}

#################################
## Target regression functions ##
#################################

# for one design point x

m <- function(j, x) { # x is a d-dim vector on the simplex
  if (j == 1) {
    # Case when j = 1
    res <- x[1] * x[2]
  } else if (j == 2) {
    # Case when j = 2
    res <- (x[1] - 0.25) ^ 2 + (x[2] - 0.75) ^ 2
  } else if (j == 3) {
    # Case when j = 3
    res <- x[1] * exp(-x[2])
  } else {
    # Default case if j is not 1, 2, or 3
    warning("Invalid value of j. Should be 1, 2, or 3.")
    res <- NULL
  }
  return(res)
}

# for a mesh of design points x_1, ..., x_n

mm <- function(j, xx) { # xx is a list of d-dim vectors on the simplex
  res <- list()
  n <- length(xx)
  for (i in 1:n) {
    res <- append(res, list(m(j, xx[[i]])))
  }
  return(res)
}

################
## Estimators ##
################

hat_m <- function(xx, b, s, j, method, y = NULL) {
  # xx is a list of d-dim vectors on the simplex, s is a d-dim vector on the simplex
  n <- length(xx)
  d <- length(xx[[1]])
  
  if (is.null(y)) {
    y <- as.numeric(mm(j, xx)) # without random noise (this is not observed)
    y <- y + 0.1 * IQR(y) * rnorm(n, 0, 1) # with random noise (this is observed)
  }
  
  u <- s / b + rep(1, d)
  v <- (1 - sum(s)) / b + 1
  
  switch(method,
         "GM" = { # GM is implemented specifically for the fixed mesh here, the general definition is different
           # GM method (Gasser-Muller)
           k <- n_to_k(n)
           w <- (k - 1 / sqrt(2)) / (k - 1)
           half_width <- w / (k + 1) / 2
           # k indices for the closest points to the line y(x) = 1 - x
           cpi <- order(sapply(xx, function(x) abs(1 - x[1] - x[2]) / sqrt(2)))[1:k]
           integrand <- function(x) {
             # Check if x is inside the simplex (sum(x) <= 1)
             if (sum(x) >= 1) {
               return(0)  # Return 0 if x is outside the simplex
             } else {
               if (abs(1 - x[1] - x[2]) / sqrt(2) >= (1 / (2 * (k + 1)))) {  # if x is regular
                 for (i in 1:n) {
                   if (abs(x[1] - xx[[i]][1]) <= half_width && abs(x[2] - xx[[i]][2]) <= half_width) {
                     return(y[i] * LaplacesDemon::ddirichlet(c(x, 1 - sum(x)), c(u, v), log = FALSE))
                   }
                 }
               } else { # if x is not regular
                 # Find the index i for which x is closest to xx[[i]] among all j in cpi
                 closest_i <- cpi[which.min(sapply(cpi, function(j) dist(rbind(x, xx[[j]]))))]
                 
                 # Return the result for the closest index
                 return(y[closest_i] * LaplacesDemon::ddirichlet(c(x, 1 - sum(x)), c(u, v), log = FALSE))
               }
             }
           }
           return(adaptIntegrate(integrand, lowerLimit = c(0, 0), upperLimit = c(1, 1), tol = tol1)$integral)
         },
         
         "LL" = {
           # LL method (Local Linear)
           kernel_vec <- rep(NA, n)
           for (i in 1:n) {
             kernel_vec[i] <- LaplacesDemon::ddirichlet(c(xx[[i]], 1 - sum(xx[[i]])), c(u, v), log = FALSE)
           }
           design_mat <- matrix(1, nrow = n, ncol = d + 1)
           for (i in 1:n) {
             design_mat[i, -1] <- xx[[i]] - s
           }
           W <- diag(kernel_vec)
           return(solve(t(design_mat) %*% W %*% design_mat, t(design_mat) %*% W %*% y)[1])
         },
         
         "NW" = {
           # NW method (Nadaraya-Watson)
           kernel_vec <- rep(NA, n)
           for (i in 1:n) {
             kernel_vec[i] <- LaplacesDemon::ddirichlet(c(xx[[i]], 1 - sum(xx[[i]])), c(u, v), log = FALSE)
           }
           return(sum(y * kernel_vec) / sum(kernel_vec))
         },
         
         stop("Invalid method. Choose either 'GM', 'LL', or 'NW'.")
  )
}

# # Set parameters
# xx <- mesh(7) # Generate the mesh of design points
# b <- 0.1
# hat_m_results <- list()
# 
# # Test hat_m for all j values and methods
# for (j in JJ) {
#   cat(paste("\nTesting hat_m for j =", j, "\n"))
#   
#   # Start timer for this j
#   start_time <- Sys.time()
#   
#   hat_m_results[[j]] <- list()
#   for (method in MM) {
#     result <- hat_m(xx, b, s = xx[[1]], j, method) # Using s = xx[[1]] as an example
#     hat_m_results[[j]][[method]] <- result
#   }
#   
#   # End timer for this j
#   end_time <- Sys.time()
#   
#   # Calculate elapsed time in seconds
#   elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
#   
#   # Print results for each j
#   for (method in MM) {
#     cat(paste("Result for method", method, "with j =", j, ":", hat_m_results[[j]][[method]], "\n"))
#   }
#   
#   # Print the time taken for this j
#   cat(paste("Time taken for j =", j, ":", round(elapsed_time, 2), "seconds\n"))
# }

#########################################################
## Least-Squares Cross-Validation (LSCV) exact version ##
#########################################################

LSCV <- function(xx, b, j, method) {
  integrand <- function(s) {
    # Check if s is inside the simplex (sum(s) <= 1)
    if (sum(s) >= 1) {
      return(0)  # Return 0 if s is outside the simplex
    } else {
      return((hat_m(xx, b, s, j, method) - m(j, s))^2)
    }
  }
  
  return(adaptIntegrate(integrand, lowerLimit = c(0, 0), upperLimit = c(1, 1), tol = tol2)$integral)
}

# # Set parameters
# xx <- mesh(7) # Generate the mesh of design points
# b <- 0.1
# lscv_results <- list()
# 
# # Test LSCV for all j values and methods
# for (j in JJ) {
#   cat(paste("\nTesting LSCV for j =", j, "\n"))
#   
#   # Start timer for this j
#   start_time <- Sys.time()
#   
#   lscv_results[[j]] <- list()
#   for (method in MM) {
#     result <- LSCV(xx, b, j, method)
#     lscv_results[[j]][[method]] <- result
#   }
#   
#   # End timer for this j
#   end_time <- Sys.time()
#   
#   # Calculate elapsed time in minutes
#   elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
#   
#   # Print results for each j
#   for (method in MM) {
#     cat(paste("LSCV value for method", method, "with j =", j, ":", lscv_results[[j]][[method]], "\n"))
#   }
#   
#   # Print the time taken for this j
#   cat(paste("Time taken for j =", j, ":", round(elapsed_time, 2), "minutes\n"))
# }

##################################################################
## Least-Squares Cross-Validation (LSCV_MC) Monte Carlo version ##
##################################################################

LSCV_MC <- function(xx, b, j, method, unif_sample = NULL) {
  d <- length(xx[[1]])
  
  # Generate unif_sample if not provided
  if (is.null(unif_sample)) {
    unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  }
  
  integ_squared <- rep(NA, MCsim)
  for (i in 1:MCsim) {
    integ_squared[i] <- (hat_m(xx, b, unif_sample[i, 1:d], j, method) - m(j, unif_sample[i, 1:d])) ^ 2
  }
  
  return(1 / factorial(d) * mean(integ_squared))
}

# # Set parameters
# xx <- mesh(7) # Generate the mesh of design points
# b <- 0.1
# MCsim <- 10^3 # Number of Monte Carlo simulations
# unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, length(xx[[1]]) + 1))
# lscv_mc_results <- list()
# 
# # Test LSCV_MC for all j values and methods
# for (j in JJ) {
#   cat(paste("\nTesting LSCV_MC for j =", j, "\n"))
#   
#   # Start timer for this j
#   start_time <- Sys.time()
#   
#   lscv_mc_results[[j]] <- list()
#   for (method in MM) {
#     result <- LSCV_MC(xx, b, j, method, unif_sample)
#     lscv_mc_results[[j]][[method]] <- result
#   }
#   
#   # End timer for this j
#   end_time <- Sys.time()
#   
#   # Calculate elapsed time in minutes
#   elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
#   
#   # Print results for each j
#   for (method in MM) {
#     cat(paste("LSCV_MC value for method", method, "with j =", j, ":", lscv_mc_results[[j]][[method]], "\n"))
#   }
#   
#   # Print the time taken for this j
#   cat(paste("Time taken for j =", j, ":", round(elapsed_time, 2), "minutes\n"))
# }

#######################################
## Optimal Bandwidth (exact version) ##
#######################################

b_opt <- function(xx, j, method) {
  objective_function <- function(b) {
    return(LSCV(xx, b, j, method))
  }
  res <- optimize(objective_function, interval = c(min(BB), max(BB)))
  return(res$minimum)
}

b_opt_MC <- function(xx, j, method, unif_sample = NULL) {
  # Generate unif_sample if not provided
  if (is.null(unif_sample)) {
    unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  }
  
  objective_function <- function(b) {
    return(LSCV_MC(xx, b, j, method, unif_sample))
  }
  res <- optimize(objective_function, interval = c(min(BB), max(BB)))
  return(res$minimum)
}

# # Set parameters for the test
# xx <- mesh(7)  # Generate the mesh of design points
# results <- list()  # To store the results
# 
# # Run b_opt_MC for each combination of j and method
# for (j in JJ) {
#   cat(paste("\nTesting b_opt_MC for j =", j, "\n"))
# 
#   results[[j]] <- list()
# 
#   for (method in MM) {
#     cat(paste("Method:", method, "\n"))
# 
#     # Generate a sample for Monte Carlo integration
#     d <- length(xx[[1]])
#     unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
# 
#     # Start the timer
#     start_time <- Sys.time()
# 
#     # Run b_opt_MC and capture the result
#     b_opt_value_MC <- b_opt_MC(xx, j, method, unif_sample)
# 
#     # End the timer
#     end_time <- Sys.time()
# 
#     # Calculate elapsed time in minutes
#     elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
#     # Store the result and time taken
#     results[[j]][[method]] <- list(
#       b_opt_value = b_opt_value_MC,
#       time_taken = elapsed_time
#     )
# 
#     # Print the result and the time taken
#     cat(paste("b_opt_MC value for method", method, "with j =", j, ":", b_opt_value_MC, "\n"))
#     cat(paste("Time taken for method", method, "with j =", j, ":", round(elapsed_time, 2), "minutes\n"))
#   }
# }

######################################
## Optimal Bandwidth (grid version) ##
######################################

# Function to find the optimal bandwidth using a grid search for LSCV
b_opt_grid <- function(xx, j, method, return_LSCV = FALSE) {
  # Determine the number of cores and set the grid size
  num_cores <- detectCores() - 1
  grid_size <- 1 * num_cores
  
  # Generate grid points between the lower and upper bounds of BB
  b_grid <- seq(min(BB), max(BB), length.out = grid_size)
  
  # Initialize parallel cluster
  cl <- setup_parallel_cluster()
  
  # Parallelize computation for all b values on the grid
  LSCV_values <- parSapply(cl, b_grid, function(b) {
    LSCV(xx, b, j, method)
  })
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  # Find the index of the b that minimizes LSCV_values
  min_index <- which.min(LSCV_values)
  
  # Get the optimal b value and the corresponding LSCV value
  b_opt_value <- b_grid[min_index]
  min_LSCV_value <- LSCV_values[min_index]
  
  # Return the desired value(s) based on the return_LSCV argument
  if (return_LSCV) {
    return(min_LSCV_value)
  } else {
    return(b_opt_value)
  }
}

# Function to find the optimal bandwidth using a grid search for LSCV_MC
b_opt_MC_grid <- function(xx, j, method, unif_sample = NULL, return_LSCV_MC = FALSE) {
  # Generate unif_sample if not provided
  if (is.null(unif_sample)) {
    d <- length(xx[[1]])
    unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  }
  
  # Determine the number of cores and set the grid size
  num_cores <- detectCores() - 1
  grid_size <- 1 * num_cores
  
  # Generate grid points between the lower and upper bounds of BB
  b_grid <- seq(min(BB), max(BB), length.out = grid_size)
  
  # Initialize parallel cluster
  cl <- setup_parallel_cluster()
  
  # Parallelize computation for all b values on the grid
  LSCV_MC_values <- parSapply(cl, b_grid, function(b) {
    LSCV_MC(xx, b, j, method, unif_sample)
  })
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  # Find the index of the b that minimizes LSCV_MC_values
  min_index <- which.min(LSCV_MC_values)
  
  # Get the optimal b value and the corresponding LSCV_MC value
  b_opt_value <- b_grid[min_index]
  min_LSCV_MC_value <- LSCV_MC_values[min_index]
  
  # Return the desired value(s) based on the return_LSCV_MC argument
  if (return_LSCV_MC) {
    return(min_LSCV_MC_value)
  } else {
    return(b_opt_value)
  }
}

# # Set parameters for the test
# xx <- mesh(7)  # Generate the mesh of design points
# results_grid <- list()  # To store the results
# 
# # Run b_opt_MC_grid for each combination of j and method
# for (j in JJ) {
#   cat(paste("\nTesting b_opt_MC_grid for j =", j, "\n"))
# 
#   results_grid[[j]] <- list()
# 
#   for (method in MM) {
#     cat(paste("Method:", method, "\n"))
# 
#     # Generate a sample for Monte Carlo integration
#     d <- length(xx[[1]])
#     unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
# 
#     # Start the timer
#     start_time <- Sys.time()
# 
#     # Run b_opt_MC_grid and capture the result
#     b_opt_value_MC_grid <- b_opt_MC_grid(xx, j, method, unif_sample)
# 
#     # End the timer
#     end_time <- Sys.time()
# 
#     # Calculate elapsed time in minutes
#     elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
# 
#     # Store the result and time taken
#     results_grid[[j]][[method]] <- list(
#       b_opt_value = b_opt_value_MC_grid,
#       time_taken = elapsed_time
#     )
# 
#     # Print the result and the time taken
#     cat(paste("b_opt_MC_grid value for method", method, "with j =", j, ":", b_opt_value_MC_grid, "\n"))
#     cat(paste("Time taken for method", method, "with j =", j, ":", round(elapsed_time, 2), "minutes\n"))
#   }
# }

#######################################
## Plotting LSCV_MC vs b_opt_MC_grid ##
#######################################

# # Set parameters for the test
# xx <- mesh(7)  # Generate the mesh of design points
# results_grid <- list()  # To store the results
# 
# # Start the timer for the entire process
# start_time_total <- Sys.time()
# 
# # Loop over all methods and j values
# for (j in JJ) {
#   cat(paste("\nPlotting LSCV_MC vs b_opt_MC_grid for j =", j, "\n"))
# 
#   results_grid[[j]] <- list()
# 
#   for (method in MM) {
#     cat(paste("Method:", method, "\n"))
# 
#     # Start the timer for the current combination of j and method
#     start_time_individual <- Sys.time()
# 
#     # Generate a sample for Monte Carlo integration
#     d <- length(xx[[1]])
#     unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
# 
#     # Initialize parallel cluster for this combination
#     cl <- setup_parallel_cluster()
# 
#     # Determine the number of cores and set the grid size
#     num_cores <- detectCores() - 1
#     grid_size <- 1 * num_cores
# 
#     # Generate grid points between the lower and upper bounds of BB
#     b_values <- seq(min(BB), max(BB), length.out = grid_size)
# 
#     # Calculate LSCV_MC for different values of b in parallel
#     LSCV_MC_values <- unlist(parLapply(cl, b_values, function(b) {
#       LSCV_MC(xx, b, j, method, unif_sample)
#     }))
# 
#     # Use the b_opt_MC_grid function to find the optimal b value
#     b_opt_value_MC_grid <- b_opt_MC_grid(xx, j, method, unif_sample)
# 
#     # Stop the cluster after the parallel computation for this combination
#     stopCluster(cl)
# 
#     # Create a data frame for plotting
#     results_df <- data.frame(b = b_values, LSCV_MC = LSCV_MC_values)
# 
#     # Plot LSCV_MC values and indicate the optimal b
#     plot <- ggplot(results_df, aes(x = b, y = LSCV_MC)) +
#       geom_line(color = "blue") +
#       geom_vline(xintercept = b_opt_value_MC_grid, color = "red", linetype = "dashed") +
#       labs(
#         title = paste("LSCV_MC for j =", j, "and method =", method),
#         x = "b",
#         y = "LSCV_MC"
#       ) +
#       theme_minimal()
# 
#     # Save the plot as a PDF file
#     pdf(file = file.path(path, paste0("LSCV_MC_vs_b_opt_MC_grid_j", j, "_", method, ".pdf")))
#     print(plot)
#     dev.off()
# 
#     # End the timer for the current combination of j and method
#     end_time_individual <- Sys.time()
# 
#     # Calculate elapsed time in minutes for the current combination
#     elapsed_time_individual <- as.numeric(difftime(end_time_individual, start_time_individual, units = "mins"))
# 
#     # Store the result and time taken
#     results_grid[[j]][[method]] <- list(
#       b_opt_value = b_opt_value_MC_grid,
#       time_taken = elapsed_time_individual
#     )
# 
#     # Print the result and the time taken
#     cat(paste("b_opt_MC_grid value for method", method, "with j =", j, ":", b_opt_value_MC_grid, "\n"))
#     cat(paste("Time taken for method", method, "with j =", j, ":", round(elapsed_time_individual, 2), "minutes\n"))
#   }
# }
# 
# # End the timer for the entire process
# end_time_total <- Sys.time()
# 
# # Calculate total elapsed time in minutes
# elapsed_time_total <- as.numeric(difftime(end_time_total, start_time_total, units = "mins"))
# 
# # Print the total elapsed time in minutes
# cat("Total elapsed time: ", round(elapsed_time_total, 2), "minutes\n")

#####################################
## Integrated Squared Errors (ISE) ##
#####################################

ISE <- function(xx, j, method) {
  return(b_opt_grid(xx, j, method, return_LSCV = TRUE))
}

ISE_MC <- function(xx, j, method) {
  d <- length(xx[[1]])
  unif_sample <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1))
  return(b_opt_MC_grid(xx, j, method, unif_sample, return_LSCV_MC = TRUE))
}

###################################
## Real-data example (Section 6) ##
###################################

# Load necessary package
library(data.table)

# Specify the path to the CSV file
file_path <- "GEMAS.csv"  # CSV file name provided by the user

# Read the CSV file
data <- fread(file_path)

# Filter out rows with missing values in 'sand_norm', 'silt_norm', 'clay_norm', and 'pH_CaCl2'
filtered_data <- data[!is.na(sand_norm) & !is.na(silt_norm) & !is.na(clay_norm) & !is.na(pH_CaCl2), 
                      .(sand_norm, silt_norm, clay_norm, pH_CaCl2)]

# Convert the filtered data to a list of 4-vectors
vectors <- split(as.matrix(filtered_data), seq(nrow(filtered_data)))

# Extract xx and y from the given vectors
xx <- lapply(vectors, function(v) v[1:2] / 100)
y <- as.numeric(lapply(vectors, function(v) v[4]))

# Define the function to calculate the LOOCV for a given bandwidth b
LOOCV <- function(xx, y, b, method) {
  n <- length(xx)
  cv <- 0
  for (i in 1:n) {
    xx_minus_i <- xx[-i]
    y_minus_i <- y[-i]
    s <- xx[[i]]
    y_i <- y[[i]]
    y_hat_i <- hat_m(xx_minus_i, b, s, j, method, y_minus_i)
    cv <- cv + (y_i - y_hat_i) ^ 2
  }
  return(cv / n)
}

# Define the function to find the optimal bandwidth by minimizing LOOCV
b_opt_LOOCV <- function(xx, y, method) {
  cl <- setup_parallel_cluster()
  num_cores <- detectCores() - 1
  b_grid <- seq(0.005, 0.10, length.out = num_cores)  # Grid size matches the number of cores
  
  # Start timing the parallel portion
  start_time <- Sys.time()
  
  LOOCV_values <- parSapply(cl, b_grid, function(b) {
    LOOCV(xx, y, b, method)
  })
  
  # End timing the parallel portion
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat("Time taken for parallel computation: ", round(time_taken, 2), "minutes\n")
  
  stopCluster(cl)
  
  min_index <- which.min(LOOCV_values)
  b_opt_value <- b_grid[min_index]
  
  return(b_opt_value)
}

# Run the function to find the optimal bandwidth
method <- "LL"
optimal_b <- b_opt_LOOCV(xx, y, method)
print(paste("Optimal bandwidth (b):", optimal_b))

# Plot LOOCV as a function of b
b_values <- seq(0.005, 0.10, length.out = num_cores)
LOOCV_values <- sapply(b_values, function(b) LOOCV(xx, y, b, method))
results_df <- data.frame(b = b_values, CV = LOOCV_values)

# Save the plot to a PDF file
plot_file <- "loocv_b.pdf"
cairo_pdf(plot_file, width = 8, height = 5)

ggplot(results_df, aes(x = b, y = CV)) +
  geom_line(color = "blue") +
  labs(x = "b", y = "LOOCV") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )

dev.off()

print(paste("LOOCV plot saved to", plot_file))

