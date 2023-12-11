# This R script aims to demonstrate how to perform comparison-of-fit tests
# of nested models using the exact methods of algstat. More precisely, we 
# demonstrate by example how to perform a comparison-of-fit test between the 
# log-linear no23 and no3w models, and compare it to the theoretical 
# asymptotics.

# Note that the method used here to compute llm_exc_comp_p could be used to
# implement exact strategy of the DIF analysis paradigm where the DIF 
# classification step is augmented with a comparison-of-fit test. We did not 
# implement this, because this paradigm seems to result in minimal differences
# and because the calculation is somewhat intensive (especially if it was being
# iterated millions of times, like it would have had to have been in our 
# simulations!). 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(doParallel)

# Pull an algstat C++ function into the namespace
computeG2sCpp <- function(x, exp) {
    .Call('_algstat_computeG2sCpp', PACKAGE = 'algstat', x, exp)
}

# Choose a model
models_file <- "models.rds"
model <- generate_simulation_models(2015, models_file = models_file)

# Choose a distribution
# Note: Index 1 is no DIF (so should usually result in a high p-value for the 
# comparison), and indices 2 and 3 are uniform DIF (so should usually result in
# low p-values for the comparison). Indices 4 and 5 are non-uniform DIF. 
dist <- model$config$dist[[1]]

# Simulate a table
# Choose a large sample size to ensure that the exact and asymptotic 
# distributions match. 
t <- simulate_table(dist, 2000)
t

# Theoretical asymptotic df for the comparison of no23 vs no3w models
# This uses loglin, but it has a closed form formula: it should always equal 
# length(model$dimnames$responses) - 1. Note that, the way that the function 
# mle() is written, it returns NA for the df if the MLE doesn't exist,
# so a larger sample size should be chosen above in that case. 
mle_no23 <- mle(t, no23, llm)
mle_no3w <- mle(t, no3w, llm)
asy_df <- mle_no23$df - mle_no3w$df
asy_df

# Turn table into a vector
v <- vectorize(t, llm)

# MCMC parameters
# The larger these are, the better the job the Markov chain does of sampling 
# the fiber if iter is bigger, but the runtime for the comparison is roughly 
# (2 + iter) * (burn + (iter * thin)), which gets large very quickly...
iter <- 10000
burn <- 1000
thin <- 10

# Sample the fiber for no23 model (ie, the smaller model, which corresponds
# to the the bigger fiber)
sample_no23 <- metropolis(v, model$moves$llm_no23, iter, burn, thin)$steps |> 
    suppressMessages()

# Compute expected table in the no23 fiber
expected_no23 <- rowMeans(sample_no23)

# Takes as input a vectorized table in the no23 fiber and then computes the 
# corresponding expected table under the no3w model
expected_no3w <- function(v, method = "exact") {
    if (method != "hyrbid") {
        # Default: Use an exact method to compute expected table (as the mean 
        # table in the fiber)
        out <- algstat::metropolis(v, model$moves$llm_no3w, 
                                   iter, burn, thin)$steps |> 
            suppressMessages() |> 
            rowMeans()
    } else {
        # If method == "hybrid", use the MLE as the expected table instead of
        # the mean table of the fiber
        t <- unvectorize(v, dimnames(t), llm)
        out <- mle(t, no3w, llm)$expected |> 
            vectorize(llm)
    }
    out
}

# Set up cluster for parallel computation
num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)
registerDoParallel(cl = cl)
clusterEvalQ(cl, library(algstat))

# Compute MCMC distribution of G using the parallel cluster
# This takes a little bit of time (on the order of 4-6 minutes on our machines,
# with the parameters that appear above). 
start_time <- Sys.time()
G <- foreach(i = 1:ncol(sample_no23), .combine = c) %dopar% {
    sample_no23[, i] |> 
        expected_no3w() |> 
        matrix() |> 
        computeG2sCpp(expected_no23)
}
Sys.time() - start_time
stopCluster(cl)

# Compute the observed value of the G test statistic (UMVUE-based method)
G_obs_exc <- expected_no3w(v) |> 
    matrix() |> 
    computeG2sCpp(expected_no23)

# Compute the observed value of the G test statistic
# This can be computed as mle_no23$lrt - mle_no3w$lrt, or equivalently, if we
# use Agresti eqn (10.3), it can be computed using the following, which
# may appear more analogous to our calculations above
G_obs_asy <- computeG2sCpp(matrix(vectorize(mle_no3w$expected, llm)), 
                           vectorize(mle_no23$expected, llm))

# Compare the outputs
G_obs_exc
G_obs_asy

# Compute exact p-value
llm_exc_comp_p <- mean(G >= G_obs_exc)

# Compute asymptotic p-values (two ways)
llm_asy_comp_p <- pchisq(G_obs_asy, asy_df, lower.tail = FALSE)

# Compare the outputs
llm_exc_comp_p
llm_asy_comp_p

# Make plot of the MCMC distribution (black) alongside the theoretical 
# asymptotic distribution (purple). The observed values of the G test statistic
# are indicated by vertical lines (faded black and faded pink)
faded <- 0.3
asy_color <- "#C699E7"
tibble(G) |> 
    ggplot(aes(G)) + 
    geom_vline(xintercept = G_obs_asy, color = asy_color, alpha = faded) + 
    stat_function(fun = dchisq, args = list(df = asy_df), color = asy_color) + 
    geom_vline(xintercept = G_obs_exc, alpha = faded) + 
    geom_density() + 
    theme_minimal()
