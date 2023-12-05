library(tidyverse)
library(detectseparation)
library(algstat)

################################################################################
# Miscellaneous
################################################################################

# Pull an algstat C++ functions into the namespace
computeUProbsCpp <- function(x) {
    .Call("_algstat_computeUProbsCpp", PACKAGE = "algstat", x)
}

# Facet specifications of models
no23 <- list(c(1,2), c(1,3))
no3w <- list(c(1,2), c(1,3), c(2,3))
full <- list(c(1,2,3))

# Other
FT <- c(FALSE, TRUE)

################################################################################
# Compute Markov moves
################################################################################

# Compute configuration matrices for relevant models. The input should be
# a list of numeric vectors, where the first vector is the list of ability
# levels, the second is the list of groups levels, and the third is the list
# of response levels. The output is a named list of matrices, where the names
# are llm_no23, llm_no3w, lrm_no23, lrm_no3w, lrm_full.
configuration_matrices <- function(dimnames) {
    # Log-linear models
    d <- map_dbl(dimnames, length)
    config <- list(llm_no23 = no23, llm_no3w = no3w) |> 
        map(\(m) unname(algstat::hmat(d, m)))
    
    # Logistic regression models
    if (d[3] == 2) {
        A <- expand_grid(1, ability = dimnames[[1]], group = dimnames[[2]]) |> 
            mutate(interaction = ability * group) |> 
            t()
        config$lrm_no23 <- algstat::lawrence(A[1:2,])
        config$lrm_no3w <- algstat::lawrence(A[1:3,])
        config$lrm_full <- algstat::lawrence(A)
    } else {
        config$lrm_no23 <- NULL
        config$lrm_no3w <- NULL
        config$lrm_full <- NULL
    }
    
    # Return
    config
}

# Helper function for going from configuration matrices to Markov bases. All
# it does on top of latte:markov is return NULL if the matrix is NULL. 
markov_wrapper <- function(A) {
    if (is_null(A)) {
        out <- NULL
    } else {
        out <- latte::markov(A, p = "arb")
    }
    out
}

# Compute Markov moves for relevant models. The input should be
# a list of numeric vectors, where the first vector is the list of ability
# levels, the second is the list of groups levels, and the third is the list
# of response levels. The output is a named list of matrices whose columns
# are the moves, where the names are llm_no23, llm_no3w, lrm_no23, lrm_no3w, 
# lrm_full.
markov_moves <- function(dimnames) {
    map(configuration_matrices(dimnames), markov_wrapper)
}

################################################################################
# Generate simulation models
################################################################################

# Helper function to compute a joint distribution of ability, group, response
# using Bock's model, where the conditional probability 
# P[R = r | A = a, G = g] 
# is proportional to 
# exp(u(r,g)*a + v(r,g)) 
# where 
# u(r,g) = u(r) + u_delta(r)*g 
# v(r,g) = v(r) + v_delta(g)*g
# for real numbers u(r), u_delta(r), v(r), v_delta(r).
# WLOG, the vectors u, v, u_delta, v_delta can be taken to be in the orthogonal
# complement of the span of (1, 1, ..., 1), though the function below does not
# insist on this. (u, v) specifies the distribution of responses for group 0. 
# dif = (u_delta, v_delta) specifies how the distribution changes between the 
# groups. Function returns a three-way array representing the joint 
# distribution of (A, G, R).  
# The input base must be a list with the following named elements:
# * base$dimnames: a named list where
#       - base$dimnames$ability is the list of numerical ability levels
#       - base$dimnames$group is the list of numerical group levels 0:1
#       - base$dimnames$response is the list of response levels
# * base$group_dist: a two-element list giving the probabilies of groups 0, 1
# * base$ability_dist: a two-column matrix giving the probabilites of ability
#       conditioned on group.
# * base$uv: the vector c(u, v) described above.
# All other named elements of base are ignored. 
# Caution: If uv or dif contain values that are too large and cause floating
# point overflow, the output to this function will contain NaN values!
joint_distribution <- function(base, dif) {
    # Extract vectors u, v, u_delta, v_delta
    d <- map_dbl(base$dimnames, length)
    K <- d[3]
    if (length(base$uv) != 2*K || length(dif) != 2*K) {
        stop("Mismatch of parameter lengths.")
    }
    u <- base$uv[1:K]
    v <- base$uv[(K+1):(2*K)]
    u_delta <- dif[1:K]
    v_delta <- dif[(K+1):(2*K)]
    
    # Helper function for computing unnormalized conditional probabilities
    lambda <- function(a, g) {
        exp((u + u_delta * g) * a + (v + v_delta * g))
    }
    
    # Array of unnormalized conditional probabilities of R given (A, G)
    alpha <- expand_grid(ability = base$dimnames$ability, 
                         group = base$dimnames$group) |>
        mutate(alpha = map2(ability, group, lambda)) |> 
        unnest(alpha) |> 
        pull(alpha) |> 
        vec2tab(d)
    dimnames(alpha) <- base$dimnames
    
    # Array of conditional probabilities of R given (A, G)
    r_given_ag <- sweep(alpha, c(1,2), marginSums(alpha, c(1,2)), `/`)
    
    # Array of joint distribution of (A, G)
    ag <- sweep(base$ability_dist, 2, base$group_dist, `*`)
    
    # Array of joint distribution of (A, G, R)
    sweep(r_given_ag, c(1, 2), ag, `*`)
}

# Helper function that takes as input a joint distribution of three variables
# specified as a three-way array, and returns the conditional distribution of 
# the third variable given the first two again as a three-way array.
conditional <- function(dist) {
    sweep(dist, c(1,2), marginSums(dist, c(1,2)), `/`)
}

# Helper function to compute the Rényi alpha-divergence from q to p, where p 
# and q are numerical vectors of equal length containing probabilities that 
# sum to 1. Defaults to alpha = 1, which corresponds to Kullback-Leibler 
# divergence.
renyi <- function(p, q, alpha = 1) {
    if (alpha < 0) {
        stop("Alpha parameter must be positive.")
    } else if (length(p) != length(q) || 
        any(p < 0) || 
        any(q < 0) || 
        abs(sum(p) - 1) > .Machine$double.eps || 
        abs(sum(q) - 1) > .Machine$double.eps
        ) {
        stop("Rényi divergence can only be computed on comparable probability vectors!")
    }
    if (alpha == 0) {
        out <- keep(q, p > 0) |> 
            sum() |> 
            log()
        out <- -out
    } else if (alpha == 1) { 
        out <- map2_dbl(p, q, \(x, y) x * log(x / y)) |> 
            sum()
    } else if (alpha == Inf) {
        out <- map2_dbl(p, q, \(x, y) log(x / y)) |> 
            max()
    } else {
        out <- map2_dbl(p, q, \(x, y) p^alpha / q^(alpha - 1)) |> 
            sum() |> 
            log()
        out <- 1/(alpha - 1) * out
    }
    out
}

# Helper function to compute the KL divergence from q to p, where p and q are 
# numerical vectors of equal length containing probabilities that sum to 1. 
kl <- function(p, q) {
    renyi(p, q, alpha = 1)
}

# Helper function that, from a joint distribution of ability, group, and 
# response, computes a vector of statistical distances: for each level a of 
# ability, compute the statistical distance from conditional distribution of 
# response given A = a and G = 0 to that given A = a and G = 1. Defaults to 
# using Kullback-Leibler divergence. 
dif_size_vector <- function(dist, distance = kl) {
    r_given_ag <- conditional(dist)
    1:(dim(dist)[1]) |> 
        map_dbl(\(a) kl(r_given_ag[a,2,], r_given_ag[a,1,]))
}

# Helper function that computes a statistically meaningful "size of DIF" 
# in a joint distribution of ability, group, and response by taking the 
# L^p norm of the distance vector. Defaults to Kullback-Leibler divergence and 
# L^infty norm. 
dif_size <- function(dist, vector = dif_size_vector, p = Inf) {
    lpnorm(vector(dist), p)
}

# Helper function that normalizes a nonzero vector in R^n to unit length
# with respect to the L^p norm. Defaults to p = 2. 
normalize <- function(x, p = 2) {
    c = lpnorm(x, p)
    if (c == 0) {
        out <- x
    } else {
        out <- x/c
    }
    out
}

# Helper function that akes a vector of length K-1 and outputs the vector of 
# length K in the orthogonal complement of (1, 1, ..., 1) whose coordinate 
# vector with respect to the basis (-1, 1, 0, ..., 0), (-1, 0, 1, ..., 0), ..., 
# (-1, 0, 0, ..., 1) is the input vector. Put more simply, it outputs 
# the vector c(-sum(v), v).
to_constants_perp <- function(v) {
    # A <- rbind(rep(-1, length(v)), diag(length(v)))
    # as.vector(A %*% v)
    c(-sum(v), v)
}

# Main function: Uses an integer seed to generate a simulation model. There
# are two groups, 0 and 1. It draws from a symmetric beta distribution (with
# parameter being the input value beta) to extract the probability of group 0.
# It draws from a uniform_distribution on 1:max_dim to generate a set of 
# ability levels (integers from 0 up till one minus the number drawn). 
# For each group, the distribution of ability given group is binomial with
# parameter being drawn from a beta distribution again. 
# It then draws an integer from 1:max_dim, with the probability of 1 being
# dichotomous_prob and the others being equally likely, and uses the integers
# from 0 up to one minus the integer drawn as response levels. 
# It then draws a unit vector (u, v) as described in the comments to the helper 
# function joint_distribution above, where u and v are each in the orthogonal 
# complement of the span of (1, 1, ..., 1). 
# It then draws a uniform DIF vector (0, v_delta) of norm 1, and then a 
# separate a nonuniform DIF vector (u_delta, v_delta). 
# It finally uses these data to produce a tibble of several joint distributions
# of (A, G, R) using the joint_distribution helper function above:
# * One with no DIF.
# * Several with uniform DIF, with DIF vector being scalar multiples of the 
#       chosen uniform DIF vector as the DIF vector of the model. 
# * Several with nonuniform DIF, with DIF vector being scalar multiples of the 
#       using scalar multiples of the chosen nonuniform DIF
# If include_moves is TRUE, it also computes relevant Markov moves for these
# models and includes the results in the output. 
# The output is a named list containing all of the "movable parts" above. 
# If models_file is not NULL, it attempts to read from models_file before
# doing the computation, and it writes to that file afterwards to save results.
generate_simulation_models <- function(seed, 
                                       beta = 10, 
                                       max_dim = 7,
                                       dichotomous_prob = 0.5, 
                                       norms = c(1, 2),
                                       include_moves = TRUE, 
                                       models_file = "models.rds"
                                       ) {
    # Read file
    if (!is_null(models_file) && file.exists(models_file)) {
        models <- readRDS(models_file)
    } else {
        models <- list()
    }
    
    # Check file for model
    for (m in models) {
        if (seed == m$seed) {
            return(m)
        }
    }
    
    # Otherwise, generate model
    # Set seed
    set.seed(seed)
    
    # Generate distribution of groups
    group_levels <- 0:1
    J <- length(group_levels)
    group_dist <- rbeta(1, beta, beta) |>
        (\(p) c(p, 1-p))()
    names(group_dist) <- group_levels
    
    # Generate ability levels
    I <- sample(2:max_dim, 1)
    ability_levels <- 0:(I-1)
    
    # Generate distributions of ability conditioned on group
    # Each column is the ability distribution of a fixed group
    ability_dist <- expand_grid(p = rbeta(2, beta, beta), ability_levels) |>
        mutate(q = dbinom(ability_levels, I-1, p)) |>
        pull(q) |>
        matrix(I, J, 
               dimnames = list(ability = ability_levels, group = group_levels))
    
    # Generate response levels, with probability of dichotomous being given
    num_polytomous_options <- length(3:max_dim)
    prob <- c(dichotomous_prob, 
              rep((1-dichotomous_prob) * 1/num_polytomous_options, 
                  num_polytomous_options))
    K <- sample(2:max_dim, 1, prob = prob)
    response_levels <- 0:(K-1)
    
    # Dimension vector
    dimnames <- list(ability = ability_levels, 
                     group = group_levels, 
                     response = response_levels)
    
    # Response distribution for group 0
    u <- runif(K-1, -1, 1) |> 
        to_constants_perp()
    v <- runif(K-1, -1, 1) |>
        to_constants_perp()
    uv <- c(u, v) |> 
        normalize()
    
    # Make base list for passing to joint_distribution
    base <- list()
    base$dimnames <- dimnames
    base$group_dist <- group_dist
    base$ability_dist <- ability_dist
    base$uv <- uv
    
    # No DIF vector
    u_delta <- rep(0, K)
    v_delta <- rep(0, K)
    none <- c(u_delta, v_delta)
    
    # Uniform DIF vector
    u_delta <- rep(0, K)
    v_delta <- runif(K-1, -1, 1) |> 
        to_constants_perp()
    unif <- c(u_delta, v_delta) |> 
        normalize()
    
    # Non-uniform DIF vector
    u_delta <- runif(K-1, -1, 1) |>
        to_constants_perp()
    v_delta <- runif(K-1, -1, 1) |> 
        to_constants_perp()
    nonunif <- c(u_delta, v_delta) |> 
        normalize()
    
    # Make tibble of all DIF configurations
    config <- tibble(type = c("unif", "nonunif"), 
                     dif = list(unif, nonunif)) |> 
        expand_grid(norm = norms) |> 
        mutate(dif = map2(dif, norm, \(x, s) s * x))
    config <- tibble(type = "none", dif = list(none), norm = 0) |> 
        bind_rows(config) |> 
        mutate(dist = map(dif, \(x) joint_distribution(base, x))) |> 
        mutate(seed = seed) |> 
        select(seed, type, norm, dif, dist) |> 
        mutate(kl = round(map_dbl(dist, dif_size), 3))
    
    # Unset seed
    set.seed(NULL)
    
    # Make output
    out <- list()
    out$seed <- seed
    out$dimnames <- dimnames
    out$is_dichotomous <- (K == 2)
    out$group_dist <- group_dist
    out$ability_dist <- ability_dist
    out$uv <- uv
    out$unif <- unif
    out$nonunif <- nonunif
    out$config <- config
    
    # Calculate moves
    if (include_moves) {
        start_timer <- Sys.time()
        out$moves <- markov_moves(dimnames)
        out$moves_computation_time <- Sys.time() - start_timer
    } else {
        out$moves <- NULL
        out$moves_computation_time <- NULL
    }
    
    # Save to models file
    models <- append(models, list(out))
    reorder <- map_dbl(models, \(m) m$seed) |> 
        order()
    models <- models[reorder]
    if (!is_null(models_file)) {
        saveRDS(models, models_file)
    }
    
    # Return
    out
}

################################################################################
# Generate simulation data
################################################################################

# Generate a three-way table of observation counts following the three-way 
# input distribution, where the number of observations is the given 
# sample_size. Insist on the table passing the given predicate, but try at 
# most max times. 
simulate_table <- function(dist, 
                           sample_size, 
                           predicate = \(t) TRUE, 
                           max = 1e3
                           ) {
    count <- 0
    repeat {
        # Generate table
        t <- rmultinom(1, sample_size, algstat::tab2vec(dist)) |>
            as.vector() |> 
            algstat::vec2tab(dim(dist)) |> 
            as.table()
        dimnames(t) <- dimnames(dist)
        
        # Check whether we keep generating tables
        count <- count + 1
        if (predicate(t)) {
            break
        } else if (count > max) {
            stop("Unable to pass predicate within max iterations.")
        }
    }
    
    # Return
    t
}

################################################################################
# MLE-related functions
################################################################################

# Helper function to create a list of toggles for models. Raises error if model 
# is not one of the facet specifications no23, no3w, full or one of those 
# four-letter strings. Raises error if family is not one of the three letter 
# strings llm or lrm. 
toggler <- function(model, family) {
    # Model toggles
    out1 <- list(no23 = FALSE, no3w = FALSE, full = FALSE)
    if (is.character(model) && length(model) == 1) {
        out1 <- map2(out1, names(out1), \(x, y) y == model)
    } else {
        out1 <- map2(out1, list(no23, no3w, full), \(x, y) setequal(y, model))
    }
    
    if (sum(unlist(out1)) != 1) {
        stop("Model selection error!")
    }
    
    # Family toggles
    out2 <- list(llm = FALSE, lrm = FALSE)
    out2 <- map2(out2, names(out2), \(x, y) y == family)
    
    if (sum(unlist(out2)) != 1) {
        stop("Family selection error!")
    }
    
    append(out1, out2)
}

# Helper function that reorganize table into a tibble that's useful for 
# logistic regression For each ability-group pair, it records the proportion 
# of correct responses in a column called p, and the total number of subjects 
# in a column called n. The table must be dichotomous for this to make sense,
# but this condition is not checked. 
reorganize <- function(t) {
    as_tibble(t) |>
        mutate(across(everything(), as.numeric)) |>
        group_by(ability, group) |> 
        summarize(p = sum(keep(n, response == 1))/sum(n), n = sum(n), 
                  .groups = "drop")
}

# Helper function that reverse the process of the reorganize function above.
# Caution: unreorganize(reorganize(t)) may differ from t by small amounts
# due to floating point arithmetic. This is a useful for computing
# expected counts under logistic regression MLEs.
unreorganize <- function(u) {
    dimnames <- list(ability = unique(u$ability), 
                     group = unique(u$group), 
                     response = 0:1)
    t <- rbind((1-u$p) * u$n, u$p * u$n) |> 
        algstat::vec2tab(map_dbl(dimnames, length)) |> 
        as.table()
    dimnames(t) <- dimnames
    t
}

# Helper function that adds two ability levels, one lower than all existing
# ones and one higher, with no observations, and returns the result as a
# tibble of counts; useful for checking MLE existence for logistic regression
# no23 and no3w models. dimnames must be map(dimnames(t), as.numeric), but
# this condition is not checked. 
tibble_and_append <- function(t, dimnames) {
    w <- as_tibble(t) |> 
        mutate(across(everything(), as.numeric))
    to_add <- c(min(dimnames$ability)-1, max(dimnames$ability)+1)
    expand_grid(ability = to_add, 
                group = dimnames$group, 
                response = dimnames$response,
                n = 0) |> 
        bind_rows(w) |> 
        arrange(ability, group, response)
}

# Helper function to check for "cross zeros." Used in checking for MLE 
# existence for the no23 and no3w models for logistic regression.
cross_zeros <- function(w, a0, a1) {
    x <- w |> 
        mutate(q = group * (a1 - a0) - (ability - a0)) |> 
        filter(q != 0) |> 
        mutate(side = (q > 0)) |> 
        group_by(side, response) |> 
        summarize(n = sum(n), .groups = "drop") |> 
        pull(n)
    if (length(x) != 4) {
        stop("Error!")
    }
    (x[1] == 0 && x[4] == 0) || (x[2] == 0 && x[3] == 0)
}

# Helper function for checking if a 2x2x2 has diagonal zeros.
# Used for MLE existence check for log-linear no3w model.
diagonal_zeros <- function(t) {
    (t[1,1,1] == 0 && t[2,2,2] == 0) ||
        (t[1,1,2] == 0 && t[2,2,1] == 0) ||
        (t[1,2,1] == 0 && t[2,1,2] == 0) ||
        (t[1,2,2] == 0 && t[2,1,1] == 0)
}

# Helper function for collapsing a IxJxK table into a 2x2x2 table.
# collapsing is a list of three logical vectors, one of length I, 
# one of length J, and one of length K. None of the three vectors 
# should be all TRUE or all FALSE, but this condition is not checked!
collapse <- function(t, collapsing) {
    out <- expand_grid(a = FT, b = FT, c = FT) |> 
        mutate(n = pmap_dbl(list(a, b, c), 
                            \(x, y, z) sum(t[which(collapsing[[1]] == x), 
                                             which(collapsing[[2]] == y), 
                                             which(collapsing[[3]] == z)]))) |> 
        pull(n) |> 
        algstat::vec2tab(c(2,2,2))
    if (!is_null(dimnames(t))) {
        dimnames(out) <- map2(dimnames(t), 
                              collapsing, 
                              \(x, y) map_chr(FT, 
                                              \(z) str_flatten(keep(x, y == z))
                                              )
                              )
    }
    out
}

# Helper function to generate a list of all nontrivial binary partitions 
# of a set of of length n. A binary partition is specified as a list of n 
# booleans indicating which elements are *not* in the same class as the first.
# Output is a matrix, where each row is a binary partition. 
# There are 2^(n-1)-1 rows. Used for MLE existence check for llm no3w.
binary_partitions <- function(n) {
    rep(list(FT), n) |> 
        expand.grid() |> 
        filter(Var1 == FALSE) |> 
        slice(-1) |> 
        as.matrix() |> 
        unname()
}

# Helper function for checking if any 2x2x2 collapsing of a IxJxK table has
# diagonal zeroes. Used for MLE existence check for llm no3w.
collapsing_has_diagonal_zeros <- function(t) {
    partitions <- unname(dim(t)) |> 
        map(binary_partitions)
    lambda <- function(x, y, z) {
        collapsing <- list(partitions[[1]][x,], 
                           partitions[[2]][y,], 
                           partitions[[3]][z,])
        collapse(t, collapsing)
    }
    df <- map(partitions, \(x) seq((dim(x))[1])) |> 
        expand.grid() |>
        as_tibble() |> 
        mutate(collapsing = pmap(list(Var1, Var2, Var3), lambda)) |>
        mutate(diagonal_zeros = map_lgl(collapsing, diagonal_zeros))
    any(df$diagonal_zeros)
}

# Check if table t has MLE for the model (no23 or no3w or full) and family (llm
# or lrm) specified by the toggle. This toggle must be a valid list returned
# by the toggler function, but this condition is not checked. The function
# attempts to short circuit calculations whenever possible, especially for 
# logistic regression since the general method there is computationally 
# expensive. To skip that short circuiting, pass avoid_lp = FALSE. 
# This option is mostly there for testing purposes.
mle_exists <- function(t, toggle, avoid_lp = TRUE) {
    # No sampling zeros; sufficient, but not necessary, for any MLE existence
    # It is also necessary for MLE existence for full log-linear model, and
    # the full logistic regression model if the table is 2x2x2
    out <- all(t > 0)
    if (out || 
        (toggle$full && (toggle$llm || (toggle$lrm && all(dim(t) == 2))))) {
        return(out)
    }
    
    if (toggle$llm) {
        # Log-linear models
        # First, check that line sums n_{a,g,+} and n_{a,+,r} are nonzero
        # Necessary and sufficient for MLE existence for no23
        # Sufficient, but not necessary for MLE existence for no3w
        out <- all(marginSums(t, c(1,2)) > 0, marginSums(t, c(1,3)) > 0)
        if (!out || toggle$no23) {
            return(out)
        } 
        
        # For the no3w model, first check that all line sums are nonzero
        # This is necessary, but not sufficient, for MLE existence
        out <- all(marginSums(t, c(2,3)) > 0)
        if (!out) {
            return(out)
        }
        
        # Then check for non-existence of diagonal pattern via 2x2x2 
        # collapsings à la Eriksson et al 2006
        # Check is somewhat expensive (exponential in the dimensions of the 
        # table) and it is short circuited
        if (any(dim(t) == 2)) {
            return(!collapsing_has_diagonal_zeros(t))
        } 
        
        # Otherwise, we are in a table where no dimension is 2, in which case
        # a characterization of MLE existence is unknown
        stop("MLE existence for log-linear no3w model is unimplemented when none of the dimensions is 2!")
    } else {
        # Logistic regression models
        # First, check that we are in a dichotomous situation
        if (dim(t)[3] != 2) {
            stop("Response must be dichotomous for logistic regression models.")
        }
        
        # First check that line sums n_{a,g,+} are nonzero 
        # Necessary, but not sufficient, for MLE existence
        if (any(marginSums(t, c(1,2)) == 0)) {
            return(FALSE)
        }
        
        # Check some easier cases of separation
        # This short circuits the general method via linear programming 
        # implemented in the detectseparation package, which is computationally 
        # expensive when the sample size is large
        if (avoid_lp) {
            # First, check that there is no separation on ability
            dimnames <- map(dimnames(t), as.numeric)
            w <- tibble_and_append(t, dimnames)
            out <- all(map_lgl(dimnames$ability, \(a) !cross_zeros(w, a, a)))
            
            # This is enough to check for no23, and if it fails, the MLE for
            # no3w or full cannot exist either
            if (!out || toggle$no23) {
                return(out)
            }
            
            # If there are two groups...
            if (dim(t)[2] == 2) {
                # First, check that line sums n_{+,g,r} are nonzero
                # This is redundant, but maybe it is useful to short circuit
                # these computations whenever possible...
                if (any(marginSums(t, c(2,3)) == 0)) {
                    return(FALSE)
                }
                
                # Check for separation along any pair of ability levels
                # This is enough to check for no3w, and if it fails, the MLE for
                # no3w cannot exist either
                ability_pairs <- expand_grid(a0 = dimnames$ability, 
                                             a1 = dimnames$ability) |> 
                    filter(a0 != a1) |> 
                    mutate(x = !map2_lgl(a0, a1, 
                                         \(a0, a1) cross_zeros(w, a0, a1)))
                out <- all(ability_pairs$x)
                if (!out || toggle$no3w) {
                    return(out)
                }
            }
        }
        
        # If the short circuit failed, check via linear programming 
        # à la Konis 2007. This computation is expensive for large sample 
        # sizes, but if the above went through, it should only be called when
        # there are more than 2 groups (which we do not implement) or if
        # we're interested in the full lrm model and the no3w MLE exists
        # (and we aren't currently using anyting about the full lrm model).
        # First, convert table t to a tidy tibble
        v <- as_tibble(t) |>
            mutate(across(everything(), as.numeric)) |> 
            uncount(n)
        if (toggle$no23) {
            out <- !glm(response ~ ability, 
                        data = v, 
                        family = binomial("logit"), 
                        method = "detect_separation")$outcome
        } else if (toggle$no3w) {
            out <- !glm(response ~ ability + group, 
                        data = v, 
                        family = binomial("logit"), 
                        method = "detect_separation")$outcome
        } else if (toggle$full) {
            out <- !glm(response ~ ability + group + ability:group, 
                        data = v, 
                        family = binomial("logit"), 
                        method = "detect_separation")$outcome
        }
        return(out)
    }
}

# Check if table has any sampling zeros.
has_no_zeros <- function(t) {
    mle_exists(t, list(full = TRUE, llm = TRUE))
}

# Check if table satisfies the heuristic of having sufficiently many expected
# counts large enough, where "sufficiently many" means at least 
# min_percent of the total number (default: 80%) and "large enough" means 
# greater than or equal to min_number (default: 5).
heuristic_check <- function(expected, min_percent = 0.8, min_number = 5) {
    percent_large <- sum(expected >= min_number) / length(expected)
    percent_large >= min_percent
}

# Main function: collect all relevant information related to the MLE of the 
# input table for the input model (no23, no3w, or full) for the input family 
# ("llm" for log-linear models and "lrm" for logistic regression models).
# The input threshold is used to determine the type of MLE existence checking
# If it is <= 0, the exact test of mle_exists is used. If it is > 0, an 
# approximate test using the fit returned by loglin or glm is used to determine 
# if the fitted distribution lies in the interior of the probability simplex,
# with all entries being at least the threshold.  
# The output collects the following information:
# * mle_exists: TRUE if the MLE exists
# * expected: three-way table of expected counts
# * heuristic: TRUE if expected counts pass the heuristic check
# * lrt: G test statistic
# * df: degrees of freedom
# * p: tail area part lrt in chi-square distribution with degrees of freedom df
mle <- function(t, model, family, threshold = 0) {
    # Make toggle and error check
    toggle = toggler(model, family)
    
    # Do an exact check of MLE existence
    if (threshold <= 0) {
        mle_exists <- mle_exists(t, toggle)
    } else {
        mle_exists <- TRUE
    }
    
    # Calculate MLEs numerically
    if (mle_exists) {
        # Get expected counts, fitted probabilities, G test statistic, and df
        if (toggle$llm) {
            # Log-linear models
            # Use loglin to compute MLE
            if (toggle$no23) {
                x <- loglin(t, no23, fit = TRUE)
            } else if (toggle$no3w) {
                x <- loglin(t, no3w, fit = TRUE)
            } else if (toggle$full) {
                x <- loglin(t, full, fit = TRUE)
            }
            expected <- x$fit
            lrt <- x$lrt
            df <- x$df
        } else {
            # Logistic regression
            # First, check that we are in a dichotomous situation
            if (dim(t)[3] != 2) {
                stop("Response must be dichotomous for logistic regression models.")
            }
            # Use glm to compute MLE
            u <- reorganize(t)
            if (toggle$no23) {
                x <- glm(p ~ ability, 
                         weights = n, 
                         data = u, 
                         family = binomial, 
                         na.action = na.pass)
            } else if (toggle$no3w) {
                x <- glm(p ~ ability + group, 
                         weights = n, 
                         data = u, 
                         family = binomial, 
                         na.action = na.pass)
            } else if (toggle$full) {
                x <- glm(p ~ ability + group + ability:group, 
                         weights = n, 
                         data = u, 
                         family = binomial, 
                         na.action = na.pass)
            }
            expected <- u |> 
                mutate(p = x$fitted.values) |> 
                unreorganize()
            lrt <- x$deviance
            df <- x$df.residual
        }
        
        # Check heuristic, compute MLE distribution, and p-value
        heuristic <- heuristic_check(expected)
        mle_dist <- expected / sum(t)
        p <- pchisq(lrt, df, lower.tail = FALSE)
    }
    
    # Do a numerical check for MLE existence
    if (threshold > 0) {
        mle_exists <- all(marginSums(t, c(1,2)) > 0, mle_dist > threshold)
    }
    
    # If the MLE doesn't exist, make empty output
    if (!mle_exists) {
        heuristic <- FALSE
        mle_dist <- NA
        expected <- NA
        lrt <- NA
        df <- NA
        p <- NA
    }
    
    # Collect output
    out <- list()
    out$mle_exists <- mle_exists
    out$heuristic <- heuristic
    out$mle_dist <- mle_dist
    out$expected <- expected
    out$lrt <- lrt
    out$df <- df
    out$p <- p
    
    # Return
    out
}

################################################################################
# Task function for each table
################################################################################

# Strings to indicate DIF test results:
# * fail = test could not determine if there is DIF
# * none = test concluded there is no DIF
# * unclassifiable = test concluded there is no DIF but could not classify
# * uniform = test concluded there is uniform DIF
# * nonuniform = test concluded there is nonuniform DIF
dif_test_results <- list(
    fail = "fail", 
    none = "none", 
    unclassifiable = "unclass",
    uniform = "unif", 
    nonuniform = "nonunif"
)

# Conduct a battery of DIF tests: the main one, the swapped variant, and the 
# variant with comparison (if comp_p is provided). 
dif_tests <- function(no23_abort, 
                      no23_p, 
                      no3w_abort,
                      no3w_p, 
                      comp_p = NULL, 
                      alpha = 0.05
) {
    # Main DIF test, using the paradigm described in the paper
    if (no23_abort) {
        main <- dif_test_results$fail
    } else if (no23_p >= alpha) {
        main <- dif_test_results$none
    } else if (no3w_abort) {
        main <- dif_test_results$unclassifiable
    } else if (no3w_p >= alpha) {
        main <- dif_test_results$uniform
    } else {
        main <- dif_test_results$nonuniform
    }
    
    # DIF test with comparison
    if (is_null(comp_p)) {
        # If comp_p is NULL, return NA to indicate that the test was not run
        comp <- NA
    } else if (main == dif_test_results$uniform && comp_p >= alpha) {
        # The only case when this variant differs from the main variant is 
        # when the main variant concluded uniform DIF, but the comparison
        # has a big p-value; in this case, DIF is unclassifiable. 
        comp <- dif_test_results$unclassifiable
    } else {
        comp <- main
    }
    
    # Swapped DIF test, using the variant paradigm from the paper
    if (no23_abort || no3w_abort) {
        swap <- dif_test_results$fail
    } else if (no3w_p < alpha) {
        swap <- dif_test_results$nonuniform
    } else if (no23_p < alpha) {
        swap <- dif_test_results$uniform
    } else {
        swap <- dif_test_results$none
    }
    
    tibble(main, comp, swap)
}

# Main function: Performs tasks for each table. Information is returned as
# a tibble with a single row. It checks a number of predicates on the table
# such as MLE existence and sample size heuristics. It records a number of
# relevant p-values, and also the results of a number of DIF tests. 
table_task <- function(t, 
                       moves,
                       iter = 10000L, 
                       burn = 1000L, 
                       thin = 10L,
                       alpha = 0.05
                       ) {
    # Log-linear models
    
    # Asymptotic LLM
    llm_asy_no23 <- mle(t, no23, "llm")
    llm_asy_no3w <- mle(t, no3w, "llm")
    llm_asy_comp_p <- pchisq(llm_asy_no23$lrt - llm_asy_no3w$lrt, 
                             llm_asy_no23$df - llm_asy_no3w$df, 
                             lower.tail = FALSE)
    
    llm_asy <- dif_tests(no23_abort = !llm_asy_no23$mle_exists,
                         no23_p = llm_asy_no23$p,
                         no3w_abort = !llm_asy_no3w$mle_exists,
                         no3w_p = llm_asy_no3w$p,
                         comp_p = llm_asy_comp_p,
                         alpha = alpha) |> 
        rename_with(\(x) str_c("llm_asy_", x))
    
    # Exact LLM
    v <- as_tibble(t) |> 
        arrange(ability, group, response) |> 
        pull(n)
    pr <- computeUProbsCpp(matrix(v))
    sample <- algstat::metropolis(v, moves$llm_no23, 
                                  iter = iter, burn = burn, thin = thin)$steps
    llm_exc_no23_p <- mean(computeUProbsCpp(sample) <= pr)
    sample <- algstat::metropolis(v, moves$llm_no3w, 
                                  iter = iter, burn = burn, thin = thin)$steps
    llm_exc_no3w_p <- mean(computeUProbsCpp(sample) <= pr)
    llm_exc <- dif_tests(no23_abort = FALSE,
                         no23_p = llm_exc_no23_p,
                         no3w_abort = FALSE,
                         no3w_p = llm_exc_no3w_p,
                         alpha = alpha) |> 
        select(-comp) |> 
        rename_with(\(x) str_c("llm_exc_", x))
    
    # Information to be recorded
    llm <- tibble(
        has_no_zeros = has_no_zeros(t),
        llm_asy_no23_mle = llm_asy_no23$mle_exists,
        llm_asy_no23_heuristic = llm_asy_no23$heuristic,
        llm_asy_no3w_mle = llm_asy_no3w$mle_exists,
        llm_asy_no3w_heuristic = llm_asy_no3w$heuristic,
        llm_asy_no23_p = llm_asy_no23$p,
        llm_asy_no3w_p = llm_asy_no3w$p,
        llm_asy_comp_p = llm_asy_comp_p,
        llm_exc_no23_p = llm_exc_no23_p,
        llm_exc_no3w_p = llm_exc_no3w_p
        ) |> 
        bind_cols(llm_asy, llm_exc)
    
    # Logistic regression tests
    
    if (dim(t)[3] != 2) {
        # Record NA for LRM stuff if item is not dichotomous
        lrm <- tibble(
            lrm_asy_no23_mle = NA,
            lrm_asy_no23_heuristic = NA,
            lrm_asy_no3w_mle = NA,
            lrm_asy_no3w_heuristic = NA,
            lrm_asy_no23_p = NA,
            lrm_asy_no3w_p = NA,
            lrm_asy_comp_p = NA,
            lrm_exc_no23_p = NA,
            lrm_exc_no3w_p = NA,
            lrm_asy_main = NA,
            lrm_asy_swap = NA,
            lrm_asy_comp = NA,
            lrm_exc_main = NA,
            lrm_exc_swap = NA
        )
    } else {
        # Asymptotic LRM
        lrm_asy_no23 <- mle(t, no23, "lrm")
        lrm_asy_no3w <- mle(t, no3w, "lrm")
        lrm_asy_comp_p <- pchisq(lrm_asy_no23$lrt - lrm_asy_no3w$lrt, 
                                 lrm_asy_no23$df - lrm_asy_no3w$df, 
                                 lower.tail = FALSE)
        
        lrm_asy <- dif_tests(no23_abort = !lrm_asy_no23$mle_exists,
                             no23_p = lrm_asy_no23$p,
                             no3w_abort = !lrm_asy_no3w$mle_exists,
                             no3w_p = lrm_asy_no3w$p,
                             comp_p = lrm_asy_comp_p,
                             alpha = alpha) |> 
            rename_with(\(x) str_c("lrm_asy_", x))
        
        # Exact LRM
        v <- as_tibble(t) |> 
            arrange(desc(response), ability, group) |> 
            pull(n)
        pr <- computeUProbsCpp(matrix(v))
        sample <- algstat::metropolis(v, moves$lrm_no23, 
                                      iter = iter, burn = burn, thin = thin)$steps
        lrm_exc_no23_p <- mean(computeUProbsCpp(sample) <= pr)
        sample <- algstat::metropolis(v, moves$lrm_no3w, 
                                      iter = iter, burn = burn, thin = thin)$steps
        lrm_exc_no3w_p <- mean(computeUProbsCpp(sample) <= pr)
        lrm_exc <- dif_tests(no23_abort = FALSE,
                             no23_p = lrm_exc_no23_p,
                             no3w_abort = FALSE,
                             no3w_p = lrm_exc_no3w_p,
                             alpha = alpha) |> 
            select(-comp) |> 
            rename_with(\(x) str_c("lrm_exc_", x))
        
        # Information to be recorded
        lrm <- tibble(
            lrm_asy_no23_mle = lrm_asy_no23$mle_exists,
            lrm_asy_no23_heuristic = lrm_asy_no23$heuristic,
            lrm_asy_no3w_mle = lrm_asy_no3w$mle_exists,
            lrm_asy_no3w_heuristic = lrm_asy_no3w$heuristic,
            lrm_asy_no23_p = lrm_asy_no23$p,
            lrm_asy_no3w_p = lrm_asy_no3w$p,
            lrm_asy_comp_p = lrm_asy_comp_p,
            lrm_exc_no23_p = lrm_exc_no23_p,
            lrm_exc_no3w_p = lrm_exc_no3w_p
            ) |> 
            bind_cols(lrm_asy, lrm_exc)
    }
    
    # Return everything recorded
    bind_cols(llm, lrm)
}