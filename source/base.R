# This file defines a number of useful functions that are used across the R 
# script files in this repository. The intention was to make them convenient
# for use for the type of DIF analyses we wanted to conduct, not necessarily
# in complete generality. 

library(tidyverse)
library(detectseparation)
library(algstat)

# ==============================================================================
# Miscellaneous
# ==============================================================================

FT <- c(FALSE, TRUE)

# Pull an algstat C++ functions into the namespace
computeUProbsCpp <- function(x) {
    .Call("_algstat_computeUProbsCpp", PACKAGE = "algstat", x)
}

# ==============================================================================
# Model and Family Selection
# ==============================================================================

# Family specifications
llm <- "llm"
lrm <- "lrm"

# Facet specifications of models
no23 <- list(c(1, 2), c(1, 3))
no3w <- list(c(1, 2), c(1, 3), c(2, 3))
full <- list(c(1, 2, 3))

# Some helper variables
all_families <- list(llm = llm, lrm = lrm)
all_models <- list(no23 = no23, no3w = no3w, full = full)

# Helper function to help create togglers, which are named lists of booleans
# to help with the casework that happens in many of the functions below. The
# input "input" can be a toggler, a string, or a facet specification. The 
# second input "all" is intended to be one of the lists all_families or
# all_models that appears above.
toggler_helper <- function(input, all) {
    if (is.list(input) && 
        all(names(all) %in% names(input)) &&
        is.logical(unlist(input))) {
        # If the input is a toggler, select the booleans that correspond to
        # names in the named list all
        out <- map2(all, names(all), \(x, y) input[[y]])
    } else if (is.character(input) && length(input) == 1) {
        # If the input is a single string, check for presence of one of the 
        # names in the named list all
        out <- str_detect(input, names(all)) |> 
            as.list() |> 
            set_names(names(all))
    } else {
        # Otherwise, check for set equality
        out <- map(all, \(x) setequal(input, x))
    }
    
    # Error check
    if (reduce(out, sum) != 1) {
        stop("Family or model selection error")
    }
    
    # Return
    out
}

# Helper function to create a toggle for families
family_toggler <- function(family) {
    toggler_helper(family, all_families)
}

# Helper function to create a list of toggles for models
model_toggler <- function(model) {
    toggler_helper(model, all_models)
}

# "Main" helper function to create a list of toggles for models and families.
# It uses the first argument to select the model using the model_toggler()
# above. If the second argument is NULL, it also uses the first argument to 
# select the family using family_toggler(). Otherwise, it uses the second 
# argument to choose the family. 
toggler <- function(model, family = NULL) {
    x <- model_toggler(model)
    if (is_null(family)) {
        y <- family_toggler(model)
    } else {
        y <- family_toggler(family)
    }
    append(x, y)
}

# Helper function to extract model facet specification from a valid toggler,
# but it does not check that the toggler is valid
facets_from_toggler <- function(toggle) {
    all_models[[keep(names(all_models), \(x) toggle[[x]])]]
}

# ==============================================================================
# Configuration Matrices and Markov Moves
# ==============================================================================

# Turns a table into a vector with the property that the matrix computed by
# the configuration_matrix function below, times the output fo this function,
# gives the sufficient statistics of the table
vectorize <- function(t, family) {
    toggle <- family_toggler(family)
    
    if (toggle$llm) {
        # For log-linear models, use lex order
        out <- as_tibble(t) |> 
            arrange(ability, group, response) |> 
            pull(n)
    } else {
        # For logistic regressions, use a different order
        out <- as_tibble(t) |> 
            arrange(desc(response), ability, group) |> 
            pull(n)
    }
    
    # Name and return
    set_names(out, vectorize_names(dimnames(t), toggle))
}

# Helper function to name the output of the vectorize function above. It also
# gets used to name the columns of configuration matrices, and also the rows of
# matrices of Markov moves. 
vectorize_names <- function(dimnames, family) {
    # Currently, if any dimension name has more than 1 character, return NULL.
    # This is not a problem for the code as it stands, but probably one should
    # find a workaround eventually.
    if (any(nchar(unlist(dimnames)) != 1)) {
        return(NULL)
    }
    
    # Make toggle
    toggle <- family_toggler(family)
    
    # Prepare names
    grid <- expand.grid(dimnames) |> 
        as_tibble() |> 
        mutate(label = str_c(ability, group, response))
    
    # Make names
    if (toggle$llm) {
        out <- grid |> 
            arrange(ability, group, response) |> 
            pull(label)
    } else {
        out <- grid |> 
            arrange(desc(response), ability, group) |> 
            pull(label)
    }
    
    # Return
    out
}

# Undoes the operation of the vectorize function above
unvectorize <- function(v, dimnames, family) {
    # Make toggle
    toggle <- family_toggler(family)
    
    # Make tibble of dimnames
    grid <- expand.grid(dimnames) |>
        as_tibble()
    
    # Insert counts into tibble and sort lexicographically
    if (toggle$llm) {
        grid <- grid |> 
            arrange(ability, group, response) |> 
            mutate(n = v)
    } else {
        grid <- grid |> 
            arrange(desc(response), ability, group) |> 
            mutate(n = v) |> 
            arrange(ability, group, response)
    }
    
    # Make table
    t <- vec2tab(grid$n, map_dbl(dimnames, length)) |> 
        as.table()
    dimnames(t) <- dimnames
    
    # Return
    t
}

# Helper variable. Used in the hmat_renamer function below. Matches a list
# that is defined in algstat (namely, in hmat.R). 
hmat_glyphs <- c(1:9, letters, LETTERS)

# Helper function to rename the row names output by hmat using the names from
# dimnames instead
hmat_renamer <- function(hmat_name, dimnames) {
    # Currently, if any dimension name has more than 1 character, return NULL.
    # This is not a problem for the code as it stands, but probably one should
    # find a workaround eventually.
    if (any(nchar(unlist(dimnames)) != 1)) {
        return(NULL)
    }
    
    # Split hmat_name into individual characters
    x <- str_split(hmat_name, "")[[1]]
    
    # Error checking
    if (length(x) != length(dimnames)) {
        stop("hmat_name incompatible with dimnames")
    }
    
    # Helper function that converts the ith character in the hmat string
    lambda <- function(i) {
        # Return + as is
        if (x[i] == "+") {
            return("+")
        }
        
        # Match
        m <- match(x[i], hmat_glyphs)
        
        # Error checking
        if (is.na(m)) {
            stop("hmat_name invalid")
        } else if (m > length(dimnames[[i]])) {
            stop("hmat_name incompatible with dimnames")
        }
        
        # Return
        as.character(dimnames[[i]][m])
    }
    
    # Apply, flatten, return
    map_chr(seq_along(x), lambda) |> 
        str_flatten()
}

# Compute configuration matrix for the given model. Names the rows and columns
# of the matrix using the names from dimnames. 
# Note: Returns NULl if requesting a configuration matrix for logistic 
# regression for non-dichotomous model. 
configuration_matrix <- function(dimnames, model, family = NULL) {
    # Make toggler
    toggle <- toggler(model, family)
    
    # Get dimensions
    dimnames <- map(dimnames, as.numeric)
    d <- map_dbl(dimnames, length)
    
    if (toggle$llm) {
        # Log-linear models
        # Use hmat to make configuration matrices
        out <- hmat(d, facets_from_toggler(toggle))
        
        # Make row names
        rownames <- dimnames(out)[[1]] |> 
            map_chr(\(x) hmat_renamer(x, dimnames))
        
        # Name output
        dimnames(out) <- list(rownames, vectorize_names(dimnames, "llm"))
        
    } else if (d[3] == 2) {
        # Logistic regression
        # Make auxiliary matrix
        aux <- expand_grid(intercept = 1, 
                           ability = dimnames[[1]], 
                           group = dimnames[[2]]) |> 
            mutate(interaction = ability * group) |> 
            t()
        
        # Use lawrence to make configuration matrices
        if (toggle$no23) {
            out <- lawrence(aux[1:2, ])
        } else if (toggle$no3w) {
            out <- lawrence(aux[1:3, ])
        } else {
            out <- lawrence(aux)
        }
        
        # Make row names
        rownames <- c(str_c("++", dimnames$response[2]), "sig")
        if (toggle$no3w || toggle$full) {
            rownames <- str_c("+", dimnames$group[2], dimnames$response[2]) |> 
                append(rownames, values = _)
        }
        if (toggle$full) {
            rownames <- append(rownames, "tau")
        }
        rownames <- expand_grid(ability = dimnames$ability, 
                                group = dimnames$group) |> 
            mutate(label = str_c(ability, group, "+")) |> 
            pull(label) |> 
            append(rownames, values = _)
        
        # Name output
        dimnames(out) <- list(rownames, vectorize_names(dimnames, "lrm"))
        
    } else {
        # Logistic regression in the nondichotomous setting
        out <- NULL
    }
    
    # Return
    out
}

# Compute sufficient statistics of the given table under the given model
sufficient_statistics <- function(t, model, family) {
    A <- configuration_matrix(dimnames(t), model, family)
    v <- vectorize(t, family)
    (A %*% v)[, 1]
}

# Compute configuration matrices for all relevant models. The input should be
# a list of numeric vectors, where the first vector is the list of ability
# levels, the second is the list of groups levels, and the third is the list
# of response levels. The output is a named list of matrices, where the names
# are llm_no23, llm_no3w, llm_full, lrm_no23, lrm_no3w, lrm_full.
configuration_matrices <- function(dimnames) {
    # Make list of labels
    labels <- expand_grid(model = names(all_models),
                        family = names(all_families)) |> 
        mutate(label = str_c(family, "_", model)) |> 
        pull(label)
    
    # Compute and return
    map(labels, \(x) configuration_matrix(dimnames, x)) |> 
        set_names(labels)
}

# Compute Markov moves for relevant models. The input should be a list of 
# numeric vectors, where the first vector is the list of ability levels, the 
# second is the list of groups levels, and the third is the list of response 
# levels. The output is a named list of matrices whose columns are the moves, 
# where the names are llm_no23, llm_no3w, llm_full, lrm_no23, lrm_no3w, and 
# lrm_full. Requires 4ti2. 
markov_moves <- function(dimnames) {
    # Helper function to compute Markov moves
    lambda <- function(A) {
        if (is_null(A)) {
            return(NULL)
        }
        markov(A, p = "arb")
    }
    
    # Note that, for saturated models (like llm_full, and also lrm_full in
    # certain cases), the Markov basis is empty, and markov() returns a matrix
    # with a single column of 0s. It also throws a warning in this case, and 
    # this is suppressed.
    out <- map(configuration_matrices(dimnames), lambda) |> 
        suppressWarnings()
    
    # Helper function for attaching names to rows
    lambda <- function(x, y) {
        if (is_null(x)) {
            return(NULL)
        }
        family <- ifelse(str_detect(y, "llm"), "llm", "lrm")
        dimnames(x) <- list(vectorize_names(dimnames, family))
        x
    }
    
    # Attach names and return
    map2(out, names(out), lambda)
}

# Count tables in the fiber of the given table. Requires LattE.
fiber_size <- function(t, model, family = NULL) {
    count_tables(t, configuration_matrix(dimnames(t), model, family))
}

# ==============================================================================
# MLE-Related Functions
# ==============================================================================

# Helper function that reorganize table into a tibble that's useful for
# logistic regression For each ability-group pair, it records the proportion
# of correct responses in a column called p, and the total number of subjects
# in a column called n. The table must be dichotomous for this to make sense, 
# but this condition is not checked.
reorganize <- function(t) {
    as_tibble(t) |> 
        mutate(across(everything(), as.numeric)) |> 
        group_by(ability, group) |> 
        summarize(p = sum(keep(n, response == 1)) / sum(n), 
                  n = sum(n), 
                  .groups = "drop")
}

# Helper function that reverse the process of the reorganize function above.
# Caution: unreorganize(reorganize(t)) may differ from t by small amounts
# due to floating point arithmetic. This is a useful for computing
# expected counts under logistic regression MLEs (or rather, for putting those
# expected counts into a table that matches the original table).
unreorganize <- function(u) {
    dimnames <- list(ability = unique(u$ability), 
                     group = unique(u$group), 
                     response = 0:1)
    t <- rbind((1 - u$p) * u$n, u$p * u$n) |> 
        vec2tab(map_dbl(dimnames, length)) |> 
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
    to_add <- c(min(dimnames$ability) - 1, max(dimnames$ability) + 1)
    expand_grid(ability = to_add, 
                group = dimnames$group, 
                response = dimnames$response, 
                n = 0) |> 
        bind_rows(w) |> 
        arrange(ability, group, response)
}

# Helper function to check for "cross zeros." In other words, consructs a line
# ell passing through the two points (a0, 0) and (a1, 1) and then checks to see
# if every point (a, g) strictly to the left got the question right and every
# point strictly to the right got the question wrong, or vice versa. 
# This condition appears in the concrete characterization of separation for
# the no23 and no3w logistic regression models. This input w should be the
# result of tibble_and_append above, but this condition is not checked. This
# means that it has 4 columns (ability, group, response, and n for counts)
# and its been extended with an extra ability level on either side with
# no observations. 
cross_zeros <- function(w, a0, a1) {
    # The column q has the property that q == 0 iff (a, g) lies on the line
    # ell passing through (a0, 0) and (a1, 1). We then produce a vector x
    # of length 4 such that x[1], x[2] are the counts of 0, 1 responses 
    # strictly on one side of the line, and x[3], x[4] are the corresponding 
    # counts strictly on the other side of the line. Note that the summarize 
    # operation arranges the table in lex order (side, response) and this 
    # arrangement is preserved after groups are dropped.
    x <- w |> 
        mutate(q = group * (a1 - a0) - (ability - a0)) |> 
        filter(q != 0) |> 
        mutate(side = (q > 0)) |> 
        group_by(side, response) |> 
        summarize(n = sum(n), .groups = "drop") |> 
        pull(n)
    
    # This should never happen...
    if (length(x) != 4) {
        stop("Error")
    }
    
    # Return
    (x[1] == 0 && x[4] == 0) || (x[2] == 0 && x[3] == 0)
}

# Helper function for checking if a 2x2x2 has "diagonal zeros," ie, if it has
# zeroes in opposite corners of the cube. There are 4 ways of having zeroes
# in opposite corners, and it just checks for any one of these conditions. 
# Used for MLE existence check for log-linear no3w model.
# The input must be a 2x2x2 table, but this condition is not checked. 
diagonal_zeros <- function(t) {
    (t[1, 1, 1] == 0 && t[2, 2, 2] == 0) ||
        (t[1, 1, 2] == 0 && t[2, 2, 1] == 0) ||
        (t[1, 2, 1] == 0 && t[2, 1, 2] == 0) ||
        (t[1, 2, 2] == 0 && t[2, 1, 1] == 0)
}

# Helper function for collapsing a IxJxK table into a 2x2x2 table. The input 
# collapsing should be a list of three boolean vectors of lengths I, J, K,
# indicating how the levels should be grouped. In order go guarantee that the
# output is indeed 2x2x2, none of the three boolean vectors can be all TRUE or 
# all FALSE, but this condition is not checked. 
collapse <- function(t, collapsing) {
    # Helper function: Takes as input three boolean values x, y, z, and outputs
    # the sum of the corresponding entries in the table, where "corresponding
    # entries" refers to the entries t[i, j, k] where collapsing[[1]][i] == x, 
    # collapsing[[2]][j] == y, and collapsing[[3]][k] == z. 
    lambda <- function(x, y, z) {
        sum(t[which(collapsing[[1]] == x), 
              which(collapsing[[2]] == y), 
              which(collapsing[[3]] == z)])
    }
    
    # Make 2x2x2 table. Note that expand_grid constructs the tibble pre-sorted
    # in lexicographic order. 
    out <- expand_grid(a = FT, b = FT, c = FT) |> 
        mutate(n = pmap_dbl(list(a, b, c), lambda)) |> 
        pull(n) |> 
        vec2tab(c(2, 2, 2))
    
    # Attach readable names for checking
    if (!is_null(dimnames(t))) {
        # Helper function: Takes as input a string vector x from dimnames(t)
        # and a boolean vector y from collapsing, and returns a vector of two 
        # strings: one concatenating the names from x which correspond to 
        # FALSE in y, and the other concatenating the names from x which
        # correspond to TRUE in y. 
        lambda <- function(x, y) {
            map_chr(FT, \(z) str_flatten(keep(x, y == z)))
        }
        dimnames(out) <- map2(dimnames(t), collapsing, lambda)
    }
    
    # Return
    out
}

# Helper function to generate a list of all nontrivial binary partitions
# of a set of of size n. A binary partition is specified as a list of n
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
# "diagonal zeroes" (cf. diagonal_zeros above). Used for MLE existence check 
# for llm no3w if at least one of I, J, K is 2 (but this is not checked).
# It would be possible to make this function sligtly more efficient by short
# circuiting the calculation: stop and return TRUE as soon as one collapsing
# with diagonal zeros is found. We have not done this. 
collapsing_has_diagonal_zeros <- function(t) {
    # Get partitions of each dimension
    partitions <- map(unname(dim(t)), binary_partitions)
    
    # Helper function to construct collapsings. Takes as input three integers
    # x, y, z, and outputs the collapsing obtained by using the collapsing
    # specified below. 
    lambda <- function(x, y, z) {
        collapsing <- list(partitions[[1]][x, ], 
                           partitions[[2]][y, ], 
                           partitions[[3]][z, ])
        collapse(t, collapsing)
    }
    
    # Make tibble of partitions, collapsings, and checks for diagonal zeros
    df <- map(partitions, \(x) seq((dim(x))[1])) |> 
        expand.grid() |> 
        as_tibble() |> 
        mutate(collapsing = pmap(list(Var1, Var2, Var3), lambda)) |> 
        mutate(diagonal_zeros = map_lgl(collapsing, diagonal_zeros))
    
    # Return
    any(df$diagonal_zeros)
}

# Check if table t has MLE for the model (no23 or no3w or full) and family (llm
# or lrm) specified by the toggle. This toggle must be a valid toggler returned
# by the toggler function, but this condition is not checked. 
# The function attempts to short circuit calculations whenever possible, 
# especially for logistic regression since the implementation of 
# detectseparation is computationally expensive. To skip that short circuiting, 
# pass avoid_lp = FALSE, but this option mostly exists for testing purposes.
mle_exists <- function(t, toggle, avoid_lp = TRUE) {
    # No sampling zeros; sufficient, but not necessary, for any MLE existence.
    # It is also necessary for MLE existence for full log-linear model, and
    # the full logistic regression model if the table is 2x2x2, because those
    # models are saturated.
    out <- all(t > 0)
    is_saturated <- toggle$full &&
        (toggle$llm || (toggle$lrm && all(dim(t) == 2)))
    if (out || is_saturated) {
        return(out)
    }
    
    if (toggle$llm) {
        # Log-linear models
        # First, check that line sums n_{a,g,+} and n_{a,+,r} are nonzero.
        # Necessary and sufficient for MLE existence for no23.
        # Sufficient, but not necessary for MLE existence for no3w.
        out <- all(marginSums(t, c(1, 2)) > 0, marginSums(t, c(1, 3)) > 0)
        if (!out || toggle$no23) {
            return(out)
        }
        
        # For the no3w model, first check that all line sums are nonzero
        # This is necessary, but not sufficient, for MLE existence
        out <- all(marginSums(t, c(2, 3)) > 0)
        if (!out) {
            return(out)
        }
        
        # Check for non-existence of diagonal pattern via 2x2x2 collapsings
        # à la Eriksson et al 2006. This check is somewhat expensive
        # (exponential in the dimensions of the table) so it's short circuited
        # by the above. 
        if (any(dim(t) == 2)) {
            return(!collapsing_has_diagonal_zeros(t))
        }
        
        # Otherwise, we are in a table where no dimension is 2, in which case
        # a characterization of MLE existence is unimplemented
        stop("LLM no3w MLE existence unimplemented when no dimension is 2")
    } else {
        # Logistic regression models
        # First, check that we are in a dichotomous situation
        if (dim(t)[3] != 2) {
            stop("Response must be dichotomous for logistic regressions")
        }
        
        # Check that line sums n_{a,g,+} are nonzero.
        # Necessary, but not sufficient, for MLE existence.
        if (any(marginSums(t, c(1, 2)) == 0)) {
            return(FALSE)
        }
        
        # Check some easier cases of separation.
        # This short circuits the expensive method via linear programming
        # implemented in the detectseparation package.
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
                # First, check that line sums n_{+,g,r} are nonzero.
                # This is redundant, but maybe it is useful to short circuit
                # the computations below...?
                if (any(marginSums(t, c(2, 3)) == 0)) {
                    return(FALSE)
                }
                
                # Check for separation along any pair of (distinct) ability 
                # levels; the case of identical ability levels has already
                # been covered by the above. 
                # This is enough to check for no3w, and if it fails, the MLE 
                # for full cannot exist either.
                # It is possible to make this more efficient by short 
                # circuiting in the sense that we return as soon as a pair
                # with cross zeros is found, but we have not done this. 
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
        
        # If the above short circuit failed, check via linear programming
        # à la Konis 2007. This computation is expensive for large sample
        # sizes, but if the above went through, it should only be called when
        # there are more than 2 groups (which we do not implement) or if
        # we're interested in the full lrm model and the no3w MLE exists
        # (and we aren't currently using anyting about the full lrm model).
        
        # Convert table t to a tidy tibble
        v <- as_tibble(t) |> 
            mutate(across(everything(), as.numeric)) |> 
            uncount(n)
        
        # Run detect separation
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
        
        # Return
        return(out)
    }
}

# Check if table has any sampling zeros
has_no_zeros <- function(t) {
    mle_exists(t, toggler(full, llm))
}

# Check if table satisfies the heuristic of having sufficiently many expected
# counts large enough, where "sufficiently many" means at least
# min_percent of the total number (default: 80%) and "large enough" means
# greater than or equal to min_number (default: 5)
heuristic_check <- function(expected, min_percent = 0.8, min_number = 5) {
    percent_large <- sum(expected >= min_number) / length(expected)
    percent_large >= min_percent
}

# Main function: collect all relevant information related to the MLE of the
# input table for the input model (as specified by the model and family 
# arguments, which get passed to the toggler() function to choose a model
# and family; see there for format requirements). 
# The input threshold is used to determine the type of MLE existence checking
# If it is <= 0, the exact test of mle_exists is used. If it is > 0, an
# approximate test using the fit returned by loglin or glm is used to determine
# if the fitted distribution lies in the interior of the probability simplex,
# with all entries being at least the threshold. This option mostly exists
# for error checking. 
# The output collects the following information:
# * mle_exists: TRUE if the MLE exists
# * heuristic: TRUE if expected counts pass the heuristic check
# * expected: three-way table of expected counts
# * lrt: G test statistic
# * df: degrees of freedom
# * p: tail area part lrt in chi-square distribution with degrees of freedom df
# If the MLE does not exist, the mle_exists and heuristic values are FALSE,
# and all of the remaining values are NA. 
mle <- function(t, model, family = NULL, threshold = 0) {
    # Make toggle and error check
    toggle <- toggler(model, family)
    
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
            # Log-linear models; use loglin
            x <- loglin(t, facets_from_toggler(toggle), 
                        fit = TRUE, 
                        print = FALSE)
            expected <- x$fit
            lrt <- x$lrt
            df <- x$df
        } else {
            # Logistic regression; check that we are dichotomous
            if (dim(t)[3] != 2) {
                stop("Response must be dichotomous for logistic regressions")
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
        mle_exists <- all(marginSums(t, c(1, 2)) > 0, mle_dist > threshold)
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

# ==============================================================================
# Construct Simulation Models
# ==============================================================================

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
    if (length(base$uv) != (2 * K) || length(dif) != (2 * K)) {
        stop("Mismatch of parameter lengths")
    }
    u <- base$uv[1:K]
    v <- base$uv[(K + 1):(2 * K)]
    u_delta <- dif[1:K]
    v_delta <- dif[(K + 1):(2 * K)]
    
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
    r_given_ag <- sweep(alpha, c(1, 2), marginSums(alpha, c(1, 2)), `/`)
    
    # Array of joint distribution of (A, G)
    ag <- sweep(base$ability_dist, 2, base$group_dist, `*`)
    
    # Array of joint distribution of (A, G, R)
    sweep(r_given_ag, c(1, 2), ag, `*`)
}

# Helper function that takes as input a joint distribution of three variables
# specified as a three-way array, and returns the conditional distribution of
# the third variable given the first two again as a three-way array.
conditional <- function(dist) {
    sweep(dist, c(1, 2), marginSums(dist, c(1, 2)), `/`)
}

# Helper function to compute the Rényi alpha-divergence from q to p, where p
# and q are numerical vectors of equal length containing probabilities that
# sum to 1. Defaults to alpha = 1, which corresponds to Kullback-Leibler
# divergence.
renyi <- function(p, q, alpha = 1) {
    if (alpha < 0) {
        stop("Alpha parameter must be nonnegative")
    } else if (length(p) != length(q) ||
               any(p < 0) ||
               any(q < 0) ||
               abs(sum(p) - 1) > .Machine$double.eps ||
               abs(sum(q) - 1) > .Machine$double.eps) {
        stop("Input vectors must be comparable probability vectors")
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
        out <- 1 / (alpha - 1) * out
    }
    
    # Return
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
        map_dbl(\(a) kl(r_given_ag[a, 2, ], r_given_ag[a, 1, ]))
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
    c <- lpnorm(x, p)
    if (c == 0) {
        out <- x
    } else {
        out <- x / c
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
# It then draws an integer from 2:max_dim, with the probability of 2 being
# dichotomous_prob and the others being equally likely, and uses the integers
# from 0 up to one minus the integer drawn as response levels. If max_dim is
# 2, dichotomous_prob is ignored (ie, it is automatically equal to 1). 
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
                                       models_file = NULL, 
                                       beta = 10, 
                                       max_dim = 7, 
                                       dichotomous_prob = 0.5, 
                                       norms = c(1, 2), 
                                       include_moves = TRUE) {
    # Enforce max_dim is integer at least 2
    max_dim <- as.integer(max_dim)
    if (max_dim < 2) {
        stop("Input max_dim must be an integer greater than or equal to 2")
    }
    
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
        (\(p) c(p, 1 - p))()
    names(group_dist) <- group_levels
    
    # Generate ability levels
    I <- sample(2:max_dim, 1)
    ability_levels <- 0:(I - 1)
    
    # Generate distributions of ability conditioned on group
    # Each column is the ability distribution of a fixed group
    ability_dist <- expand_grid(p = rbeta(2, beta, beta), ability_levels) |> 
        mutate(q = dbinom(ability_levels, I - 1, p)) |> 
        pull(q) |> 
        matrix(I, J, 
               dimnames = list(ability = ability_levels, group = group_levels))
    
    # Generate number of response levels
    if (max_dim == 2) {
        # Nothing to do in this case
        K <- 2
    } else {
        # Generate number of response levels, with prob of 2 being given
        num_polytomous_options <- max_dim - 2
        prob <- c(dichotomous_prob, 
                  rep((1 - dichotomous_prob) / num_polytomous_options, 
                      num_polytomous_options))
        K <- sample(2:max_dim, 1, prob = prob)
    }
    
    # Generate response levels
    response_levels <- 0:(K - 1)
    
    # Dimension vector
    dimnames <- list(ability = ability_levels, 
                     group = group_levels, 
                     response = response_levels)
    
    # Response distribution for group 0
    u <- runif(K - 1, -1, 1) |> 
        to_constants_perp()
    v <- runif(K - 1, -1, 1) |> 
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
    v_delta <- runif(K - 1, -1, 1) |> 
        to_constants_perp()
    unif <- c(u_delta, v_delta) |> 
        normalize()
    
    # Non-uniform DIF vector
    u_delta <- runif(K - 1, -1, 1) |> 
        to_constants_perp()
    v_delta <- runif(K - 1, -1, 1) |> 
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

# ==============================================================================
# Simulate Data
# ==============================================================================

# Generate a three-way table of observation counts following the three-way
# input distribution, where the number of observations is the given
# sample_size. Insist on the table passing the given predicate, but try at
# most max times.
simulate_table <- function(dist, 
                           sample_size, 
                           predicate = \(t) TRUE, 
                           max = 1e3) {
    count <- 0
    repeat {
        # Generate table
        t <- rmultinom(1, sample_size, tab2vec(dist)) |> 
            as.vector() |> 
            vec2tab(dim(dist)) |> 
            as.table()
        dimnames(t) <- dimnames(dist)
        
        # Check whether we keep generating tables
        count <- count + 1
        if (predicate(t)) {
            break
        } else if (count > max) {
            stop("Unable to pass predicate within max iterations")
        }
    }
    
    # Return
    t
}

# ==============================================================================
# Task Function
# ==============================================================================

# Strings to indicate DIF test results:
# * fail = test could not determine if there is DIF
# * none = test concluded there is no DIF
# * unclassifiable = test concluded there is no DIF but could not classify
# * uniform = test concluded there is uniform DIF
# * nonuniform = test concluded there is nonuniform DIF
dif_test_results <- list(fail = "fail", 
                         none = "none", 
                         unclassifiable = "unclass", 
                         uniform = "unif", 
                         nonuniform = "nonunif")

# Conduct a battery of DIF tests: the main one, the swapped variant, and the
# variant with comparison (if comp_p is provided).
dif_tests <- function(no23_abort, 
                      no23_p, 
                      no3w_abort, 
                      no3w_p, 
                      comp_p = NULL, 
                      alpha = 0.05) {
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
                       alpha = 0.05) {
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
    v <- vectorize(t, "llm")
    pr <- computeUProbsCpp(matrix(v))
    sample <- metropolis(v, moves$llm_no23, 
                                  iter = iter, 
                                  burn = burn, 
                                  thin = thin)$steps |> 
        suppressMessages()
    llm_exc_no23_p <- mean(computeUProbsCpp(sample) <= pr)
    sample <- metropolis(v, moves$llm_no3w, 
                                  iter = iter, 
                                  burn = burn, 
                                  thin = thin)$steps |> 
        suppressMessages()
    llm_exc_no3w_p <- mean(computeUProbsCpp(sample) <= pr)
    llm_exc <- dif_tests(no23_abort = FALSE, 
                         no23_p = llm_exc_no23_p, 
                         no3w_abort = FALSE, 
                         no3w_p = llm_exc_no3w_p, 
                         alpha = alpha) |> 
        select(-comp) |> 
        rename_with(\(x) str_c("llm_exc_", x))
    
    # Information to be recorded
    llm <- tibble(has_no_zeros = has_no_zeros(t), 
                  llm_asy_no23_mle = llm_asy_no23$mle_exists, 
                  llm_asy_no23_heuristic = llm_asy_no23$heuristic, 
                  llm_asy_no3w_mle = llm_asy_no3w$mle_exists, 
                  llm_asy_no3w_heuristic = llm_asy_no3w$heuristic, 
                  llm_asy_no23_p = llm_asy_no23$p, 
                  llm_asy_no3w_p = llm_asy_no3w$p, 
                  llm_asy_comp_p = llm_asy_comp_p, 
                  llm_exc_no23_p = llm_exc_no23_p, 
                  llm_exc_no3w_p = llm_exc_no3w_p) |> 
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
        v <- vectorize(t, "lrm")
        pr <- computeUProbsCpp(matrix(v))
        sample <- metropolis(v, moves$lrm_no23, 
                                      iter = iter, 
                                      burn = burn, 
                                      thin = thin)$steps |> 
            suppressMessages()
        lrm_exc_no23_p <- mean(computeUProbsCpp(sample) <= pr)
        sample <- metropolis(v, moves$lrm_no3w, 
                                      iter = iter, 
                                      burn = burn, 
                                      thin = thin)$steps |> 
            suppressMessages()
        lrm_exc_no3w_p <- mean(computeUProbsCpp(sample) <= pr)
        lrm_exc <- dif_tests(no23_abort = FALSE, 
                             no23_p = lrm_exc_no23_p, 
                             no3w_abort = FALSE, 
                             no3w_p = lrm_exc_no3w_p, 
                             alpha = alpha) |> 
            select(-comp) |> 
            rename_with(\(x) str_c("lrm_exc_", x))
        
        # Information to be recorded
        lrm <- tibble(lrm_asy_no23_mle = lrm_asy_no23$mle_exists, 
                      lrm_asy_no23_heuristic = lrm_asy_no23$heuristic, 
                      lrm_asy_no3w_mle = lrm_asy_no3w$mle_exists, 
                      lrm_asy_no3w_heuristic = lrm_asy_no3w$heuristic, 
                      lrm_asy_no23_p = lrm_asy_no23$p, 
                      lrm_asy_no3w_p = lrm_asy_no3w$p, 
                      lrm_asy_comp_p = lrm_asy_comp_p, 
                      lrm_exc_no23_p = lrm_exc_no23_p, 
                      lrm_exc_no3w_p = lrm_exc_no3w_p) |> 
            bind_cols(lrm_asy, lrm_exc)
    }
    
    # Return everything recorded
    bind_cols(llm, lrm)
}
