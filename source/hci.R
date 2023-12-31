# This script is used to conduct the analysis of the HCI dataset that appears
# in the paper.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(xtable)
library(ShinyItemAnalysis)

# Pull an algstat C++ functions into the namespace
computeUProbsCpp <- \(x) .Call("_algstat_computeUProbsCpp", PACKAGE = "algstat", x)

# ==============================================================================
# Initialize
# ==============================================================================

# Load data
data(HCI)

# Choose grouping variable: either "major" or "gender"
grouper <- "major"

# Discretizes a numerical vector x into discrete levels 0, 1, ..., k-1
discretize <- function(x, k) {
    as.integer(cut_interval(x, n = k)) - 1
}

# ==============================================================================
# Focus Item
# ==============================================================================

# This code is intended to match the snippets that appear in our paper closely, 
# but not exactly. It is usually not the quicket code or the code that takes
# the most advantage of other functions that we've defined in this repository,
# but the point is for readers of the paper to be able to follow along with the
# code snippets in the paper using this section, if they'd like. 

# Choose focus item: Item 1, Item 2, ..., Item 20
focus_item <- 17

# Construct three-way table
t <- tibble(ability = discretize(HCI$total, 6), 
            group = HCI[[grouper]], 
            response = HCI[[paste("Item", focus_item)]]) |> 
    table()
t

# Asymptotic log-linear
no23 <- list(c(1, 2), c(1, 3))                                  # Facet specification of model
x <- loglin(t, no23, fit = TRUE)                                # Fit the model
sum(x$fit >= 5) / length(x$fit)                                 # Percent large expected counts
pchisq(x$lrt, x$df, lower.tail = FALSE)                         # p-value

no3w <- list(c(1, 2), c(1, 3), c(2, 3))                         # Facet specification of model
x <- loglin(t, no3w, fit = TRUE)                                # Fit the model
sum(x$fit >= 5) / length(x$fit)                                 # Percent large expected counts
pchisq(x$lrt, x$df, lower.tail = FALSE)                         # p-value

# Exact log-linear
x <- loglinear(no23, t)                                         # Sample the fiber
x$p.value[["PR"]]                                               # p-value
count_tables(t, hmat(dim(t), no23))                             # Count tables in fiber

x <- loglinear(no3w, t)                                         # Sample the fiber
x$p.value[["PR"]]                                               # p-value
count_tables(t, hmat(dim(t), no3w))                             # Count tables in fiber

# Asymptotic logistic regression
u <- as_tibble(t) |> 
    mutate(across(everything(), as.numeric)) |> 
    group_by(ability, group) |> 
    summarize(p = sum(keep(n, response == 1)) / sum(n), 
              n = sum(n), 
              .groups = "drop")                                 # Transform

x <- glm(p ~ ability, 
         binomial, weights = n, data = u)                       # Fit model
sum(mle(t, no23, "lrm")$expected >= 5) / length(t)              # Percent large expected counts
pchisq(x$deviance, x$df.residual, lower.tail = FALSE)           # p-value

x <- glm(p ~ ability + group, 
         binomial, weights = n, data = u)                       # Fit model
sum(mle(t, no3w, "lrm")$expected >= 5) / length(t)              # Percent large expected counts
pchisq(x$deviance, x$df.residual, lower.tail = FALSE)           # p-value

x <- glm(p ~ ability + group + ability:group, 
         binomial, weights = n, data = u)                       # Fit model
sum(mle(t, full, "lrm")$expected >= 5) / length(t)              # Percent large expected counts
pchisq(x$deviance, x$df.residual, lower.tail = FALSE)           # p-value

# Exact logistic regression
d <- map(dimnames(t), as.numeric)                               # Dimension names
aux <- expand_grid(intercept = 1, 
                   ability = d$ability, 
                   group = d$group) |> 
    mutate(interaction = ability * group) |> 
    t()                                                         # Auxiliary matrix
v <- as_tibble(t) |> 
    arrange(desc(response), ability, group) |> 
    pull(n)                                                     # Table as a vector
pr <- computeUProbsCpp(matrix(v))                               # Probability of table

config_no23 <- lawrence(aux[1:2, ])                             # Configuration matrix
moves <- markov(config_no23, p = "arb")                         # Markov basis for fiber
sample <- metropolis(v, moves, 
                     iter = 10000, 
                     burn = 1000, 
                     thin = 10)$steps                           # Sample fiber
mean(computeUProbsCpp(sample) <= pr)                            # p-value
count_tables(v, config_no23)                                    # Count tables in fiber

config_no3w <- lawrence(aux[1:3, ])                             # Configuration matrix
moves <- markov(config_no3w, p = "arb")                         # Markov basis for fiber
sample <- metropolis(v, moves, 
                     iter = 10000, 
                     burn = 1000, 
                     thin = 10)$steps                           # Sample fiber
mean(computeUProbsCpp(sample) <= pr)                            # p-value
count_tables(v, config_no3w)                                    # Count tables in fiber

config_full <- lawrence(aux)                                    # Configuration matrix
moves <- markov(config_full, p = "arb")                         # Markov basis for fiber
sample <- metropolis(v, moves, 
                     iter = 10000, 
                     burn = 1000, 
                     thin = 10)$steps                           # Sample fiber
mean(computeUProbsCpp(sample) <= pr)                            # p-value
count_tables(v, config_no3w)                                    # Count tables in fiber

# ==============================================================================
# All Items
# ==============================================================================

# MCMC parameters
iter <- 10000
burn <- 1000
thin <- 10

# Significance level
alpha <- 0.05

# Helper function: for item x and with k ability levels, construct the relevant
# three-way contingency table and perform table_task on the resulting table.
# It computes Markov moves if needed. If k is large, this computation can take
# a little while.
lambda <- function(x, k) {
    # Read moves if possible, otherwise compute and save
    moves_file <- str_c("hci-moves", k, ".rds")
    dimnames <- list(ability = 0:(k - 1), group = 0:1, response = 0:1)
    if (file.exists(moves_file)) {
        moves <- readRDS(moves_file)
    } else {
        moves <- markov_moves(dimnames)
        saveRDS(moves, moves_file)
    }
    
    # Construct table
    t <- tibble(ability = discretize(HCI$total, k), 
                group = HCI[[grouper]], 
                response = HCI[[paste("Item", x)]]) |> 
        table()
    
    # Compute exact p-value for full logistic regression model
    v <- vectorize(t, "lrm")
    pr <- computeUProbsCpp(matrix(v))
    sample <- metropolis(v, moves$lrm_full, iter, burn, thin)$steps |> 
        suppressMessages()
    lrm_exc_full_p <- mean(computeUProbsCpp(sample) <= pr)
    
    # Perform table tasks
    bind_cols(table_task(t, moves, iter, burn, thin, alpha), 
              lrm_asy_full_p = mle(t, full, "lrm")$p, 
              lrm_exc_full_p = lrm_exc_full_p)
}

# Run the above helper function on all items with 6 and 9 ability levels. 
# Note that there are a lot of Markov chains to generate, so this computation
# will take a little bit of time (even after Markov moves have been computed).
df <- expand_grid(k = c(6L, 9L), item = 1:20) |> 
    mutate(out = pmap(list(item, k), lambda)) |> 
    unnest(out) |> 
    arrange(k, item)

# Add adjusted p-values for multiple testing. Group p-values into groups of 20.
# In other words, there is one group for each tuple (binning, strategy, model)
# and each group contains a p-value for each item. Use Benjamini-Hochberg
# adjustments on each group.
df <- df |> 
    group_by(k) |> 
    mutate(across(ends_with("_p"), 
                  \(x) p.adjust(x, method = "BH"), 
                  .names = "{.col}_adj")) |> 
    ungroup()

# Helper function for TeX output: Removes k labels after the first item
ability_label_remover <- function(k, item) {
    ifelse(item == 1, as.character(k), "")
}

# Helper function for TeX output: Encloses the first string x in \emph{...}
# if the first and second strings differ
empher <- function(x, y) {
    ifelse(x == y, x, str_c("\\emph{", x, "}"))
}

# Helper function for TeX output: Adds decorations to indicate various
# subtleties that should be considered with the DIF analysis
decorator <- function(x, no23_p_adj, full_p = 1) {
    # Make a list of decorations to add
    decorations <- character()
    if (!is.na(no23_p_adj) && 
               !str_detect(x, "none") && 
               no23_p_adj >= alpha) {
        # If DIF was detected, but the p-value for the no23 model becomes large
        # after BH adjustments, it means that the detected DIF can probably be
        # ignored and we mark it with a dagger
        decorations <- c(decorations, "\\textdagger")
    } else if (full_p < alpha) {
        # If full_p is small, it means that the full model does not fit well, 
        # so DIF analysis results should be taken with a grain of salt and
        # we mark it with an asterisk
        decorations <- c(decorations, "*")
    }
    
    # Add decorations and return
    decorations <- str_flatten_comma(decorations)
    if (decorations != "") {
        decorations <- str_c("\\textsuperscript{", decorations, "}")
    }
    str_c(x, decorations)
}

# Helper function for TeX output: Unabbreviates result strings
unabbreviator <- function(x) {
    case_match(x, 
               "none" ~ "none", 
               "unif" ~ "uniform", 
               "nonunif" ~ "nonuniform", 
               "fail" ~ "failure", 
               "unclass" ~ "unclassifiable")
}

# Tibble for TeX output
df_tex <- df |> 
    mutate(across(ends_with("_main"), unabbreviator)) |> 
    mutate(k = map2_chr(k, item, ability_label_remover)) |> 
    mutate(llm_asy = map2_chr(llm_asy_main, llm_exc_main, empher)) |> 
    mutate(llm_exc = map2_chr(llm_exc_main, llm_asy_main, empher)) |> 
    mutate(lrm_asy = map2_chr(lrm_asy_main, lrm_exc_main, empher)) |> 
    mutate(lrm_exc = map2_chr(lrm_exc_main, lrm_asy_main, empher)) |> 
    mutate(llm_asy = pmap_chr(list(llm_asy, llm_asy_no23_p_adj), decorator)) |> 
    mutate(llm_exc = pmap_chr(list(llm_exc, llm_exc_no23_p_adj), decorator)) |> 
    mutate(lrm_asy = pmap_chr(list(lrm_asy, 
                                   lrm_asy_no23_p_adj, 
                                   lrm_asy_full_p_adj), decorator)) |> 
    mutate(lrm_exc = pmap_chr(list(lrm_exc, 
                                   llm_exc_no23_p_adj, 
                                   lrm_exc_full_p_adj), decorator)) |> 
    select(k, item, llm_asy, llm_exc, lrm_asy, lrm_exc)

# Preview
view(df_tex)

# TeX header for table
header <- "& & \\multicolumn{2}{l}{Log-Linear Models} & \\multicolumn{2}{l}{Logistic Regressions} \\\\
\\cmidrule{3-6}
$\\#\\cA$ & Item & Asymptotic & Exact & Asymptotic & Exact \\\\"

# Save TeX output to file
df_tex |> 
    xtable(align = c(rep("c", 3), rep("l", 4))) |> 
    print(floating = FALSE, 
          include.rownames = FALSE, 
          include.colnames = FALSE, 
          booktabs = TRUE, 
          add.to.row = list(pos = list(0), command = header), 
          hline.after = c(-1, 0, 20, 40), 
          sanitize.text.function = identity, 
          comment = FALSE, 
          file = "../paper/table-hci.tex")
