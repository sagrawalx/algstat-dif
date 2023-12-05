setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(xtable)
library(ShinyItemAnalysis)

################################################################################
# Initialize
################################################################################

# Load data
data(HCI)

# Choose grouping variable: either "major" or "gender"
grouper <- "major"

# Discretizes a numerical vector x into discrete levels 0, 1, ..., k-1 
discretize <- function(x, k) {
    y <- cut_interval(x, n = k)
    as.integer(y) - 1
}

################################################################################
# Focus item
################################################################################

# Choose focus item: Item 1, Item 2, ..., Item 20
focus_item <- 17

# Construct three-way table
t <- tibble(ability = discretize(HCI$total, 6),
            group = HCI[[grouper]],
            response = HCI[[paste("Item", focus_item)]]) |> 
    table()
t

# Asymptotic log-linear
no23 <- list(c(1,2), c(1,3))                    # Facet specification of model
x <- loglin(t, no23, fit = TRUE)                # Fit the model
sum(x$fit >= 5) / length(x$fit)                 # Percent large expected counts
pchisq(x$lrt, x$df, lower.tail = FALSE)         # p-value

no3w <- list(c(1,2), c(1,3), c(2,3))            # Facet specification of model
x <- loglin(t, no3w, fit = TRUE)                # Fit the model
sum(x$fit >= 5) / length(x$fit)                 # Percent large expected counts
pchisq(x$lrt, x$df, lower.tail = FALSE)         # p-value

# Exact log-linear
x <- loglinear(no23, t)                         # Sample the fiber
x$p.value[["PR"]]                               # p-value
count_tables(t, hmat(dim(t), no23))             # Count tables in fiber

x <- loglinear(no3w, t)                         # Sample the fiber
x$p.value[["PR"]]                               # p-value
count_tables(t, hmat(dim(t), no3w))             # Count tables in fiber

# Asymptotic logistic regression
u <- as_tibble(t) |>
    mutate(across(everything(), as.numeric)) |>
    group_by(ability, group) |> 
    summarize(p = sum(keep(n, response == 1))/sum(n), 
              n = sum(n), 
              .groups = "drop")                                 # Transform

x <- glm(p ~ ability, binomial, weights = n, data = u)          # Fit model
sum(mle(t, no23, "lrm")$expected >= 5)/length(t)                # Percent large expected counts
pchisq(x$deviance, x$df.residual, lower.tail = FALSE)           # p-value

x <- glm(p ~ ability + group, binomial, weights = n, data = u)  # Fit model
sum(mle(t, no3w, "lrm")$expected >= 5)/length(t)                # Percent large expected counts
pchisq(x$deviance, x$df.residual, lower.tail = FALSE)           # p-value

# Exact logistic regression
d <- map(dimnames(t), as.numeric)                               # Dimension names
A <- t(expand_grid(1, d$ability, d$group))                      # Auxiliary matrix
v <- as_tibble(t) |>
    arrange(desc(response), ability, group) |>
    pull(n)                                                     # Table as a vector
pr <- computeUProbsCpp(matrix(v))                               # Probability of table

config_no23 <- lawrence(A[1:2,])                                # Configuration matrix
moves <- markov(config_no23, p = "arb")                         # Markov basis for fiber
sample <- metropolis(v, moves, 
                     iter = 10000, 
                     burn = 1000, 
                     thin = 10)$steps                           # Sample fiber
mean(computeUProbsCpp(sample) <= pr)                            # p-value
count_tables(v, config_no23)                                    # Count tables in fiber

config_no3w <- lawrence(A)                                      # Configuration matrix
moves <- markov(config_no3w, p = "arb")                         # Markov basis for fiber
sample <- metropolis(v, moves,
                     iter = 10000,
                     burn = 1000,
                     thin = 10)$steps                           # Sample fiber
mean(computeUProbsCpp(sample) <= pr)                            # p-value
count_tables(v, config_no3w)                                    # Count tables in fiber

################################################################################
# All Items
################################################################################

# Helper function: for item x and with k ability levels, construct the relevant
# three-way contingency table and perform table_task on the resulting table.
# It computes Markov moves if needed. If k is large, this computation can take
# a little while. 
lambda <- function(x, k) {
    # Read moves if possible, otherwise compute and save
    moves_file <- paste("hci-moves", k, ".rds", sep = "")
    dimnames <- list(ability = 0:(k-1), group = 0:1, response = 0:1)
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
    
    # Perform table tasks
    table_task(t, moves)
}

# Run the helper function on all items with 6 and 9 ability levels
df <- expand_grid(k = c(6L, 9L), item = 1:20) |> 
    mutate(out = pmap(list(item, k), lambda)) |>
    unnest(out) |> 
    arrange(k, item)

# Helper function for TeX output in the paper: encloses the first string x
# in a TeX \emph{...} if the first and second strings differ.
empher <- function(x, y) {
    if (x == y) {
        out <- x
    } else {
        out <- str_c("\\emph{", x, "}")
    }
}

# Tibble for TeX output
x <- df |> 
    select(k, item, llm_asy_main, llm_exc_main, lrm_asy_main, lrm_exc_main) |> 
    mutate(across(contains("_"), 
                  \(x) case_match(x, 
                                  "none" ~ "none", 
                                  "unif" ~ "uniform", 
                                  "nonunif" ~ "nonuniform", 
                                  "fail" ~ "failure",
                                  "unclass" ~ "unclassifiable"))) |> 
    mutate(k = map2_chr(k, item, 
                        \(k, i) if (i == 1) { as.character(k) } else { "" })) |> 
    mutate(llm_asy = map2_chr(llm_asy_main, llm_exc_main, empher)) |> 
    mutate(llm_exc = map2_chr(llm_exc_main, llm_asy_main, empher)) |> 
    mutate(lrm_asy = map2_chr(lrm_asy_main, lrm_exc_main, empher)) |> 
    mutate(lrm_exc = map2_chr(lrm_exc_main, lrm_asy_main, empher)) |> 
    select(-ends_with("_main"))

# Save TeX output to file; copied into the main paper and mildly edited there
header <- c("& & \\multicolumn{2}{l}{Log-Linear Models} & \\multicolumn{2}{l}{Logistic Regressions} \\\\",
            "\\cmidrule{3-6}",
            "$\\#\\cA$ & Item & Asymptotic & Exact & Asymptotic & Exact \\\\") |> 
    str_flatten(collapse = "\n")
xtable(x, align = c(rep("c",3), rep("l",4))) |> 
    print(floating = FALSE,
          include.rownames = FALSE, 
          include.colnames = FALSE, 
          booktabs = TRUE, 
          add.to.row = list(pos = list(0), command = header),
          hline.after = c(-1,0,20,40),
          sanitize.text.function = identity,
          comment = FALSE,
          file = "../paper/hci-table.tex")
