setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(colorspace)
library(xtable)

################################################################################
# Read Files
################################################################################

seeds_file <- "seeds.csv"
models_file <- "models.rds"
results_file <- "results.csv"

seeds <- read_csv(seeds_file, col_types = "ic") |> 
    pull(seed)
models <- readRDS(models_file)
results <- read_csv(results_file)

################################################################################
# Themes and Colors
################################################################################

theme_set(theme_minimal())

short_names <- c("mle", 
                 "heuristic", 
                 "asy_main", 
                 "asy_swap", 
                 "asy_comp", 
                 "exc_main", 
                 "exc_swap")
colors <- c("#BBBBBB", "#BBBBBB", qualitative_hcl(5, palette = "set2"))
alphas = c(1, 0.5, rep(1, 5))

################################################################################
# Some Helper Functions
################################################################################

# Takes as input a "true" vector, ie, a vector of the true DIF status of some
# tables, and a "result" vector which indicates the output of some method
# of DIF analysis. Returns a logical vector which indicates whether the method
# succeeded at correctly *detecting* DIF. 
detect <- function(true, result) {
    map2_lgl(true, 
             result, 
             \(x, y) if (x == dif_test_results$none) { y == dif_test_results$none } 
             else { y %in% c(dif_test_results$uniform, 
                             dif_test_results$nonuniform, 
                             dif_test_results$unclassifiable) })
}

# Like above, but now returns a logical vector which indicates whether the 
# method succeeded at correctly *classifying* DIF. 
classify <- function(true, result) {
    map2_lgl(true, 
             result, 
             \(x, y) if (x == dif_test_results$none) { NA } else { x == y })
}

# For labeling plots
step_labels <- c(detect = "Detect",
                 classify = "Classify")
type_labels <- c(none = "No DIF",
                 unif = "Uniform DIF",
                 nonunif = "Nonuniform DIF")
norm_labels <- c("0" = "Norm 0",
                 "1" = "Norm 1",
                 "2" = "Norm 2")
thing_labels <- c(mle = "MLE Exists",
                  heuristic = "Heuristic Passes",
                  asy_main = "Asymptotic Succeeds",
                  exc_main = "Exact Succeeds",
                  asy_swap = "Swapped Asymptotic Succeeds",
                  exc_swap = "Swapped Exact Succeeds",
                  asy_comp = "Asymptotic with Comparison Succeeds")
family_labels <- c(llm = "Log-Linear Models",
                   lrm = "Logistic Regressions")

names(colors) <- map_chr(short_names, \(x) thing_labels[[x]])
names(alphas) <- map_chr(short_names, \(x) thing_labels[[x]])

# For naming plots
main_method_pair <- c("asy_main", "exc_main")
method_pairs <- list(main_method_pair,
                     c("asy_main", "asy_swap"),
                     c("exc_main", "exc_swap"),
                     c("asy_main", "asy_comp"))

method_pair_suffix <- function(x) {
    if (setequal(x, method_pairs[[1]])) {
        out <- ""
    } else if (setequal(x, method_pairs[[2]])) {
        out <- "-swap-asy"
    } else if (setequal(x, method_pairs[[3]])) {
        out <- "-swap-exc"
    } else if (setequal(x, method_pairs[[4]])) {
        out <- "-comp-asy"
    }
    return(out)
}

filename <- function(s, 
                     family, 
                     methods = main_method_pair, 
                     format = "svg") {
    str_c(s, "-", family, method_pair_suffix(methods), ".", format)
}

################################################################################
# Main Plotting Function
################################################################################

plot <- function(s, 
                 family, 
                 methods = main_method_pair,
                 filename = NULL) {
    # Return NULL if requesting LRM plot for non-dichotomous seed
    if (family == "lrm" && !models[[match(s, seeds)]]$is_dichotomous) {
        return("")
    }
    
    # Relevant column names
    method_cols <- str_c(family, "_", methods)
    mle_cols <- str_c(family, "_asy_", c("no23", "no3w"), "_mle")
    heuristic_cols <- str_c(family, "_asy_", c("no23", "no3w"), "_heuristic")
    
    # Shrink tibble
    df <- results |> 
        filter(seed == s) |> 
        select(type, norm, sample_size, 
               detect_mle = mle_cols[1], 
               classify_mle = mle_cols[2],
               detect_heuristic = heuristic_cols[1],
               classify_heuristic = heuristic_cols[2],
               all_of(method_cols))
    
    # Compare results
    for (m in method_cols) {
        df <- df |> 
            mutate(detect(df$type, df[[m]])) |> 
            rename_with(\(x) str_c("detect_", m), last_col()) |> 
            mutate(classify(df$type, df[[m]])) |> 
            rename_with(\(x) str_c("classify_", m), last_col())
            
    }
    
    # Aggregate
    df <- df |> 
        select(-all_of(method_cols)) |> 
        group_by(type, norm, sample_size) |> 
        summarize(across(everything(), mean), .groups = "keep") |> 
        pivot_longer(cols = -group_cols(), 
                     names_to = "thing", 
                     values_to = "proportion") |> 
        filter(!is.na(proportion)) |> 
        ungroup() |> 
        mutate(step = map_chr(thing, \(x) 
                              if (str_detect(x, "detect_")) { "detect" } 
                              else { "classify" })) |> 
        mutate(thing = str_replace(thing, "detect_", "")) |> 
        mutate(thing = str_replace(thing, "classify_", "")) |> 
        mutate(thing = str_replace(thing, str_c(family, "_"), ""))
    
    # Remove MLE/heuristic if no asymptotics are being plotted
    include_mle <- any(str_detect(methods, "asy"))
    if (!include_mle) {
        df <- filter(df, !(thing %in% c("mle", "heuristic")))
    }
    
    # Clean up labels
    df <- df |> 
        mutate(type = factor(type, 
                             levels = names(type_labels), 
                             labels = unname(type_labels),
                             ordered = TRUE)) |> 
        mutate(norm = factor(norm, 
                             levels = names(norm_labels), 
                             labels = unname(norm_labels),
                             ordered = TRUE)) |> 
        mutate(thing = factor(thing, 
                              levels = names(thing_labels), 
                              labels = unname(thing_labels),
                              ordered = TRUE)) |> 
        mutate(step = factor(step, 
                             levels = names(step_labels), 
                             labels = unname(step_labels),
                             ordered = TRUE))
    
    # Make plot    
    out <- df |> 
        ggplot(aes(x = sample_size, 
                   y = proportion, 
                   color = thing,
                   alpha = thing)) + 
        facet_grid(step ~ type + norm, 
                   labeller = labeller(Helper = function(x) { rep("", length(x)) })) + 
        geom_line() +
        geom_hline(yintercept = 0.95, linetype = "dotted", color="lightgray") + 
        expand_limits(y = c(0, 1)) + 
        scale_x_continuous(breaks = c(0, 1000, 2000)) + 
        scale_colour_manual(values = colors, name = "") + 
        scale_alpha_manual(values = alphas, name = "") + 
        theme(legend.position="bottom") + 
        xlab("Sample Size") + 
        ylab("Proportion") + 
        ggtitle(str_c(family_labels[[family]], ", Seed ", s))
    
    # Return or save
    if (is_null(filename)) {
        return(out)
    } else {
        ggsave(filename, out, 
               width = 10, height = 6, units = "in", 
               bg = 'white')
        return(filename)
    }
}

################################################################################
# Some Visual Tests
################################################################################

plot(1948, "llm")
plot(1947, "lrm", method_pairs[[2]])
plot(1947, "lrm", method_pairs[[3]])
plot(1948, "llm", method_pairs[[4]])

################################################################################
# Make and save PDF plots for the paper
################################################################################

paper_folder <- "../paper/"
tibble(seed = c(1947, 1947, 1948, 1930), 
       family = c("llm", "lrm", "llm", "lrm")) |> 
    mutate(filename = map2_chr(seed, family,
                               \(x, y) filename(x, y, format = "pdf"))) |> 
    mutate(filename = str_c(paper_folder, filename)) |> 
    mutate(filename = pmap_chr(list(seed, family, filename), 
                               \(x, y, z) plot(x, y, filename = z)))

################################################################################
# Make and save all plots
################################################################################

fig_folder <- "../figs/" 
grid <- expand_grid(seed = seeds, 
            family = c("llm", "lrm"),
            method_pair = method_pairs) |> 
    mutate(filename = pmap_chr(list(seed, family, method_pair), filename)) |> 
    mutate(filename = str_c(fig_folder, filename)) |> 
    mutate(filename = pmap_chr(list(seed, family, method_pair, filename), 
                               plot))

################################################################################
# Reorganize for output
################################################################################

# Column width for seed column as percentage of total page width
seed_col_width <- 4

cell_maker <- function(filename) {
    if (filename == "") {
        return("")
    }
    x <- str_replace(filename, "../", "")
    str_c("<img src=\"./", x, "\" />")
}

grid <- grid |> 
    mutate(cell = map_chr(filename, cell_maker)) |> 
    mutate(name = str_c(family, map_chr(method_pair, method_pair_suffix))) |> 
    select(seed, name, cell)

primary <- grid |> 
    filter(name %in% c("llm", "lrm")) |> 
    pivot_wider(names_from = name, values_from = cell)

secondary <- grid |> 
    filter(!(name %in% c("llm", "lrm"))) |> 
    pivot_wider(names_from = name, values_from = cell)

################################################################################
# Primary Comparison Plots
################################################################################

titles <- c("Seed",
            "Log-Linear Models",
            "Logistic Regressions")
n <- length(titles)
primary <- primary |> 
    rename_with(\(.) titles)

primary_output <- primary |> 
    xtable(align = rep("l", n+1)) |> 
    print(type = "html",
          include.rownames = FALSE,
          sanitize.text.function = identity,
          html.table.attributes = "style=\"width: 100%;\"",
          comment = FALSE,
          print.results = FALSE) |> 
    str_replace_all("\n  ", "\n") |> 
    str_replace_all("\n ", "\n") |> 
    str_replace("<th> Seed", str_c("<th style=\"width: ", 
                                   seed_col_width, "%;\"> Seed")) |> 
    str_replace_all("<th>", str_c("<th style=\"width: ", 
                                  (100 - seed_col_width)/(n-1), 
                                  "%;\">"))

# Check
cat(primary_output)

################################################################################
# Secondary Comparison Plots
################################################################################

titles <- c("Seed",
            "LLM Asymptotic Swapped",
            "LLM Exact Swapped",
            "LLM Asymptotic Comparison",
            "LRM Asymptotic Swapped",
            "LRM Exact Swapped",
            "LRM Asymptotic Comparison")
n <- length(titles)
secondary <- secondary |> 
    rename_with(\(.) titles)

header <- str_c("<tr> ",
                "<th rowspan=2 style=\"width: ", 
                seed_col_width, 
                "%;\"> Seed </th> ",
                "<th colspan=3> Log-Linear Models </th> ",
                "<th colspan=3> Logistic Regressions </th> ",
                "</tr>\n")

secondary_output <- secondary |> 
    xtable(align = rep("l", n+1)) |> 
    print(type = "html",
          include.rownames = FALSE,
          sanitize.text.function = identity,
          html.table.attributes = "style=\"width: 100%\"",
          comment = FALSE,
          print.results = FALSE) |> 
    str_replace_all("\n  ", "\n") |> 
    str_replace_all("\n ", "\n") |> 
    str_replace("\">\n", str_c("\">\n", header)) |>
    str_replace_all("LLM ", "") |> 
    str_replace_all("LRM ", "") |> 
    str_replace("<th> Seed </th>", "") |> 
    str_replace_all("<th>", str_c("<th style=\"width: ", 
                                  (100 - seed_col_width)/(n-1), 
                                  "%;\">"))

# Check
cat(secondary_output)

################################################################################
# Make README
################################################################################

readme_source <- "README.md"
primary_placeholder <- "<!--PRIMARY-->"
secondary_placeholder <- "<!--SECONDARY-->"

readme_output <- "../README.md"

note <- "<!-- 
Do not edit README.md in the root directory of this repository! 

This file is generated using:

* source/README.md for the text and structure, and 
* source/sim-plot.R for the tables of plots. 

To make changes, edit the above instead. 
-->\n\n"

out <- read_file(readme_source) |> 
    str_replace(primary_placeholder, primary_output) |> 
    str_replace(secondary_placeholder, secondary_output)
out <- str_c(note, out)
write_file(out, readme_output)
