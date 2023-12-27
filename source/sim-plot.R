# This script file generates plots of the results produced by sim-run.csv.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(colorspace)
library(xtable)

theme_set(theme_minimal())

# ==============================================================================
# Read Files
# ==============================================================================

seeds_file <- "seeds.csv"
models_file <- "models.rds"
results_file <- "results.csv"

seeds <- read_csv(seeds_file, col_types = "ic") |> 
    pull(seed)
models <- readRDS(models_file)
results <- read_csv(results_file)

# ==============================================================================
# Helper Functions
# ==============================================================================

# Takes as input a "true" vector, ie, a vector of the true DIF status of some
# tables, and a "result" vector which indicates the output of some method
# of DIF analysis. Returns a logical vector which indicates whether the method
# succeeded at correctly *detecting* DIF. We first define a helper function, 
# and then a vectorized version.
detect_helper <- function(true, result) {
    ifelse(true == dif_test_results$none, 
           result == dif_test_results$none, 
           result %in% c(dif_test_results$uniform, 
                         dif_test_results$nonuniform, 
                         dif_test_results$unclassifiable))
}

detect <- function(true, result) {
    map2_lgl(true, result, detect_helper)
}

# Like above, but now returns a logical vector which indicates whether the
# method succeeded at correctly *classifying* DIF. Again, we first define a
# helper function, and then a vectorized version.
classify_helper <- function(true, result) {
    ifelse(true == dif_test_results$none, NA, true == result)
}

classify <- function(true, result) {
    map2_lgl(true, result, classify_helper)
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
                  asy_swap = "Swapped Asymptotic Succeeds", 
                  asy_comp = "Asymptotic with Comparison Succeeds",
                  exc_main = "Exact Succeeds", 
                  exc_swap = "Swapped Exact Succeeds")

family_labels <- c(llm = "Log-Linear Models", 
                   lrm = "Logistic Regressions")

# Pairs of DIF analysis methods to compare in plots
method_pairs <- list(c("asy_main", "exc_main"), 
                     c("asy_main", "asy_swap"), 
                     c("exc_main", "exc_swap"), 
                     c("asy_main", "asy_comp"))

main_method_pair <- method_pairs[[1]]

# For naming plots
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

# ==============================================================================
# Line Styling
# ==============================================================================

# For arXiv and GitHub
colors <- c("#BBBBBB", "#BBBBBB", qualitative_hcl(5, palette = "set2")) |> 
    set_names(thing_labels)
alphas <- c(1, 0.5, rep(1, 5)) |>
    set_names(thing_labels)

# For journal: variant 1
linetypes_var <- c("dashed", "dashed", "solid", NA, NA, "solid", NA) |> 
    set_names(thing_labels)
alphas_var <- c(0.5, 0.2, 0.5, NA, NA, 1, NA) |> 
    set_names(thing_labels)

# For journal: variant 2
colors_var <- c("#BBBBBB", "#BBBBBB", "#0d0887", NA, NA, "#cc4778", NA) |> 
    set_names(thing_labels)

# ==============================================================================
# Main Plotting Functions
# ==============================================================================

# Build a relevant data frame for plotting
plot_df <- function(s, family, methods = main_method_pair) {
    # Return NULL if requesting LRM plot for non-dichotomous seed
    if (family == "lrm" && !models[[match(s, seeds)]]$is_dichotomous) {
        return(NULL)
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
        mutate(step = map_chr(thing, \(x) ifelse(str_detect(x, "detect_"), 
                                                 "detect", 
                                                 "classify"))) |> 
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
    
    # Return
    df
}
                    
# Main plotting function (for arXiv and GitHub, uses color)
plot <- function(s, family, methods = main_method_pair, filename = NULL) {
    # Build data frame
    df <- plot_df(s, family, methods)
    
    # Return empty string if requesting LRM plot for non-dichotomous seed
    if (is_null(df)) {
        return("")
    }
    
    # Make plot
    out <- df |> 
        ggplot(aes(x = sample_size, 
                   y = proportion, 
                   color = thing, 
                   alpha = thing)) + 
        facet_grid(step ~ type + norm, 
                   labeller = labeller(Helper = \(x) rep("", length(x)))) + 
        geom_hline(yintercept = 0.95, 
                   linetype = "dotted", 
                   color = "lightgray") + 
        geom_line() + 
        expand_limits(y = c(0, 1)) + 
        scale_x_continuous(breaks = c(0, 1000, 2000)) + 
        scale_colour_manual(values = colors, name = "") + 
        scale_alpha_manual(values = alphas, name = "") + 
        theme(legend.position = "bottom") + 
        xlab("Sample Size") + 
        ylab("Proportion") + 
        ggtitle(str_c(family_labels[[family]], ", Seed ", s))
    
    # Return or save
    if (is_null(filename)) {
        return(out)
    } else {
        ggsave(filename, out, 
               width = 10, height = 6, units = "in", 
               bg = "white")
        return(filename)
    }
}

# Variant plotting function for grayscale plots
plot_var1 <- function(s, family, filename = NULL) {
    # Build data frame
    df <- plot_df(s, family)
    
    # Return empty string if requesting LRM plot for non-dichotomous seed
    if (is_null(df)) {
        return("")
    }
    
    # Make plot
    out <- df |> 
        ggplot(aes(x = sample_size, 
                   y = proportion,
                   linetype = thing,
                   alpha = thing)) + 
        facet_grid(step ~ type + norm, 
                   labeller = labeller(Helper = \(x) rep("", length(x)))) + 
        geom_hline(yintercept = 0.95, 
                   linetype = "dotted", 
                   color = "lightgray") + 
        geom_line() + 
        expand_limits(y = c(0, 1)) + 
        scale_x_continuous(breaks = c(0, 1000, 2000)) + 
        scale_linetype_manual(values = linetypes_var, name = "") + 
        scale_alpha_manual(values = alphas_var, name = "") + 
        theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
        xlab("Sample Size") + 
        ylab("Proportion") + 
        ggtitle(str_c(family_labels[[family]], ", Seed ", s))
    
    # Return or save
    if (is_null(filename)) {
        return(out)
    } else {
        ggsave(filename, out, 
               width = 10, height = 6, units = "in", 
               bg = "white")
        return(filename)
    }
}

# Variant plotting function for color plots that are readable in grayscale
plot_var2 <- function(s, family, filename = NULL) {
    # Build data frame
    df <- plot_df(s, family)
    
    # Return empty string if requesting LRM plot for non-dichotomous seed
    if (is_null(df)) {
        return("")
    }
    
    # Make plot
    out <- df |> 
        ggplot(aes(x = sample_size, 
                   y = proportion, 
                   color = thing, 
                   alpha = thing)) + 
        facet_grid(step ~ type + norm, 
                   labeller = labeller(Helper = \(x) rep("", length(x)))) + 
        geom_hline(yintercept = 0.95, 
                   linetype = "dotted", 
                   color = "lightgray") + 
        geom_line() + 
        expand_limits(y = c(0, 1)) + 
        scale_x_continuous(breaks = c(0, 1000, 2000)) + 
        scale_colour_manual(values = colors_var, name = "") + 
        scale_alpha_manual(values = alphas, name = "") + 
        theme(legend.position = "bottom") + 
        xlab("Sample Size") + 
        ylab("Proportion") + 
        ggtitle(str_c(family_labels[[family]], ", Seed ", s))
    
    # Return or save
    if (is_null(filename)) {
        return(out)
    } else {
        ggsave(filename, out, 
               width = 10, height = 6, units = "in", 
               bg = "white")
        return(filename)
    }
}

# ==============================================================================
# Some Visual Tests
# ==============================================================================

plot(1948, "llm")
plot(1947, "lrm", method_pairs[[2]])
plot(1947, "lrm", method_pairs[[3]])
plot(1948, "llm", method_pairs[[4]])

plot_var1(1947, "llm")
plot_var1(1948, "llm")

plot_var2(1947, "llm")
plot_var2(1948, "llm")

# ==============================================================================
# Make Plots for Paper (arXiv)
# ==============================================================================

paper_folder <- "../paper/arXiv/"

tibble(seed = c(1947, 1947, 1948, 1930), 
       family = c("llm", "lrm", "llm", "lrm")) |> 
    mutate(filename = map2_chr(seed, family, 
                               \(x, y) filename(x, y, format = "pdf"))) |> 
    mutate(filename = str_c(paper_folder, filename)) |> 
    mutate(filename = pmap_chr(list(seed, family, filename), 
                               \(x, y, z) plot(x, y, filename = z)))

# ==============================================================================
# Make Plots for Paper (Journal)
# ==============================================================================

paper_folder <- "../paper/journal/"

tibble(seed = c(1947, 1947, 1948, 1930), 
       family = c("llm", "lrm", "llm", "lrm")) |> 
  mutate(filename = map2_chr(seed, family, 
                             \(x, y) filename(x, y, format = "pdf"))) |> 
  mutate(filename = str_c(paper_folder, filename)) |> 
  mutate(filename = pmap_chr(list(seed, family, filename), 
                             \(x, y, z) plot_var2(x, y, filename = z)))

# ==============================================================================
# All Plots for GitHub
# ==============================================================================

fig_folder <- "../figs/"
grid_file <- "tests/grid.rds"

# Uncomment (Ctrl+Shift+C) to regenerate; this takes a few minutes!
# grid <- expand_grid(seed = seeds,
#                     family = c("llm", "lrm"),
#                     method_pair = method_pairs) |>
#     mutate(filename = pmap_chr(list(seed, family, method_pair), filename)) |>
#     mutate(filename = str_c(fig_folder, filename)) |>
#     mutate(filename = pmap_chr(list(seed, family, method_pair, filename),
#                                plot))
# saveRDS(grid, grid_file)

# ==============================================================================
# README Table Output Prep
# ==============================================================================

if (file.exists(grid_file)) {
    grid <- readRDS(grid_file)
}

# Restructure the grid
grid <- grid |> 
    mutate(name = str_c(family, map_chr(method_pair, method_pair_suffix))) |> 
    select(seed, name, filename)

# Turns filename into HTML image tag
cell_maker <- function(filename, img_width) {
    if (filename == "") {
        return("")
    }
    x <- str_replace(filename, "../", "")
    str_c("<img width=\"", img_width, "\" src=\"./", x, "\" />")
}

# ==============================================================================
# README Table for Primary Comparison Plots
# ==============================================================================

titles <- c("Seed", 
            "Log-Linear Models", 
            "Logistic Regressions")

primary <- grid |> 
    filter(name %in% c("llm", "lrm")) |> 
    mutate(filename = map_chr(filename, \(x) cell_maker(x, 300))) |> 
    pivot_wider(names_from = name, values_from = filename) |> 
    rename_with(\(.) titles)

primary_output <- primary |> 
    xtable(align = rep("l", ncol(primary) + 1)) |> 
    print(type = "html", 
          include.rownames = FALSE, 
          sanitize.text.function = identity, 
          html.table.attributes = "", 
          comment = FALSE, 
          print.results = FALSE) |> 
    str_replace_all("\n  ", "\n") |> 
    str_replace_all("\n ", "\n") |> 
    str_replace("<table >", "<table>") |> 
    str_replace_all("<th>", "<th align=\"center\">") |> 
    str_replace_all("<td>", "<td align=\"center\">")

# Check
cat(primary_output)

# ==============================================================================
# README Table for Secondary Comparison Plots
# ==============================================================================

titles <- c("Seed", 
            "LLM Asy Swap", 
            "LLM Exc Swap", 
            "LLM Asy Comp", 
            "LRM Asy Swap", 
            "LRM Exc Swap", 
            "LRM Asy Comp")

secondary <- grid |> 
    filter(!(name %in% c("llm", "lrm"))) |> 
    mutate(filename = map_chr(filename, \(x) cell_maker(x, 95))) |> 
    pivot_wider(names_from = name, values_from = filename) |> 
    rename_with(\(.) titles)

header <- str_c("<tr> ", 
                "<th rowspan=2> Seed </th> ", 
                "<th colspan=3> Log-Linear Models </th> ", 
                "<th colspan=3> Logistic Regressions </th> ", 
                "</tr>\n")

secondary_output <- secondary |> 
    xtable(align = rep("l", ncol(secondary) + 1)) |> 
    print(type = "html", 
          include.rownames = FALSE, 
          sanitize.text.function = identity, 
          html.table.attributes = "", 
          comment = FALSE, 
          print.results = FALSE) |> 
    str_replace_all("\n  ", "\n") |> 
    str_replace_all("\n ", "\n") |> 
    str_replace("<table >\n", str_c("<table>\n", header)) |> 
    str_replace("<th> Seed </th>", "") |> 
    str_replace_all("LLM ", "") |> 
    str_replace_all("LRM ", "") |> 
    str_replace_all("<th>", "<th align=\"center\">") |> 
    str_replace_all("<td>", "<td align=\"center\">")

# Check
cat(secondary_output)

# ==============================================================================
# README Output
# ==============================================================================

readme_output <- "../README.md"

readme_source <- "README.md"
primary_placeholder <- "<!--PRIMARY-->"
secondary_placeholder <- "<!--SECONDARY-->"
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
