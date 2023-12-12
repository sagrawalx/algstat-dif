# This makes the TeX table in the paper that describes parameters chosen by the
# seeds in seeds.csv. It's set up to make a table just for the three seeds that
# we focus on in the paper, but it is easily modified to generate a longer 
# table for all 16 seeds. 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(xtable)

# Read models file
models_file <- "models.rds"
models <- readRDS(models_file)

# Filter models; comment this out to generate a table for all seeds!
select_seeds <- c(1930, 1947, 1948)
models <- models |> 
    keep(\(m) m$seed %in% select_seeds)

# Some helper functions to convert various types of parameters to strings
num_str <- function(x) {
    sprintf("%.3f", x)
}

set_str <- function(x) {
    str_c("\\{", str_flatten_comma(x), "\\}")
}

vec_str <- function(x) {
    str_c("(", str_flatten_comma(num_str(x)), ")")
}

# Make list of strings for cells in the table
s <- map(models, \(m) m$seed) |> 
    map_chr(as.character)
a <- map(models, \(m) m$dimnames$ability) |> 
    map_chr(set_str)
r <- map(models, \(m) m$dimnames$response) |> 
    map_chr(set_str)
p <- map(models, \(m) m$group_dist[1]) |> 
    map_chr(num_str)
q0 <- map(models, \(m) m$ability_dist[, 1]) |> 
    map_chr(vec_str)
q1 <- map(models, \(m) m$ability_dist[, 2]) |> 
    map_chr(vec_str)
uv <- map(models, \(m) m$uv) |> 
    map_chr(vec_str)
unif <- map(models, \(m) m$unif) |> 
    map_chr(vec_str)
nonunif <- map(models, \(m) m$nonunif) |> 
    map_chr(vec_str)

# Helper function to clear repeated entries in the seeds column
clear <- function(x, y) {
    ifelse(y, x, "")
}

# Make TeX tibble
df_tex <- tibble(s, p, a, q0, q1, r, uv, unif, nonunif) |> 
    pivot_longer(-s) |> 
    mutate(s = map2_chr(s, name == "p", clear)) |> 
    mutate(name = case_when(name == "p" ~ "$\\Pr[G = 0]$",
                            name == "a" ~ "$\\cA$", 
                            name == "q0" ~ "$\\Pr[A \\mid G = 0]$", 
                            name == "q1" ~ "$\\Pr[A \\mid G = 1]$",
                            name == "r" ~ "$\\cR$", 
                            name == "uv" ~ "$(u_0, v_0)$",
                            name == "unif" ~ "Uniform $(u_\\Delta, v_\\Delta)$",
                            name == "nonunif" ~ "Nonuniform $(u_\\Delta, v_\\Delta)$")) |> 
    rename(Seed = s, 
           Parameter = name, 
           Value = value)

# Check
view(df_tex)

# Output
df_tex |> 
    xtable(align = c(rep("l", ncol(df_tex) + 1))) |> 
    print(floating = FALSE, 
          include.rownames = FALSE, 
          include.colnames = TRUE, 
          booktabs = TRUE, 
          hline.after = c(-1, seq(0, nrow(df_tex), by = nrow(df_tex) / length(s))), 
          sanitize.text.function = identity, 
          comment = FALSE,
          file = "../paper/table-params.tex")
