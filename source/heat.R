setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("base.R")
library(latex2exp)  # For LaTeX in ggplot
library(ordr)       # For vectors in ggplot
library(colorspace) # For colors

# The following color is qualitative_hcl(3, palette = "set2")[1]
colors <- c("lightgray", "#ED90A4")

# Read dichotomous models
models <- readRDS("models.rds") |> 
    keep(\(m) m$is_dichotomous)

# Helper function that takes a model m and (x, y) in R^2, maps (x, y) into 
# (C^perp)^2 and then computes the joint distribution that would result by 
# taking that vector in (C^perp)^2 as a DIF vector
lambda <- function(m, x, y) {
    dif <- c(to_constants_perp(x), to_constants_perp(y))
    joint_distribution(m, dif)
}

# Tibble to record data about models
model_df <- tibble(model = models) |>
    mutate(seed = map_dbl(model, \(m) m$seed)) |> 
    mutate(I = map_dbl(model, \(m) length(m$dimnames$ability))) |> 
    mutate(x = map_dbl(model, \(m) m$uv[2]), 
           y = map_dbl(model, \(m) m$uv[4])) |> 
    mutate(unif_x = map_dbl(model, \(m) m$unif[2]), 
           unif_y = map_dbl(model, \(m) m$unif[4])) |> 
    mutate(nonunif_x = map_dbl(model, \(m) m$nonunif[2]), 
           nonunif_y = map_dbl(model, \(m) m$nonunif[4]))

# Mesh to use on each acix
mesh <- seq(-2.5, 2.5, by = 0.25)

# Compute grid of Kullback-Leibler divergences
df <- expand_grid(model_df, u_prime = mesh, v_prime = mesh) |>
    mutate(dist = pmap(list(model, u_prime, v_prime), lambda)) |> 
    mutate(KL = map_dbl(dist, dif_size))

# Make plot
x <- df |> 
    ggplot(aes(u_prime, v_prime, fill = KL)) + 
    facet_wrap(~ seed, nrow = 2) + 
    geom_tile() + 
    geom_vector(aes(x = 0.5*x, y = 0.5*y), color = "gray") + 
    geom_vector(aes(x = unif_x, y = unif_y)) + 
    geom_vector(aes(x = nonunif_x, y = nonunif_y)) + 
    ylab(TeX("$v_\\Delta$")) +
    xlab(TeX("$u_\\Delta$")) +
    scale_fill_gradient(low = colors[1], high = colors[2]) + 
    theme_minimal() + 
    theme(line = element_blank())

# Visual test
x

ggsave("../figs/heat.svg", x, 
       width = 10, height = 4, units = "in", 
       bg = "white")
