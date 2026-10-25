## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(sltree)

## -----------------------------------------------------------------------------
model <- latent_tree_models(1)
head(model)

## ----eval = FALSE-------------------------------------------------------------
# view_ltree(1)                        # direction = "LR" suits trees with many leaves

## ----fig.width = 7, fig.height = 5--------------------------------------------
preset_tree <- view_ltree(1, "plot")

## ----eval = FALSE-------------------------------------------------------------
# set.seed(1234)
# dat <- simulate_ltree_data(n = 1e5, model = model)

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 400
h <- rnorm(n)
gaussian_x <- as.data.frame(replicate(4, 0.8 * h + 0.6 * rnorm(n)))
names(gaussian_x) <- paste0("V", 1:4)
x <- data.frame(
  V1 = exp(gaussian_x$V1),
  V2 = gaussian_x$V2^3,
  V3 = as.integer(gaussian_x$V3 > 0),
  V4 = findInterval(gaussian_x$V4, c(-0.5, 0.5))
)
head(x)

true_edges <- data.frame(
  node = c(paste0("V", 1:4), "H1"),
  parent = c(rep("H1", 4), "")
)
true_types <- c("continuous", "continuous", "binary", "ordinal", "latent")
true_tree <- ltree_from_edge(true_edges, true_types)

## -----------------------------------------------------------------------------
fit <- learn_sltree(x, tau1 = 3, tau2 = 5, use_near_pd = FALSE,
                   return_intermediate = TRUE)
est_tree <- fit$tree
fit$variable_types
round(fit$distance, 3)
edges_from_ltree(est_tree)

## ----fig.width = 7, fig.height = 5--------------------------------------------
view_ltree(est_tree, "plot")

## ----eval = FALSE-------------------------------------------------------------
# view_ltree(est_tree)
# save_ltree(est_tree, file.path(tempdir(), "est_tree"))  # PNG, EPS and HTML

## -----------------------------------------------------------------------------
tree_gaussian <- learn_sltree(gaussian_x, model = "gaussian")
h_binary <- rbinom(n, 1, 0.5)
discrete_x <- as.data.frame(replicate(4, rbinom(n, 1, 0.1 + 0.8 * h_binary)))
names(discrete_x) <- paste0("V", 1:4)
tree_discrete <- learn_sltree(discrete_x, model = "discrete")
mixed_x <- data.frame(
  V1 = rnorm(n, mean = 2 * h_binary - 1, sd = 0.6),
  V2 = rnorm(n, mean = 2 * h_binary - 1, sd = 0.6),
  V3 = discrete_x$V3,
  V4 = discrete_x$V4
)
tree_mixed <- learn_sltree(mixed_x, model = "mixed")
c(stm = is_ltree(est_tree), gaussian = is_ltree(tree_gaussian),
  discrete = is_ltree(tree_discrete), mixed = is_ltree(tree_mixed))

## -----------------------------------------------------------------------------
variable_types <- infer_variable_types(x)
correlation <- latent_correlation(x, variable_types)
distance <- -log(abs(correlation)); diag(distance) <- 0
tree_from_steps <- sls(distance, tau1 = 3, tau2 = 5,
                       node_types = variable_types, node_names = colnames(x))

## -----------------------------------------------------------------------------
robinson_foulds_distance(true_tree, est_tree)   # 0 if fully recovered

