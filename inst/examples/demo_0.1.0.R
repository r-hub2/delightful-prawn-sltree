# sltree demo

library(sltree)

# 1. View a preset tree.
view_ltree(1)

# 2. Latent root H_0 with three latent children, each with four observed variables.
set.seed(2026)
n <- 12000
H_0 <- rnorm(n)
latent_child <- function(parent, rho) {
  rho * parent + sqrt(1 - rho^2) * rnorm(length(parent))
}
H_A <- latent_child(H_0, 0.78)
H_B <- latent_child(H_0, 0.72)
H_C <- latent_child(H_0, 0.68)

observed_child <- function(parent, rho = 0.82) {
  rho * parent + sqrt(1 - rho^2) * rnorm(length(parent))
}
A <- replicate(4, observed_child(H_A))
B <- replicate(4, observed_child(H_B))
C <- replicate(4, observed_child(H_C))

x <- data.frame(
  A_C1 = exp(A[, 1]),
  A_B  = as.integer(A[, 2] > 0),
  A_O  = findInterval(A[, 3], c(-0.70, 0.50)),
  A_C2 = A[, 4]^3,
  B_C1 = exp(B[, 1]),
  B_B  = as.integer(B[, 2] > 0.10),
  B_O  = findInterval(B[, 3], c(-0.60, 0.60)),
  B_C2 = B[, 4]^3,
  C_C1 = exp(C[, 1]),
  C_B  = as.integer(C[, 2] > -0.10),
  C_O  = findInterval(C[, 3], c(-0.80, 0.40)),
  C_C2 = C[, 4]^3
)

print(head(x))
print(infer_variable_types(x))

# 3. Learn the tree (model = "stm").
tree <- learn_sltree(x)
stopifnot(length(tree$latent_nodes) == 4L)
print(edges_from_ltree(tree))

# 4. Summary and plots.
view_ltree(tree, "summary")
view_ltree(tree, "plot")
view_ltree(tree)

# 5. Save PNG, EPS and HTML (HTML needs visNetwork).
save_ltree(tree, "recovered_tree")

# 6. Keep the intermediate correlation and distance matrices.
fit <- learn_sltree(x, return_intermediate = TRUE)
print(fit$model)
print(fit$allow_observed_internal)
print(fit$run_time)

# 7. Gaussian, discrete and parametric mixed models.
x_gaussian <- data.frame(
  G1 = H_0 + rnorm(n, sd = 0.60),
  G2 = 0.9 * H_0 + rnorm(n, sd = 0.65),
  G3 = 0.8 * H_0 + rnorm(n, sd = 0.70),
  G4 = 0.7 * H_0 + rnorm(n, sd = 0.75),
  G5 = 0.6 * H_0 + rnorm(n, sd = 0.80)
)
tree_gaussian <- learn_sltree(x_gaussian, model = "gaussian")

h <- rbinom(n, size = 1, prob = 0.5)
x_discrete <- data.frame(
  D1 = rbinom(n, 1, 0.15 + 0.70 * h),
  D2 = rbinom(n, 1, 0.20 + 0.60 * h),
  D3 = rbinom(n, 1, 0.25 + 0.50 * h),
  D4 = rbinom(n, 1, 0.10 + 0.75 * h),
  D5 = rbinom(n, 1, 0.30 + 0.45 * h)
)
tree_discrete <- learn_sltree(x_discrete, model = "discrete")

x_mixed <- data.frame(
  M_C1 = rnorm(n, mean = -1.0 + 2.0 * h, sd = 0.8),
  M_C2 = rnorm(n, mean = -0.8 + 1.6 * h, sd = 0.9),
  M_C3 = rnorm(n, mean = -0.6 + 1.2 * h, sd = 1.0),
  M_B1 = rbinom(n, 1, 0.15 + 0.70 * h),
  M_B2 = rbinom(n, 1, 0.20 + 0.60 * h)
)
tree_mixed <- learn_sltree(x_mixed, model = "mixed")

print(c(stm = is_ltree(tree), gaussian = is_ltree(tree_gaussian),
        discrete = is_ltree(tree_discrete), mixed = is_ltree(tree_mixed)))
