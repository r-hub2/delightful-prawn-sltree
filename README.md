# sltree

Structure learning for semi-parametric/latent tree models from mixed
continuous, binary and ordinal data.

Authors: Can Zhou, Shuai Wang, Xiaofei Wang, Jianhua Guo.

Latent correlations are estimated under a Gaussian copula by inverting the
bridge functions between Kendall's tau and the latent Pearson correlation.
They are converted to the additive information distance `-log|rho|`, and the
SLS algorithm then recovers the tree topology, including the number and
placement of latent nodes. The bridge inversion is written in C++
(`src/latent_cpp.cpp`).

The package also covers Gaussian, discrete and parametric mixed latent tree
models through the same interface.

## Installation

Install from CRAN with:

```r
install.packages("sltree")
```

CRAN binary packages for Windows and macOS do not require a local compiler.
Installing from source requires a C++ compiler and the R development tools
for your platform (Rtools on Windows or the Xcode command line tools on macOS):

```r
install.packages("sltree_0.1.0.tar.gz", repos = NULL, type = "source")
```

Imports: igraph, mvtnorm, pcaPP, Rcpp.

Suggests: Matrix (`nearPD`), pcalg and graph (data simulation; graph is on
Bioconductor, `BiocManager::install("graph")`), visNetwork (interactive
plots), testthat, knitr, rmarkdown. None of these is needed for
`infer_variable_types()`, `latent_correlation()`, `sls()` or
`robinson_foulds_distance()`.

## Example

```r
library(sltree)
model     <- latent_tree_models(1)
dat       <- simulate_ltree_data(2e4, model)      # needs pcalg and graph
est_tree  <- learn_sltree(dat$data, tau1 = 3, tau2 = 5)
true_tree <- ltree_from_edge(model[, -2], setNames(model[, 2], model[, 1]))
robinson_foulds_distance(true_tree, est_tree)     # 0 if fully recovered
```

See `vignette("sltree-intro", package = "sltree")` for a longer example, and
`system.file("examples", "demo_0.1.0.R", package = "sltree")` for a
standalone script.

## Main functions

- `learn_sltree()`: estimate the information distance and learn the tree.
- `latent_correlation()`, `information_distance()`, `sls()`: the individual
  steps of `learn_sltree()`.
- `latent_tree_models()`: six preset trees.
- `generate_latent_tree()`, `simulate_ltree_data()`: random trees and mixed
  data generated from a tree.
- `robinson_foulds_distance()`: compare two trees.
- `view_ltree()`, `save_ltree()`: summaries, plots and file output.

## Models

```r
tree_stm      <- learn_sltree(x, model = "stm")            # semi-parametric (default)
tree_gaussian <- learn_sltree(x_gaussian, model = "gaussian")
tree_discrete <- learn_sltree(x_discrete, model = "discrete")
tree_mixed    <- learn_sltree(x_mixed, model = "mixed")    # Zhou et al. (2020)
```

The `stm`, `gaussian` and `discrete` models allow observed variables to be
internal nodes; `mixed` restricts them to leaves. `allow_observed_internal`
overrides this and gives a warning if it conflicts with the model.

Only `stm` projects the correlation matrix with `nearPD` by default. Use
`use_near_pd = TRUE` or `FALSE` to change this.

## Plotting

```r
view_ltree(1)                           # interactive plot of preset tree 1
view_ltree(est_tree, "summary")
view_ltree(est_tree)                    # interactive (visNetwork)
view_ltree(est_tree, "plot")            # static
save_ltree(est_tree, "estimated_tree")  # PNG, EPS and HTML
```

Without visNetwork, `view_ltree()` draws the static plot. Both plots show a
node-type legend. Options follow igraph (`vertex.*`, `edge.*`) and base
graphics (`main`, `xlab`, ...); see `?plot.ltree`, `?view_ltree` and
`?save_ltree`.

| Option | Static (`"plot"`) | Interactive (default) |
|---|---|---|
| `direction` | `"LR"` (default), `"RL"`, `"UD"`, `"DU"` | `"UD"` (default), `"DU"`, `"LR"`, `"RL"` |
| `vertex.color`, `vertex.frame.color` | one colour, one per node, named by node, or named by type (`latent`, `continuous`, `discrete`) | same |
| `legend`, `legend.args` | `legend()` | `visLegend()` |
| `vertex.args` | `polygon()` | `visNodes()` |
| `edge.args` | `segments()` | `visEdges()` |
| `label.args`, `par.args` | `text()`, `par()` | not used |
| other arguments | `title()`: `main`, `sub`, `xlab`, `ylab`, ... | `visNetwork()`: `main`, `submain`, `footer`, ... |

`save_ltree()` also takes `plot_args` for the static images and `png_args`,
`eps_args` and `html_args`, passed to `png()`, `postscript()` and
`visSave()`. Other arguments go to the interactive plot.

```r
view_ltree(est_tree, "plot", main = "Recovered tree", direction = "UD",
           vertex.color = c(latent = "grey60"),
           edge.args = list(lty = 2),
           legend.args = list(x = "topright", title = "Node type"))
view_ltree(est_tree, legend = FALSE, vertex.args = list(font = list(size = 20)))
save_ltree(est_tree, "estimated_tree",
           plot_args = list(main = "Recovered tree"),
           png_args = list(res = 300))
```

## Acknowledgments
The research underlying this package was supported by the Key Program of 
the National Natural Science Foundation of China (Grant No. 12431009), 
the National Key Research and Development Program of China (Grant No. 
2020YFA0714100), the National Natural Science Foundation of China 
(Grant Nos. 12171076 and 12226003), and the Natural Science Foundation
of Jiangsu Province (Grant No. BK20230439).

## References

Zhou, C., Li, N., Wang, S., Wang, X. and Guo, J. (2026). Learning
semi-parametric tree models from mixed data. *Artificial Intelligence*, 353,
104499. <https://doi.org/10.1016/j.artint.2026.104499>

Zhou, C., Wang, X. and Guo, J. (2020). Learning mixed latent tree models.
*Journal of Machine Learning Research*, 21(249), 1-35.
<https://jmlr.org/papers/v21/20-365.html>
