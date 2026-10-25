.sltree_model_cache <- new.env(parent = emptyenv())

#' Preset latent tree models
#'
#' Six fixed latent tree structures. Each tree is a three-column data frame
#' \code{nodes}/\code{type}/\code{parent} (\code{""} is the root).
#'
#' @param model_id An integer from 1 to 6, or \code{NULL} to return all models.
#' @return A three-column data frame, or a list named \code{model1} through
#'   \code{model6} when \code{model_id = NULL}.
#' @seealso \code{\link{simulate_ltree_data}}, \code{\link{generate_latent_tree}}
#' @examples
#' m1 <- latent_tree_models(1)
#' all_models <- latent_tree_models()
#' @export
latent_tree_models <- function(model_id = NULL) {
  if (!is.null(model_id) &&
      (length(model_id) != 1L || !(model_id %in% 1:6))) {
    stop("`model_id` must be a single integer from 1 to 6.")
  }
  if (exists("models", envir = .sltree_model_cache, inherits = FALSE)) {
    models <- get("models", envir = .sltree_model_cache, inherits = FALSE)
    if (is.null(model_id)) return(models)
    return(models[[paste0("model", model_id)]])
  }

# Models 1-2
# Model 1
node_ids <- c(paste("V", 1:16, sep = ""), paste("H", 1:8, sep = ""))
parent <- c(rep(paste("H", 1:8, sep = ""), each = 2),
                    "","H1","H1","H1","H2","H3","H3","H4")
type <- rep(c("observed_binary", "observed_Gaussian", "observed_binary",
                      "observed_Gaussian", "observed_binary", "latent_Gaussian"),
                    c(2, 4, 2, 4, 4, 8))
model1 <- data.frame("nodes" = node_ids, type, parent, stringsAsFactors = FALSE)
row.names(model1) <- node_ids

# Model 2
node_ids <- c(paste("V", 1:24, sep = ""), paste("H", 1:8, sep = ""))
parent <- c(rep(paste("H", 1:8, sep = ""), each = 3),
                    "","H1","H1","H1","H2","H3","H3","H4")
type <- rep(c("observed_binary", "observed_Gaussian", "observed_binary",
                      "observed_Gaussian", "observed_binary", "latent_Gaussian"),
                    c(3, 6, 3, 6, 6, 8))
model2 <- data.frame("nodes" = node_ids, type, parent, stringsAsFactors = FALSE)
row.names(model2) <- node_ids

# Models 3-5
# Model 3: Gaussian and binary observations
node_ids <- c(paste("V", 1:14, sep = ""), paste("H", 1:5, sep = ""))
parent <- c("","V1","H1","H1","H2","V2","H3","H3","V5","V5","H4","H4","H5","H5",
            "V1","V1","H1","H2","V2")
type <- rep(c("observed_Gaussian", "observed_binary", "observed_Gaussian",
              "observed_binary", "latent_Gaussian"),
            c(1, 5, 6, 2, 5))
model3 <- data.frame("nodes" = node_ids, type, parent, stringsAsFactors = FALSE)
row.names(model3) <- node_ids

# Model 4: Gaussian and ordinal observations
node_ids <- c(paste("V", 1:14, sep = ""), paste("H", 1:5, sep = ""))
parent <- c("","V1","H1","H1","H2","V2","H3","H3","V5","V5","H4","H4","H5","H5",
            "V1","V1","H1","H2","V2")
type <- rep(c("observed_Gaussian", "observed_ordinal", "observed_Gaussian",
              "observed_ordinal", "latent_Gaussian"),
            c(1, 5, 6, 2, 5))
model4 <- data.frame("nodes" = node_ids, type, parent, stringsAsFactors = FALSE)
row.names(model4) <- node_ids

# Model 5: non-Gaussian continuous and binary observations
node_ids <- c(paste("V", 1:14, sep = ""), paste("H", 1:5, sep = ""))
parent <- c("","V1","H1","H1","H2","V2","H3","H3","V5","V5","H4","H4","H5","H5",
            "V1","V1","H1","H2","V2")
type <- rep(c("observed_continuous", "observed_binary", "observed_continuous",
              "observed_binary", "latent_Gaussian"),
            c(1, 5, 6, 2, 5))
model5 <- data.frame("nodes" = node_ids, type, parent, stringsAsFactors = FALSE)
row.names(model5) <- node_ids

# Model 6: Gaussian and binary observations; 0 latent nodes
node_ids <- paste("V", 1:500, sep = "")
parent <- rep(c("","V1", paste0("V",2:10)), c(1, 49, rep(50, 9)))
type <- rep(c("observed_Gaussian", "observed_binary"), c(250, 250))
model6 <- data.frame("nodes" = node_ids, type, parent, stringsAsFactors = FALSE)
row.names(model6) <- node_ids

  models <- lapply(mget(paste0("model", 1:6)), function(model) {
    names(model) <- c("nodes", "type", "parent")
    model
  })
  assign("models", models, envir = .sltree_model_cache)
  if (is.null(model_id)) return(models)
  models[[paste0("model", model_id)]]
}
