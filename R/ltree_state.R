.LTREE_FIELDS <- c(
  "observed_nodes", "latent_nodes", "internal_nodes", "node_metadata",
  "parent_by_node", "children_by_node",
  "observed_descendants", "active_nodes", "discovered_nodes", "distance",
  "node_types", "node_names"
)


# Empty ltree object.
.new_ltree <- function() {
  structure(
    stats::setNames(vector(mode = "list", length = length(.LTREE_FIELDS)),
                    .LTREE_FIELDS),
    class = "ltree"
  )
}
