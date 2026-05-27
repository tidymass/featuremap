library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)

set.seed(1)

n_feat   <- 8
n_annot  <- 15
n_hidden <- 10

feat_ids   <- paste0("F", 1:n_feat)
annot_ids  <- paste0("A", 1:n_annot)
hidden_ids <- paste0("H", 1:n_hidden)

nodes <- tibble(
  name = c(feat_ids, annot_ids, hidden_ids),
  type = c(rep("feature", n_feat),
           rep("annotated", n_annot),
           rep("hidden", n_hidden))
)

edge_list <- list()
for (f in feat_ids)
  edge_list[[length(edge_list)+1]] <- c(f, sample(annot_ids, 1))
for (h in hidden_ids) {
  k <- sample(3:5, 1)
  for (a in sample(annot_ids, k))
    edge_list[[length(edge_list)+1]] <- c(a, h)
}
for (i in 1:20)
  edge_list[[length(edge_list)+1]] <- sample(annot_ids, 2)

edges <- do.call(rbind, edge_list) %>%
  as.data.frame() %>% setNames(c("from","to")) %>% distinct()

g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

# --- Layout: tighter, less radial separation ---
g_ann <- induced_subgraph(as.igraph(g), annot_ids)
lay_ann <- layout_with_fr(g_ann)
rownames(lay_ann) <- annot_ids
lay_ann <- lay_ann / max(abs(lay_ann)) * 0.50        # annotated core (tight)

# hidden: interspersed close to annotated, not in a separate ring
ah <- runif(n_hidden, 0, 2*pi)
rh <- runif(n_hidden, 0.30, 0.55)                    # overlaps annotated zone -> mixed in
lay_hid <- cbind(rh*cos(ah), rh*sin(ah)); rownames(lay_hid) <- hidden_ids

# features: scattered just outside the core, much closer now
af <- runif(n_feat, 0, 2*pi)
rf <- runif(n_feat, 0.50, 0.68)                      # was 1.0-1.35, pulled in
lay_feat <- cbind(rf*cos(af), rf*sin(af)); rownames(lay_feat) <- feat_ids

lay <- matrix(0, nrow(nodes), 2); rownames(lay) <- nodes$name
lay[annot_ids, ]  <- lay_ann
lay[hidden_ids, ] <- lay_hid
lay[feat_ids, ]   <- lay_feat

# --- Plot ---
p <- ggraph(g, layout = "manual", x = lay[nodes$name,1], y = lay[nodes$name,2]) +
  geom_edge_link(colour = "grey75", alpha = 0.55, width = 0.4) +
  geom_node_point(aes(colour = type, shape = type, size = type)) +
  scale_shape_manual(values = c(feature = 15, annotated = 16, hidden = 16)) +
  scale_colour_manual(values = c(feature   = "#ba3339",
                                 annotated = "#397fc3",
                                 hidden    = "#9e9e9e")) +
  scale_size_manual(values = c(feature = 4, annotated = 5.5, hidden = 4)) +
  theme_void() +
  theme(legend.position = "none",                    # no legend
        plot.background = element_rect(fill = "white", colour = NA))

pdf("./Figures/Figure1/pcsf_eg2.pdf", width = 2, height = 2)
p
dev.off()

# plot final pcsf sub-network --------------------------------------------------

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)
set.seed(101)
n_feat <- 12; n_annot <- 18; n_hidden <- 8
feat_ids   <- paste0("F", 1:n_feat)
annot_ids  <- paste0("A", 1:n_annot)
hidden_ids <- paste0("H", 1:n_hidden)
nodes <- tibble(
  name = c(feat_ids, annot_ids, hidden_ids),
  type = c(rep("feature", n_feat), rep("annotated", n_annot), rep("hidden", n_hidden))
)
edge_list <- list()
for (f in feat_ids) edge_list[[length(edge_list)+1]] <- c(f, sample(annot_ids, 1))
for (h in hidden_ids) { k <- sample(3:5,1); for (a in sample(annot_ids,k)) edge_list[[length(edge_list)+1]] <- c(a,h) }
for (i in 1:20) edge_list[[length(edge_list)+1]] <- sample(annot_ids, 2)
edges <- do.call(rbind, edge_list) %>% as.data.frame() %>% setNames(c("from","to")) %>% distinct()

set.seed(21)
edges$weight <- runif(nrow(edges), 0.3, 2.2)
base <- c(feature = 4, annotated = 5.5, hidden = 4)
mult <- numeric(nrow(nodes))
mult[nodes$type=="feature"]   <- runif(sum(nodes$type=="feature"),   0.6, 1.8)
mult[nodes$type=="annotated"] <- runif(sum(nodes$type=="annotated"), 0.7, 1.4)
mult[nodes$type=="hidden"]    <- runif(sum(nodes$type=="hidden"),    0.7, 1.4)
nodes$node_size <- base[nodes$type] * mult
g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

scatter_ring <- function(ids, r_lo, r_hi, seed) {
  set.seed(seed); k <- length(ids)
  ang <- seq(0, 2*pi, length.out = k + 1)[1:k] + runif(k, -0.25, 0.25)
  rad <- runif(k, r_lo, r_hi)
  m <- cbind(rad*cos(ang), rad*sin(ang)); rownames(m) <- ids; m
}
core_ids <- sample(c(annot_ids, hidden_ids))
lay_core <- scatter_ring(core_ids, 0.0, 0.62, seed = 7)
lay_feat <- scatter_ring(feat_ids, 0.72, 0.92, seed = 13)
lay <- matrix(0, nrow(nodes), 2); rownames(lay) <- nodes$name
lay[rownames(lay_core), ] <- lay_core
lay[rownames(lay_feat), ] <- lay_feat

# ---- Collision relaxation: push apart any points closer than min_dist ----
relax_overlap <- function(xy, min_dist = 0.16, iters = 200, step = 0.5) {
  for (it in 1:iters) {
    moved <- FALSE
    for (i in 1:(nrow(xy)-1)) for (j in (i+1):nrow(xy)) {
      d <- xy[j,] - xy[i,]
      dist <- sqrt(sum(d^2))
      if (dist < min_dist && dist > 1e-9) {
        push <- (min_dist - dist) / 2 * step
        u <- d / dist
        xy[i,] <- xy[i,] - u * push
        xy[j,] <- xy[j,] + u * push
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  xy
}
lay <- relax_overlap(lay[nodes$name, ], min_dist = 0.18, iters = 300, step = 0.6)

p <- ggraph(g, layout = "manual", x = lay[nodes$name,1], y = lay[nodes$name,2]) +
  geom_edge_link(aes(width = weight), colour = "grey75", alpha = 1) +
  geom_node_point(aes(colour = type, shape = type, size = node_size)) +
  scale_edge_width(range = c(0.2, 1.8), guide = "none") +
  scale_size(range = c(2, 7), guide = "none") +
  scale_shape_manual(values = c(feature = 15, annotated = 16, hidden = 16)) +
  scale_colour_manual(values = c(feature = "#ba3339", annotated = "#397fc3", hidden = "#9e9e9e")) +
  theme_void() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA))

pdf("./Figures/Figure1/pcsf_final.pdf", width = 4, height = 4)
print(p)
dev.off() 



