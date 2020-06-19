library(ape)
library(pvclust)
library(tidytree)
library(treeio)
# using the example from pvclust
example(pvclust) 

## 1st try with pvclust example----
### tranforming the pvclust tree object to tibble 
phylo_res <- ape::as.phylo(result$hclust) %>% 
  tidytree::as_tibble() 
### tranforming the pvclust bootstraps values to tibble with key column:"node"
tree_boots <- (round(result$edges[, c("si","au", "bp")],2)*100) %>% 
  as_tibble() %>%
  mutate(node= seq(max(phylo_res$node), treeio::Ntip(ape::as.phylo(result$hclust))+1))
### resulting tree after tibble transformation to tree_data object
new_tree <- phylo_res %>% 
full_join(tree_boots) %>% 
as.treedata()
### checking resulting plot to if it is the same as the one resulted from example(pvclust) 
new_tree %>% 
  ggtree() + 
  geom_tiplab() +
  geom_label(aes(label=bp), fill="green")

## 2nd try with pvclust example and slightly different function of ape::as.phylo----
### ape::as.phylo.hclust function downloaded from https://stackoverflow.com/questions/22749634/how-to-append-bootstrapped-values-of-clusters-tree-nodes-in-newick-format-in
### and edited for the purpose of as_tree_data.pvclust 
### modification at the commented lines
as.phylo.hclust_node <- function(x){
   N <- dim(x$merge)[1]
   edge <- matrix(0L, 2 * N, 2)
   edge.length <- numeric(2 * N)
   node <- integer(N)
   node[N] <- N + 2L
   cur.nod <- N + 3L
   j <- 1L
   for (i in N:1) {
     edge[j:(j + 1), 1] <- node[i]
     for (l in 1:2) {
       k <- j + l - 1L
       y <- x$merge[i, l]
       if (y > 0) {
         edge[k, 2] <- node[y] <- cur.nod
         cur.nod <- cur.nod + 1L
         edge.length[k] <- x$height[i] - x$height[y]
       }
       else {
         edge[k, 2] <- -y
         edge.length[k] <- x$height[i]
       }
     }
     j <- j + 2L
   }
   if (is.null(x$labels)) 
     x$labels <- as.character(1:(N + 1))  
 node.lab <- order(node) # here we keep the order for the edges
 obj <- list(edge = edge, edge.length = edge.length/2, 
   tip.label = x$labels, 
   Nnode = N, 
   node.label = paste(node.lab,"_edge",sep = "")) # export it to the final object 
 class(obj) <- "phylo"
 reorder(obj)
 }
### tranforming the pvclust tree object to tibble 
phylo_res <- as.phylo.hclust_node(result$hclust) %>% 
  tidytree::as_tibble() 
### tranforming the pvclust bootstraps values to tibble with key column:"node"
tree_boots <- (round(result$edges[, c("si","au", "bp")],2)*100) %>% 
        as_tibble() %>%
        mutate(label = paste(seq(1 , treeio::Ntip(ape::as.phylo(result$hclust))-1),"_edge",sep = ""))
### resulting tree after tibble transformation to tree_data object
new_tree <- phylo_res %>% 
full_join(tree_boots) %>% 
as.treedata()
### checking resulting plot to if it is the same as the one resulted from example(pvclust) 
new_tree %>% 
  ggtree() + 
  geom_tiplab() +
  geom_label(aes(label=bp), fill="green")

## Function as_tree_data.pvclust ----
pvclust.as_tree_data <- function(pvclust_obj){
  phylo_res <- as.phylo.hclust_node(pvclust_obj$hclust) %>% 
    tidytree::as_tibble() 
  new_tree <- phylo_res %>% 
    full_join(
      (round(pvclust_obj$edges[, c("si","au", "bp")],2)*100) %>% 
        as_tibble() %>%
        mutate(label = paste(seq(1,treeio::Ntip(ape::as.phylo(pvclust_obj$hclust))-1),"_edge",sep = ""))
    ) %>% 
    as.treedata()
  new_tree
}

pvclust_tidytree <- as_tree_data.pvclust(result)
write.beast(pvclust_tidy_tree)
