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

## Function as_tree_data.pvclust ----
as_tree_data.pvclust <- function(pvclust_obj){
  phylo_res <- ape::as.phylo(pvclust_obj$hclust) %>% 
    tidytree::as_tibble() 
  new_tree <- phylo_res %>% 
    full_join(
      (round(pvclust_obj$edges[, c("si","au", "bp")],2)*100) %>% 
        as_tibble() %>%
        mutate(node= seq(max(phylo_res$node), 
          treeio::Ntip(ape::as.phylo(pvclust_obj$hclust))+1))) %>% 
    as.treedata()
  new_tree
}

pvclust_tidytree <- as_tree_data.pvclust(result)
write.beast(pvclust_tidy_tree)
