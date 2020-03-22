# =============================================================================
# title:        select-k.R
# summary:      selecting the 'optimal' clustering configuration
# dependencies: tidyverse
# =============================================================================

# clear environment
rm(list = ls(all.names = TRUE))

# dependencies
require(tidyverse)

# source `clustering.R`
source('scripts/clustering.R')

# load cluster assignments
load("rdata/cluster_assignments.RData")

# average within-cluster pairwise dissimilarity
#####
compute_W_k <- function(k, linkage){
  # extract cluster assignments
  cluster <- assignments[[linkage]][[paste0("k_", k)]]$cluster
  
  # create a data frame with cluster
  df <- as.data.frame(cbind(site_df, cluster))
  
  # function to divide by 2
  half <- function(x){x/2}
  
  w <- rep(NA, k)
  for(C in 1:k){
    proximity <- df %>% filter(cluster==C) %>% select(-cluster)
    
    # sanity check: singleton
    if(nrow(proximity)==1){
      w[C] <- 0 # singletons have 0 distance
    }
    else{
      # compute w
      w[C] <- proximity %>%
        jaccard_dist_mat %>%
        as.vector %>%
        mean %>%
        half
    }
  }
  return(mean(w))
}
#####

# table of cluster sizes per k and linkage
#####
cluster_size_table <- function(k, linkage)
{
  # extract cluster assignments
  cl <- assignments[[linkage]][[paste0("k_", k)]]$cluster
  
  # make a data frame of cluster sizes
  cluster_sizes <- cl %>% table %>% data.frame
  names(cluster_sizes) <- c("cluster", "size")
  
  return(cluster_sizes)
}
#####

# penalized W(k) for each linkage method
#####
penalized_W_k <- function(linkage, k_range) #k-range is a vector
{
  n_sing <- list()
  W <- list()
  
  for(k in as.character(k_range))
  {
    f_table <- cluster_size_table(k, linkage)
    n_sing[[k]] <- sum(f_table$size == 1)
    W[[k]] <- compute_W_k(k, linkage)
  }
  
  output <- list("n_sing" = unlist(n_sing),
                 "W_k" = unlist(W))
  
  # penalized W_k can be:
  # n_sing + W_k
  # n_sing/k + W_k
  
  return(output)
}
#####


# average cluster size for each linkage and k
#####
all_W_k <- list(
  "average" = list(),
  "complete" = list(),
  "single" = list()
)


k_range <- 2:22

for(linkage in c("average",
                 "complete",
                 "single"))
{
  all_W_k[[linkage]] = penalized_W_k(linkage = linkage, 
                                     k_range = 2:22)
}



all_cl_size <- list(
  "average" = list(),
  "complete" = list(),
  "single" = list()
)

for(k in k_range)
{
  for(linkage in c("average",
                   "complete",
                   "single"))
  {
    all_cl_size[[linkage]][[as.character(k)]] = median(cluster_size_table(k, linkage)$size)
  }
}

sd_cl_size <- list(
  "average" = list(),
  "complete" = list(),
  "single" = list()
)

for(k in k_range)
{
  for(linkage in c("average",
                   "complete",
                   "single"))
  {
    sd_cl_size[[linkage]][[as.character(k)]] = sd(cluster_size_table(k, linkage)$size)
  }
}
#####

# what is the cluster size of target?
#####
#target_cl_size <- list(
#  "average" = list(),
#  "complete" = list(),
#  "single" = list()
#)

#for(k in k_range)
#{
#  for(linkage in c("average",
#                   "complete",
#                   "single"))
#  {
#    cl_assignment <- assignments[[linkage]][[paste0("k_",k)]]
#    target_cl_size[[linkage]][[as.character(k)]] <- cl_assignment[cl_assignment$site == "target",]$cluster_size
#  }
#}
#####

# clustering critereon
#####
k_range = 2:22

cluster_critereon <- data.frame(
  "linkage" = rep(c("average",
                    "complete",
                    "single"),
                  each = length(k_range)),
  "W" = unlist(lapply(all_W_k,
                        function(x){x$W_k})),
  "n_sing" = unlist(lapply(all_W_k,
                           function(x){x$n_sing})),
  "k" = rep(k_range, 3),
  "med_cl_size" = unlist(all_cl_size),
  "sd_cl_size" = unlist(sd_cl_size)#,
  #"target_cl_size" = unlist(target_cl_size)
  )

rownames(cluster_critereon) <- NULL

# best 10
cluster_critereon %>%
  arrange(W+(n_sing/k)) %>%
  head(10)

# we select 'complete' with k = 7
split(assignments$complete$k_7$site,
      assignments$complete$k_7$cluster)

