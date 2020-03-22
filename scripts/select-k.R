# =============================================================================
# title:        select-k.R
# summary:      selecting the 'optimal' clustering configuration
# dependencies: tidyverse, ggplot2
# =============================================================================

# clear environment
rm(list = ls(all.names = TRUE))

# dependencies
require(tidyverse)
require(ggplot2)
require(exPrior)

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
#####

# visualization of selection procedure
#####

head(cluster_critereon)

# W(k), not penalized
W_k_plot <- ggplot(cluster_critereon,
                   aes(x = k, 
                       y = W, 
                       group = linkage, 
                       color = linkage)) + 
  geom_line() + 
  geom_point(size = 1.7) +
  #geom_vline(xintercept = 8, linetype = "dashed") + 
  labs(y = "W(k)") + 
  scale_color_manual(values = c("tomato", "blue", "orange")) + 
  scale_x_continuous(breaks = k_range) +
  #geom_rect(aes(xmin=7, xmax=9, ymin=-inf, ymax=Inf,
  #              size = 0,
  #              alpha = 0.7))+
  
  #geom_vline(xintercept = 7, linetype = "dashed") +
  #geom_vline(xintercept = 8, linetype = "dashed") +
  #geom_vline(xintercept = 9, linetype = "dashed") +
  #geom_rect(mapping=aes(xmin=6.8, xmax=11.2, ymin=0.05, ymax=0.2), 
  #          color="black", alpha=0.0) +
  
  ggtitle("W(k), not penalized") + theme_classic() #+
  #geom_vline(xintercept = 11) +
  #geom_vline(xintercept = 14) + 
  #geom_vline(xintercept = 4) + geom_vline(xintercept = 7)

W_k_plot

# number of singletons
n_sing_plot <- ggplot(cluster_critereon,
                      aes(x = k, 
                          y = n_sing, 
                          group = linkage, 
                          color = linkage)) + 
  geom_line() + 
  geom_point(size = 1.7) +
  #geom_vline(xintercept = 8, linetype = "dashed") + 
  labs(y = "n_sing") + 
  scale_color_manual(values = c("tomato", "blue", "orange")) + 
  scale_x_continuous(breaks = k_range) +
  #geom_rect(aes(xmin=7, xmax=9, ymin=-inf, ymax=Inf,
  #              size = 0,
  #              alpha = 0.7))+
  
  #geom_vline(xintercept = 7, linetype = "dashed") +
  #geom_vline(xintercept = 8, linetype = "dashed") +
  #geom_vline(xintercept = 9, linetype = "dashed") +
  #geom_rect(mapping=aes(xmin=6.8, xmax=11.2, ymin=0.05, ymax=0.2), 
  #          color="black", alpha=0.0) +
  
  ggtitle("Number of Singleton Clusters") + theme_classic()

n_sing_plot

# W(k) + n_sing
W_k_n_sing_plot <- ggplot(cluster_critereon,
                          aes(x = k, 
                              y = W+n_sing, 
                              group = linkage, 
                              color = linkage)) + 
  geom_line() + 
  geom_point(size = 1.7) +
  #geom_vline(xintercept = 8, linetype = "dashed") + 
  labs(y = "W(k) + n_sing") + 
  scale_color_manual(values = c("tomato", "blue", "orange")) + 
  scale_x_continuous(breaks = k_range) +
  #geom_rect(aes(xmin=7, xmax=9, ymin=-inf, ymax=Inf,
  #              size = 0,
  #              alpha = 0.7))+
  
  #geom_vline(xintercept = 7, linetype = "dashed") +
  #geom_vline(xintercept = 8, linetype = "dashed") +
  #geom_vline(xintercept = 9, linetype = "dashed") +
  #geom_rect(mapping=aes(xmin=6.8, xmax=11.2, ymin=0.05, ymax=0.2), 
  #          color="black", alpha=0.0) +
  
  ggtitle("W(k), penalized") + theme_classic()


# W(k) + (n_sing/k)
penalized_plot <- ggplot(cluster_critereon,
                         aes(x = k, 
                             #y = mean_cluster_size,
                             y = W+(n_sing/k), 
                             group = linkage, 
                             color = linkage)) + 
  geom_line() + 
  geom_point(size = 1.7) +
  #geom_vline(xintercept = 8, linetype = "dashed") + 
  labs(y = "W(k) + n_sing/k") + 
  scale_color_manual(values = c("tomato", "blue", "orange")) + 
  scale_x_continuous(breaks = k_range) +
  #geom_rect(aes(xmin=7, xmax=9, ymin=-inf, ymax=Inf,
  #              size = 0,
  #              alpha = 0.7))+
  
  #geom_vline(xintercept = 7, linetype = "dashed") +
  #geom_vline(xintercept = 8, linetype = "dashed") +
  #geom_vline(xintercept = 9, linetype = "dashed") +
  #geom_rect(mapping=aes(xmin=6.8, xmax=11.2, ymin=0.05, ymax=0.2), 
  #          color="black", alpha=0.0) +
  
  ggtitle("W(k), penalized") + theme_classic()

penalized_plot

# find the minimizer

cluster_critereon[,"penalized"] <- cluster_critereon$W + (cluster_critereon$n_sing/cluster_critereon$k)

cluster_critereon %>%
  arrange(penalized)

#####

