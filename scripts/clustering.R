# =============================================================================
# title:        scripts/clustering.R
# summary:      script to cluster the data with hierarchical clustering.
#               leaves choice of linkage method and # clusters to user.
# dependencies: cluster
# =============================================================================

# clear environment
rm(list = ls(all.names = TRUE))

# dependencies
require(cluster)

# load site_df
load('rdata/site_df.RData')

# Create a Jaccard proximity matrix, `site_prox_mat`
#####
# Matching coefficient, SMC
#SMC <- function(A, B){(sum((A==B)))/(sum(A!=B)+sum((A==B)))}

# Jaccard similarity, J
J <- function(A, B){(sum((A==B)&(A==1)))/(sum(A!=B)+sum((A==B)&(A==1)))}

jaccard_dist_mat <- function(df)
{
  # input: df, a data frame of binary attributes, with rows as points
  # output: a jaccard distance proximity matrix
  
  # initialize proximity matrix
  proximity_matrix <- matrix(NA, nrow = nrow(df), ncol = nrow(df))
  
  # find all pairwise combinations
  pairs <- combn(1:nrow(df), 2)
  
  # compute jaccard distance (1 - Jaccard coefficient) for each pair
  dists <- apply(pairs,2,function(x){1-(J(df[x[1],], df[x[2],]))}) # using 1 - Jaccard for dissimilarity
  
  # arrange dists on a matrix
  for(i in 1:ncol(pairs)){
    proximity_matrix[pairs[,i][1], pairs[,i][2]] <- dists[i]
    proximity_matrix[pairs[,i][2], pairs[,i][1]] <- dists[i]
  }
  
  # make the diagonal 0
  diag(proximity_matrix) <- rep(0, nrow(proximity_matrix))
  
  # return proximity matrix
  return(proximity_matrix)
}


site_prox_mat <- jaccard_dist_mat(site_df)
rownames(site_prox_mat) <- rownames(site_df)
#####

# cluster the sites using hclust(). create and save assignments.
#####

hclust_average <- hclust(as.dist(site_prox_mat), method = "average")
hclust_single <- hclust(as.dist(site_prox_mat), method = "single")
hclust_complete <- hclust(as.dist(site_prox_mat), method = "complete")

clust_sil <- function(hc, k, proximity_matrix){
  sil <- cluster::silhouette(cutree(hc, k), dmatrix = proximity_matrix)
  return(sil)
}


output_clusters <- function(hclust_obj, k){
  cluster_neighbor <- clust_sil(hclust_obj, k, site_prox_mat)
  
  clustered_sites <- data.frame(
    
    "site" = rownames(site_df),
    "cluster" = cluster_neighbor[,1],
    "neighbor" = cluster_neighbor[,2],
    
    stringsAsFactors = F)
  
  cluster_size <- as.data.frame(table(clustered_sites$cluster))
  colnames(cluster_size) <- c("cluster", "size")
  
  clustered_sites$"cluster_size" = cluster_size[cluster_neighbor[,1],2]
  clustered_sites$"neighbor_size" = cluster_size[cluster_neighbor[,2],2]
  clustered_sites$"total_size" = clustered_sites$cluster_size+clustered_sites$neighbor_size
  
  return(clustered_sites)
}


assignments <- list(
  "average" = list(),
  "single" = list(),
  "complete" = list()
)

for (linkage in c("average", "single", "complete")){
  choose_obj <- function(linkage){
    switch(linkage,
           "average" = hclust_average,
           "single" = hclust_single,
           "complete" = hclust_complete)
  }
  for(k in 2:24){print(paste("linkage", linkage, "k=", k))
    assignments[[linkage]][[paste0("k_", k)]] = output_clusters(choose_obj(linkage), k) 
  }
}


#save(assignments,
#     file = "rdata/cluster-assignments.RData")
#####


