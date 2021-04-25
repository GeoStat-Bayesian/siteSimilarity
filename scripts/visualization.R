# =============================================================================
# title:        scripts/visualization.R
# summary:      script to visualize clustering. we use complete linkage, k=7
# dependencies: ggplot2
# =============================================================================
# clear environment
rm(list = ls(all.names = TRUE))

# dependencies
require(ggplot2)
require(tidyverse)
# source
source('scripts/clustering.R')

# extract cluster assignments
k = 7
cluster_output <- output_clusters(hclust_complete, k)

# plot the dendrogram
par(cex=0.75)
plot(as.dendrogram(hclust_complete,
                   hang = 0.05),
     type = 'triangle',
     center = TRUE,
     main = 'Cluster Dendrogram of Hydrogelogical Sites')


# create a data frame with cluster
df <- as.data.frame(cbind(site_df, cluster_output$cluster))
colnames(df)[51] <- 'cluster'

# function to divide by 2
half <- function(x){x/2}

cluster_info <- matrix(NA,nrow=k,ncol=50)
colnames(cluster_info) <- colnames(df)[1:50]
for(C in 1:k)
{
  proximity <- df %>% filter(cluster==C) %>% select(-cluster)
  cluster_info[C,] = round(colMeans(proximity),2)
}

cluster_sizes <- as.data.frame(table(df$cluster))

cols <- c("key_rt","rt_id_parent","id_env", "env_id_parent","key_Fract")

feature <- c(rep("fractology", sum(grepl(cols[5], colnames(cluster_info)))),
             rep("env.type", sum(grepl(cols[3], colnames(cluster_info)))),
             rep("env.parent", sum(grepl(cols[4], colnames(cluster_info)))),
             rep("rock.type", sum(grepl(cols[1], colnames(cluster_info)))),
             rep("rt.parent", sum(grepl(cols[2], colnames(cluster_info)))))

feature <- factor(feature)



cl_k_df <- data.frame(
  "name" = rep(colnames(cluster_info),k),
  "feature" = rep(feature, k),
  "k" = rep(1:k, each=50),
  "freq" = as.vector(t(cluster_info)),
  "size" = rep(cluster_sizes$Freq, each=50)
)

cluster_purity <- ggplot(cl_k_df,
                         aes(
                           x = name,
                           y = freq,
                           fill = feature
                         )) + 
  geom_bar(stat = "identity") + 
  facet_grid(k~.) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  #theme_bw() +
  ggtitle("Categorical Features per Cluster")

gg_cluster_sizes <- ggplot(cluster_sizes,
                           aes(x = Var1,
                               y = Freq,
                               fill = Freq)) + geom_bar(stat="identity") +
  geom_hline(yintercept = median(cluster_sizes$Freq),
             color = "red") +
  ggtitle("Cluster Sizes")


cluster_purity
gg_cluster_sizes




