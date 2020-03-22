# siteSimilarity

scripts
  - preprocessing.R:            data preprocessing
  - clustering.R:               clustering procedure
  - select-k.R                  selecting best clustering of data 
  - visualization.R             visualizing results of clustering
rdata
  - site_df.RData:              output of `preprocessing.R`
  - cluster_assignments.RData:  output of `clustering.R`
output


Workflow  
1. preprocessing.R --> site_df.RData
2. clustering.R --> cluster_assignments.RData
3. select-k.R
4. visualization.R
