# =============================================================================
# title:        scripts/preprocessing.R
# summary:      preprocesses data before clustering. output 'site_df'.
# dependencies: geostatDB, tidyverse
# =============================================================================

# dependencies
require(tidyverse)
require(geostatDB)

# load the wwhypda from geostatDB
wwhypda <- geostatDB::getData()

# clean the data: obtain a dataframe of sites and qualities in order to cluster

# for reference, make a data frame of rock types and parents
rocks <- select(wwhypda, contains("rt")) %>% distinct %>% arrange(rt_id_parent)

# function trim.trailing: remove the white space at the end of the name
trim.trailing <- function (x) sub("\\s+$", "", x)

sites <- wwhypda %>% select(id_Measure,val,id_smpl,
                            site_id, site_name, #ISO_code, country_name,
                            key_Fract, id_env, env_id_parent,
                            key_rt, rt_id_parent,
                            param_name) %>%  # this selects the columns we want
  filter(param_name == "hydraulic conductivity") %>% # this selects the sites with hyd. cond.
  mutate(param_name  = NULL) %>% # removes column param_name
  arrange(site_name) %>% # arranges by site name
  mutate(site_name = tolower(site_name)) %>% # makes sites lower case
  mutate(site_name = trim.trailing(site_name)) # removes extra spaces

# correct spelling of sites
site_names <- sites %>% select(site_name) %>% 
  #mutate(site_name = tolower(site_name)) %>% 
  distinct

site_names

# inspect by eye which are misspelled
misspelled <- c("helsby sandstone formation",
                "penruith sandstone",
                "permo-triassic sandstones ",
                "sbb-graulztunnel bei bern",
                "sherwood sandstone group",
                "the chichester chalk block, swanbourne south n?,??3",
                "the nortern great plains",
                "the northern great plains, emphasis on north dakota",
                "the northern great plaint")

correction <- c("helsby-ormskirk sandstone formation",
                "penrith sandstone",
                "permo-triassic sandstones",
                "sbb-grauholztunnel",
                "sherwood sandstone",
                "the chichester chalk block, swanbourne south",
                "the northern great plains",
                "the northern great plains",
                "the northern great plains")

correct_spelling <- function(x)
{
  if(sum(grepl(x, misspelled)) > 0){x <- correction[grepl(x, misspelled)==TRUE]}
  else{x <- x}
  return(x)
}

# apply this to a vector - sapply is returning a list, for some reason
correct_spelling_vector <- function(x){
  y <- rep(NA, length(x))
  for(j in 1:length(x))
  {y[j] <- correct_spelling(x[j])}
  return(y)}



# fix the names of sites
sites <- sites %>%  mutate(site_name = correct_spelling_vector(sites$site_name)) %>%
  arrange(site_name) %>%
  filter(site_name != "unknown") %>% filter(site_name != "na") %>% 
  mutate(site_name = trim.trailing(site_name))



# select data frame for clustering sites
clustering_sites <- sites %>% select(site_name,
                                     #ISO_code,
                                     key_Fract,
                                     id_env, env_id_parent,
                                     key_rt,rt_id_parent) %>% distinct

#clustering_sites # notice that some sites have multiple rock types and environments
# so we need to make binary matrices showing whether or not a site contains a value


# notice that some sites have multiple rock types and environments
# so we need to make binary matrices showing whether or not a site contains 
# a value

create_binary_predictors <- function(x)
{
  the_table <- table(clustering_sites$site_name,
                     clustering_sites[,x])
  colnames(the_table) <- paste0(x, "_", colnames(the_table))
  return(the_table)
}

# doing this for all 6 columns

mat_list <- list(list())
for(i in 1:length(colnames(clustering_sites)[2:6])){
  mat_list[[i]] <- create_binary_predictors(colnames(clustering_sites)[2:6][i])
}

num_rows <- length(unique(clustering_sites$site_name))

pred_mats <- as.data.frame(matrix(unlist(mat_list),
                                  nrow = num_rows,
                                  byrow = F))

colnames(pred_mats) <- unlist(lapply(mat_list, colnames))

pred_mats <- (pred_mats >= 1) + 0 # making it binary for containing, not counts
site_names <- clustering_sites %>% select(site_name) %>% distinct
rownames(pred_mats) <- site_names$site_name

site_df <- pred_mats

save(site_df,
     file = 'rdata/site_df.RData')
