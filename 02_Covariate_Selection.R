################# Covariate Selection by Site ######################

library(here) #file input/output
library(vegan) #beta diversity and permanova function (adonis)

##### Data Input #####
df_mothers <- read.csv(here("Data", "Correct_Data_2022_02_19.csv"))

#Which columns are OTUs:
otus <- readRDS(here("Data", "otu_list_2022_02_19.rds"))

#Which columns are Covariates:
covars <- colnames(df_mothers)[c(194:197, 200:214, 217:310, 326:374)] #162

#need Batch and Cluster covariates! Make sure that no batch or cluster has less than 10 observations:
table(df_mothers$Batch) #Batch 11 is too small, combine it with batch 15:
df_mothers$BatchAdj <- ifelse(df_mothers$Batch == "WFSTL11", "WFSTL15", df_mothers$Batch)

table(df_mothers$Cluster.x) #Combine: 619 with 627, 815 with 826, and 911 915 916 921 924 931 932 and 934
df_mothers$ClusterAdj <- ifelse(df_mothers$Cluster.x == "619", "627",
                                ifelse(df_mothers$Cluster.x == "815", "826",
                                       ifelse(df_mothers$Cluster.x %in% c("911", "915", "916", "921", "924", "931", "932"), "934", df_mothers$Cluster.x)))

table(df_mothers$ClusterAdj)
table(df_mothers$BatchAdj)

covars <- c(covars, "ClusterAdj", "BatchAdj")

#Make standardized reads column as well:
df_mothers$std_reads <- scale(df_mothers$Reads)

write.csv(df_mothers, file = here("Data", "Correct_Data_2022_02_19.csv"), row.names = F)

#Which columns are covariates that MUST be included:
covars_include <- c("Supp_Status", "std_reads")


##### Covariate Potential #####
#by timepoint and site, test each potential covariates association with betadiversity (BC or Jaccard), store covariate
# a 'keeps' list by site if association of covariate and beta diversity is less than 0.1

sites <- unique(df_mothers$SiteLocation.x)
times <- unique(df_mothers$Timepoint)

covs_keeps_bc <- covs_keeps_jac <- list(DRC = NULL,
                   Guatemala = NULL,
                   India = NULL,
                   Pakistan = NULL)

df_mothers$Supp_Status <- as.factor(df_mothers$Supp_Status)

set.seed(43221) #Ensures reproducible results which is needed when working with permanova!

for(i in sites){
  for(j in times){
    
    if(j == "34Weeks" & i == "India"){
      next
    }
    
    df_tmp <- df_mothers[df_mothers$SiteLocation.x == i & df_mothers$Timepoint == j,]
    
    for(k in 1:length(covars)){
      if(k %in% c(1:2,9:16,21:22,26:27,30:33,36,74:83,145:161,163:164)){
        df_tmp[,covars[k]] <- as.factor(df_tmp[,covars[k]])
      }
      
      if(length(levels(df_tmp[,covars[k]])) ==1){
        next
      }
      
      if(sum(complete.cases(df_mothers[df_mothers$SiteLocation.x == i,covars[k]])) != length(df_mothers[df_mothers$SiteLocation.x == i,covars[k]])){
        next
      }
      
      
      cts <- df_tmp[,otus]
      perm_num <- 1000
      
      
      tryCatch(
      {
        tmp_mod_bc <- adonis2(as.formula(paste("cts ~ ",covars[k],"+ Supp_Status + Reads", sep = "")), 
                          data = df_tmp, permutations = perm_num, method = "bray")
    
        tmp_mod_jac <- adonis2(as.formula(paste("cts ~ ",covars[k],"+ Supp_Status + Reads", sep = "")), 
                           data = df_tmp, permutations = perm_num, method = "jaccard")
    
        if(tmp_mod_bc$`Pr(>F)`[1] < 0.01){
          covs_keeps_bc[[i]] <- c(covs_keeps_bc[[i]], covars[k])
        }
    
        if(tmp_mod_jac$`Pr(>F)`[1] < 0.01){
          covs_keeps_jac[[i]] <- c(covs_keeps_jac[[i]], covars[k])
        }
      })
    }
  }
}


####Identify groups of covariates that are highly correlated across all samples and timepoints
covar_keep_overall_bc <- c(unique(unlist(covs_keeps_bc)))
covar_keep_overall_jac <- c(unique(unlist(covs_keeps_jac)))

cor_values_bc <- cor_vals_jac <- list()
df_mothers_num <- apply(df_mothers[,covars], 2, function(x) as.numeric(x)) 
cor_mat_bc <- as.matrix(cor(df_mothers_num[,covar_keep_overall_bc], use = "pairwise.complete.obs"))
cor_mat_jac <- as.matrix(cor(df_mothers_num[,covar_keep_overall_jac], use = "pairwise.complete.obs"))

for(i in 1:length(covar_keep_overall_bc)){
  cor_values_bc[[covar_keep_overall_bc[i]]] <- names(which(abs(cor_mat_bc[covar_keep_overall_bc[i],])>0.7))
}

for(i in 1:length(covar_keep_overall_jac)){
  cor_vals_jac[[covar_keep_overall_jac[i]]] <- names(which(abs(cor_mat_jac[covar_keep_overall_jac[i],])>0.7))
}


##If overlap, which OTUs to keep?
#From Jaccard and Bray, we will keep SES_Common and electricity

covs_keeps_jac <- lapply(covs_keeps_jac, function(x) unique(ifelse(x == "electricity_01", "SES_common", x)))
covs_keeps_bc <- lapply(covs_keeps_bc, function(x) unique(ifelse(x == "electricity_01" , "SES_common", x)))

#From Jaccard and Bray, we will keep F1AA16 as a numeric and remove medu
covs_keeps_jac <- lapply(covs_keeps_jac, function(x) unique(ifelse(x == "medu", "F1AA16", x)))
covs_keeps_bc <- lapply(covs_keeps_bc, function(x) unique(ifelse(x == "medu" , "F1AA16", x)))

#ClusterAdj is going to be too many parameters for our model... we could include it as a randome effect that subject is nested within?
# We will just remove the variable for our analysis!
covs_keeps_bc2 <- covs_keeps_bc
covs_keeps_jac2 <- covs_keeps_jac
for(i in sites){
  covs_keeps_bc2[[i]] <- dplyr::setdiff(covs_keeps_bc[[i]],c("ClusterAdj"))
  covs_keeps_jac2[[i]] <- dplyr::setdiff(covs_keeps_jac[[i]],c("ClusterAdj"))
}

#Lists are almost identical!

#Remove any covariates that do not change (subset and view all matrices by site and selected otus. Remove if constant):
#Nothing to remove for India, DRC, Guatemala, or Pakistan


#Now add all covariates that need to be adjusted for in our models (Reads). We will use standardized reads to keep estimates on similar scales:
covs_keeps_jac2 <- lapply(covs_keeps_jac2,function(x) c(x,covars_include))
covs_keeps_bc2 <- lapply(covs_keeps_bc2,function(x) c(x,covars_include))

write.csv(df_mothers, file = here("Data", "Correct_Data_2022_02_19.csv"), row.names = F)
##Save covariate lists for future work!
saveRDS(covs_keeps_jac2, file=here("Data", "Covariate", "jaccard_covariates.rds"))
saveRDS(covs_keeps_bc2, file=here("Data", "Covariate", "brayCurtis_covariates.rds"))
