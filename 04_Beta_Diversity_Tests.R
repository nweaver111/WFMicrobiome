############### Beta Diversity Analysis #################

library(here)
library(vegan)

##### Set Seed #####
set.seed(1234)

##### Data Input #####

#Dataset:
df_mothers <- read.csv(here("Data", "Correct_Data_2022_02_19.csv"))
rownames(df_mothers) <- df_mothers$Lib

otus <- readRDS(here("Data", "otu_list_2022_02_19.rds"))

#Covariate Lists:
cov_list_bc <- readRDS(here("Data", "Covariate", "brayCurtis_covariates.rds"))
cov_list_jac <- readRDS(here("Data", "Covariate", "jaccard_covariates.rds"))


#make sure all covariates are correct type:
df_mothers$SES_common <- as.numeric(df_mothers$SES_common)
df_mothers$BatchAdj <- as.factor(df_mothers$BatchAdj)
df_mothers$F1AA16 <- as.numeric(df_mothers$F1AA16)
df_mothers$impr_water_01 <- as.factor(df_mothers$impr_water_01)
df_mothers$deliv_date <- as.numeric(df_mothers$deliv_date)
df_mothers$Supp_Status <- as.factor(df_mothers$Supp_Status)
df_mothers$std_reads <- as.numeric(df_mothers$std_reads)
df_mothers$s1_comp <- as.numeric(df_mothers$s1_comp)
df_mothers$impr_cook_fuel_01 <- as.factor(df_mothers$impr_cook_fuel_01)


##### Models #####
#set number of permutations:
perm = how(nperm = 1000)

#Establish list of sites to loop through:
sites = unique(df_mothers$SiteLocation.x)

beta_bc <- list()
beta_jac <- list()

for(i in sites){
  df <- df_mothers[df_mothers$SiteLocation.x == i,]
  df <- df[complete.cases(df[,unique(c("Timepoint", "Supp_Status", cov_list_bc[[i]], cov_list_jac[[i]]))]),]
  setBlocks(perm) <- with(df, SubjectID)
  cts <- df[,otus]
  
  beta_bc[[i]] <- adonis2(as.formula(paste("cts ~ Timepoint*Supp_Status + ", 
                                           paste(unique(cov_list_bc[[i]]),collapse=" + "),sep="")), 
                         data = df, permutations = perm, method = "bray")
  
  beta_jac[[i]] <- adonis2(as.formula(paste("cts ~ Timepoint*Supp_Status + ", 
                                            paste(unique(cov_list_bc[[i]]),collapse=" + "),sep="")), 
                     data = df, permutations = perm, method = "jaccard")
}

#Determine p-values:
lapply(beta_bc, function(x) x["Timepoint:Supp_Status", "Pr(>F)"]) 
lapply(beta_jac, function(x) x["Timepoint:Supp_Status", "Pr(>F)"])

#No sites have enough evidence to conclude that interaction between time and supp is associated with beta diversity (either metric)

##No interaction models:
beta_bc_2 <- list()
beta_jac_2 <- list()

for(i in sites){
  df <- df_mothers[df_mothers$SiteLocation.x == i,]
  df <- df[complete.cases(df[,unique(c("Timepoint", "Supp_Status", cov_list_bc[[i]], cov_list_jac[[i]]))]),]
  setBlocks(perm) <- with(df, SubjectID)
  cts <- df[,otus]
  
  beta_bc_2[[i]] <- adonis2(as.formula(paste("cts ~ Timepoint+Supp_Status + ", paste(unique(cov_list_bc[[i]]),collapse=" + "),sep="")), 
                          data = df, permutations = perm, method = "bray")
  
  beta_jac_2[[i]] <- adonis2(as.formula(paste("cts ~ Timepoint+Supp_Status + ", paste(unique(cov_list_bc[[i]]),collapse=" + "),sep="")), 
                           data = df, permutations = perm, method = "jaccard")
}

#Determine time pvalues:
lapply(beta_bc_2, function(x) x["Timepoint", "Pr(>F)"]) 
lapply(beta_jac_2, function(x) x["Timepoint", "Pr(>F)"])

#Determine supp pvalues:
lapply(beta_bc_2, function(x) x["Supp_Status", "Pr(>F)"]) 
lapply(beta_jac_2, function(x) x["Supp_Status", "Pr(>F)"])


#Results output
beta_res <- data.frame(Interaction_Pvalue_BrayCurtis = unlist(lapply(beta_bc, function(x) x["Timepoint:Supp_Status", "Pr(>F)"])),
                       Time_Pvalue_BrayCurtis = unlist(lapply(beta_bc_2, function(x) x["Timepoint", "Pr(>F)"])),
                       Supplement_Pvalue_BrayCurtis = unlist(lapply(beta_bc_2, function(x) x["Supp_Status", "Pr(>F)"])),
                       Interaction_Pvalue_Jaccard = unlist(lapply(beta_jac, function(x) x["Timepoint:Supp_Status", "Pr(>F)"])),
                       Time_Pvalue_Jaccard = unlist(lapply(beta_jac_2, function(x) x["Timepoint", "Pr(>F)"])),
                       Supplement_Pvalue_Jaccard = unlist(lapply(beta_jac_2, function(x) x["Supp_Status", "Pr(>F)"])))

#Save the df for reference
write.csv(beta_res, file = here("Results","Beta_Diversity_2022_02_19.csv"))
