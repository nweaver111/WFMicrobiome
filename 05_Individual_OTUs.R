#Load Library
library(here)
library(lme4)

##Load important information:
df_mothers <- read.csv(here("Data", "Correct_Data_2022_02_19.csv"))
covars_bc <- readRDS(here("Data", "Covariate", "brayCurtis_covariates.rds"))
covars_jac <- readRDS(here("Data", "Covariate", "jaccard_covariates.rds"))

otus <- readRDS(here("Data", "otu_list_2022_02_19.rds"))


##create lists to store values
time_est <- c()
supp_est <- c()
interaction_est <- c()

## Make sure covariates are factor or numeric for all 4 site models
df_mothers$SES_common <- as.numeric(df_mothers$SES_common)
df_mothers$BatchAdj <- as.factor(df_mothers$BatchAdj)
df_mothers$F1AA16 <- as.numeric(df_mothers$F1AA16)
df_mothers$impr_water_01 <- as.factor(df_mothers$impr_water_01)
df_mothers$deliv_date <- as.numeric(df_mothers$deliv_date)
df_mothers$Supp_Status <- as.factor(df_mothers$Supp_Status)
df_mothers$std_reads <- as.numeric(df_mothers$std_reads)
df_mothers$s1_comp <- as.numeric(df_mothers$s1_comp)
df_mothers$impr_cook_fuel_01 <- as.factor(df_mothers$impr_cook_fuel_01)


## Run model on data for all sites x OTU pairings using bray-curtis covariates:
sites <- unique(df_mothers$SiteLocation.x)
tracker <- 0
for(i in sites){
  for(j in otus){
    tmp_dat <- df_mothers[df_mothers$SiteLocation.x == i,]
    
    tryCatch({
      tmp_mod <- glmer.nb(as.formula(paste(j, " ~ Supp_Status + Timepoint + ",paste(covars_bc[[i]], collapse = " + "),
                                           "+ (1|SubjectID)", sep = "")), data = tmp_dat)
      sum_mod <- summary(tmp_mod)
      
      time_est <- c(time_est, rbind(i, j, sum_mod$coefficients["Timepoint34Weeks",]))
      supp_est <- c(supp_est, rbind(i, j, sum_mod$coefficients["Supp_Status1",]))
    }, error = function(e){
      time_est <- cbind(time_est, rbind(i, j, ""))
      supp_est <- cbind(supp_est, rbind(i, j, ""))
    })
    
    tracker <- tracker + 1
    print(tracker)
  }
}

##Save results to specific directory
write.csv(time_est, file = here("Results", "OTU", "Time_Estimates.csv"))
write.csv(supp_est, file = here("Results", "OTU", "Supplement_Estimates.csv"))
