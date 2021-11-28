##### Models #####

##Function to perform a bootstrap model
bootstrap <- function(dat, otu, site) {
  #Resample (with replacement) to construct a new data set of same size:
  m <- length(dat[,1])
  tmp_index <- sample(1:m, size = m, replace = T) #Each index is equally likely, and resampling is allowed
  dat_tmp <- dat[tmp_index,] 
  
  if(site == "Pakistan"){
    tmp <- glmer.nb(as.formula(paste(otu, " ~ Supp_Status + Timepoint + BatchAdj + F1AA15 +
                      floor_01 + impr_cook_fuel_01 + stdreads + (1|SubjectID)")), data = dat_tmp)
  }
  
  if(site == "Guatemala"){
    tmp <- glmer.nb(as.formula(paste(otu, " ~ BatchAdj + Timepoint + Supp_Status + F1AA16 + F1AA17 + SES_common + stdreads + (1|SubjectID)")),
                    data = dat_tmp)
  }
  
  if(site == "India"){
    tmp <- glmer.nb(as.formula(paste(otu, " ~ Supp_Status + Timepoint + BatchAdj + parity_continuous + stdreads + 
                      (1|SubjectID)")), data = dat_tmp)
  }
  
  if(site == "DRC"){
    tmp <- glmer.nb(as.formula(paste(otu, " ~ Supp_Status + Timepoint + BatchAdj + F1AA16 + SES_common + stdreads + (1|SubjectID)")), 
                    data = dat_tmp)
  }
  a <- summary(tmp)
  
  return(c(a$coefficients["Timepoint34Weeks", "Estimate"], a$coefficients["Supp_Status1", "Estimate"]))
}

#establish parameters
B = n_samples
otu = otu
dat = tmp_dat
site = site

#create lists to store values
time_est <- c()
supp_est <- c()

#Load Library
library(lme4)

### Make sure covariates are factor or numeric for all 4 site models
dat[,"Supp_Status"] <- as.factor(dat[,"Supp_Status"])
dat[,"Timepoint"] <- as.factor(dat[,"Timepoint"])
dat[,"BatchAdj"] <- as.factor(dat[,"BatchAdj"])
dat[,"floor_01"] <- as.factor(dat[,"floor_01"])
dat[,"impr_cook_fuel_01"] <- as.factor(dat[,"impr_cook_fuel_01"])
dat[,"ClusterAdj"] <- as.factor(dat[,"ClusterAdj"])
dat[,"hh_assets_01"] <- as.factor(dat[,"hh_assets_01"])
dat[,"impr_water_01"] <- as.factor(dat[,"impr_water_01"])

###loop for sampling
tracker <- 0
while(length(supp_est) < B){
  tryCatch({
    res <- bootstrap(dat = dat, otu = otu, site = site)
    time_est <- c(time_est, res[1])
    supp_est <- c(supp_est, res[2])
    tracker <- tracker + 1
    print(tracker)
  }, error = function(e){
    #Just want to make sure I acquire exactly B models
  })
}

##Save results to specific directory
setwd(paste("/newhome/weavenic/WF_Microbiome/Oct_30/Bootstraps/TopHits/", site, sep=""))
write.csv(time_est, file = paste("Bootstrap Models for ", otu, " Time Estimate in ", site, ".csv", sep =""))
write.csv(supp_est, file = paste("Bootstrap Models for ", otu, " Supplement Estimate in ", site, ".csv", sep = ""))