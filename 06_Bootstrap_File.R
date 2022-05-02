##Function to perform a bootstrap model

bootstrap <- function(dat, otu, site, b_div) {
  #Resample (with replacement) to construct a new data set of same size:
  dat <- dat[dat$SiteLocation.x == site,]
  m_list <- unique(dat$SubjectID)
  m <- length(m_list)#How many unique subjects in original data?
  tmp_index <- sample(m_list, size = m, replace = T) #Each index is equally likely, and resampling is allowed
  
  #Need a list of rownames to be able to create bootstrapped data frame with longitudinal values!
  row_list <- c() #List of row numbers
  new_ID <- c() #List to store new subject ID number for data set!

  for(i in 1:length(tmp_index)){
    row_list <- c(row_list, rownames(dat[dat$SubjectID == tmp_index[i],]))
    new_ID <- c(new_ID, rep.int(paste("Subjet_",i,sep=""),nrow(dat[dat$SubjectID == tmp_index[i],])))
  }
  
  dat_tmp <- dat[row_list,] 
  dat_tmp$BootstrapID <- new_ID  
  
  #Run models by site:
  rse <- unique(c("Supp_Status", "Timepoint", b_div[[site]]))
  tmp <- glm.nb(as.formula(paste(otu, " ~ ",paste(rse, collapse="+"), sep = "")), data = dat_tmp)
  
  a <- summary(tmp)
  
  return(c(a$coefficients["Timepoint34Weeks", "Estimate"], a$coefficients["Supp_Status1", "Estimate"]))
}

library(lme4)
df_mothers <- read.csv("C:\\Users\\nweav\\OneDrive\\Documents\\WF Microbiome\\2022_02_15_Corrected_Analysis\\2022_WF_Microbiome_Corrected_Analysis\\Data\\Correct_Data_2022_02_19.csv")
bdiv <- readRDS("C:\\Users\\nweav\\OneDrive\\Documents\\WF Microbiome\\2022_02_15_Corrected_Analysis\\2022_WF_Microbiome_Corrected_Analysis\\Data\\Covariate\\brayCurtis_covariates.rds")
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


#establish parameters
B = 1000
otu = colnames(df_mothers)[43]
dat = df_mothers
site = "Guatemala"

#create lists to store values
time_est <- c()
supp_est <- c()


###loop for sampling
tracker <- 0
while(length(supp_est) < B){
  tryCatch({
    res <- bootstrap(dat = df_mothers, otu = otu, site = "Guatemala", b_div = bdiv)
    time_est <- c(time_est, res[1])
    supp_est <- c(supp_est, res[2])
    tracker <- tracker + 1
    print(tracker)
  }, error = function(e){
    B <- B-1
    tracker <- tracker + 1
    print(paste("No model for iteration: ", tracker,sep=""))
    #Just want to make sure I acquire exactly B models
  })
}

