############## Bootstrap Analysis ###############

#### Load Data and packages ####
library(here)
library(ggplot2)
library(gghighlight)

dat <- read.csv(file = here("Data", "Correct_Data_2022_02_19.csv"), header = T)

boots_nams <- list.files(path = here("Results", "OTU", "Bootstrap_Results"))

boots_res <- vector('list', length = length(boots_nams))
for(i in 1:length(boots_nams)){
  boots_res[[i]] <- read.csv(file = here("Results", "OTU", "Bootstrap_Results", boots_nams[i]), header = T)
  names(boots_res)[i] <- boots_nams[i]
}


#### Construct output File ####
out_file <- data.frame(OTU = NA, Site = NA,Estimate_Name = NA, Lower_CI = NA, 
                       Mean_Estimate = NA, Error_Estimate = NA,
                       Upper_CI = NA)

for(i in 1:length(boots_nams)){
  tmp_nam <- strsplit((boots_nams)[i], split = "_")
  
  tmp_OTU <- tmp_nam[[1]][4]
  tmp_site <- gsub(".csv", "",tmp_nam[[1]][8])
  tmp_est <- tmp_nam[[1]][5]
  
  tmp_lower <- quantile(boots_res[[i]]$x, 0.025)
  tmp_mean <- mean(boots_res[[i]]$x)
  tmp_se <- sd(boots_res[[i]]$x)
  tmp_upper <- quantile(boots_res[[i]]$x, 0.975)
  
  out_file <- rbind(out_file, cbind(OTU = tmp_OTU,
                                    Site = tmp_site,
                                    Estimate_Name = tmp_est,
                                    Lower_CI = tmp_lower,
                                    Mean_Estimate = tmp_mean,
                                    Error_Estimate = tmp_se,
                                    Upper_CI = tmp_upper))
}

out_file <- out_file[-1,]

#### Make Plots ####
out_file[,c("Lower_CI", "Mean_Estimate", "Error_Estimate", "Upper_CI")] <- apply(out_file[,c("Lower_CI", "Mean_Estimate","Error_Estimate", "Upper_CI")], 2, as.numeric)

ggplot(out_file[out_file$Estimate_Name == "Time",], aes(x = OTU, y = Mean_Estimate, color = Site)) + theme_bw() +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, linetype = Estimate_Name))

ggplot(out_file[out_file$Estimate_Name == "Time" & out_file$Site == "Guatemala",], aes(x = OTU, y = Mean_Estimate, color = Site)) + theme_bw() +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, linetype = Estimate_Name))

#Smaller plots:
otus <- unique(out_file$OTU)
sub_res <- out_file[out_file$OTU %in% otus[1:5],]

ggplot(sub_res, aes(x = Site, y = Mean_Estimate, color = Site)) + theme_bw() +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, linetype = Estimate_Name)) +
  facet_grid(Estimate_Name~OTU)

#Which CI's do not contain zero? Randomly select 5 of these to plot:
sub_res <- out_file[out_file$OTU %in% sample(unique(out_file[out_file$Lower_CI > 0 | out_file$Upper_CI < 0, "OTU"]), 5),]
ggplot(sub_res, aes(x = Site, y = Mean_Estimate)) + theme_bw() +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), color = ifelse(sub_res$Upper_CI < 0, "red", ifelse(sub_res$Lower_CI >0, "blue", "black")),
                lwd = ifelse(sub_res$Upper_CI <0, 1, ifelse(sub_res$Lower_CI >0, 1, 0.5))) +
  facet_grid(Estimate_Name ~ OTU)

#### Test ####
#Use mean/error to get TS#
out_file$Test_Statistics <- out_file$Mean_Estimate/out_file$Error_Estimate
out_file$P_Value <- 1-pnorm(abs(out_file$Test_Statistics))
out_file$FDR <- p.adjust(out_file$P_Value, method = 'fdr')

# Visualize results
fdr_nams <-unique(out_file[(out_file$FDR<0.05),"OTU"])
fdr_hits <- out_file[out_file$OTU %in% fdr_nams,]
ggplot(fdr_hits, aes(x = Site, y = Mean_Estimate)) + theme_bw() +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), color = ifelse(fdr_hits$Upper_CI < 0, "red", ifelse(fdr_hits$Lower_CI >0, "blue", "black")),
                lwd = ifelse(fdr_hits$Upper_CI <0, 1, ifelse(fdr_hits$Lower_CI >0, 1, 0.5))) +
  facet_grid(Estimate_Name ~ OTU)

#table:
table(fdr_hits$OTU, fdr_hits$P_Value <0.05)

#### Format Output ####
#Make it wide:
library(tidyr)
wide_table <- (pivot_wider(out_file, names_from = c("Site", "Estimate_Name"),
                 values_from = c("Lower_CI", "Mean_Estimate", "Error_Estimate", "Upper_CI",
                                 "Test_Statistics", "P_Value", "FDR")))

write.csv(wide_table, file = here("Results", "OTU", "All_Results_Wide.csv"))

#Add important information to the long format (number of non-zero samples, average abundance)
sites <- unique(dat$SiteLocation.x)
dat[,otus] <- t(apply(dat[,otus], 1, function(x) x/sum(x)))
out_file$Sample_Size <- 0
out_file$Zero_Count <- 0
out_file$Average_relAbundance <- 0
out_file$Max_relAbundance <- 0
for(i in otus){
  for(j in sites){
    samp_size <- length(dat[dat$SiteLocation.x ==j, i])
    n_zeros <- sum(dat[dat$SiteLocation.x == j, i] == 0)
    avg_rel <- mean(dat[dat$SiteLocation.x ==j, i])
    max_rel <- max(dat[dat$SiteLocation.x ==j, i])
    out_file[out_file$OTU ==i & out_file$Site == j, "Sample_Size"] <- samp_size
    out_file[out_file$OTU ==i & out_file$Site == j, "Zero_Count"] <- n_zeros
    out_file[out_file$OTU ==i & out_file$Site == j, "Average_relAbundance"] <- avg_rel
    out_file[out_file$OTU ==i & out_file$Site == j, "Max_relAbundance"] <- max_rel
  }
}

View(out_file)
#Reorder so basic data information is first:
out_file <- out_file[,c("OTU", "Site", "Sample_Size", "Zero_Count", "Average_relAbundance", "Max_relAbundance", "Estimate_Name",
            "Mean_Estimate", "FDR", "Lower_CI", "Upper_CI", "Error_Estimate","Test_Statistics", "P_Value")]
write.csv(out_file, file = here("Results", "OTU", "All_Results_long.csv"))


#### Understand Results ####
library(stringr)
wide_table[,str_detect(colnames(wide_table),pattern = "FDR")] < 0.05

#number of hits by site/time and site (True is a hit, False is not):
table(out_file$Site, out_file$FDR < 0.05, out_file$Estimate_Name)

table(fdr_hits$Site, fdr_hits$FDR < 0.05, fdr_hits$Estimate_Name) #Number of Trues should match last table!

#38 OTUs are of interest, what is overlap?
overlap_time_info <- data.frame(OTU = fdr_nams, DRC = NA, Guatemala = NA, India = NA, Pakistan = NA)
overlap_supp_info <- data.frame(OTU = fdr_nams, DRC = NA, Guatemala = NA, India = NA, Pakistan = NA)

for(i in fdr_nams){
  for(j in sites){
    overlap_time_info[overlap_time_info$OTU == i, j] <- ifelse(fdr_hits[fdr_hits$OTU == i & fdr_hits$Site == j & fdr_hits$Estimate_Name == "Time", "FDR"] < 0.05, "Yes", "No")
    overlap_supp_info[overlap_supp_info$OTU == i, j] <- ifelse(fdr_hits[fdr_hits$OTU == i & fdr_hits$Site == j & fdr_hits$Estimate_Name == "Supplement", "FDR"] < 0.05, "Yes", "No")
  }
}
#Time overlaps:
table(overlap_time_info$DRC, overlap_time_info$Guatemala) #DRC and Guatemala overlap for time
table(overlap_time_info$DRC, overlap_time_info$India) #no overlap DRC/India
table(overlap_time_info$DRC, overlap_time_info$Pakistan) #no overlap DRC/Pakistan
table(overlap_time_info$Guatemala, overlap_time_info$India) #no overlap Guatemala/India
table(overlap_time_info$Guatemala, overlap_time_info$Pakistan) #no overlap Guatemala/Pakistan
table(overlap_time_info$Pakistan, overlap_time_info$India) #no overlap Pakistan/India

#Supplement overlaps:
table(overlap_supp_info$DRC, overlap_supp_info$Guatemala) #no overlap DRC/Guatemala
table(overlap_supp_info$DRC, overlap_supp_info$India) #no overlap DRC/India
table(overlap_supp_info$DRC, overlap_supp_info$Pakistan) #no overlap DRC/Pakistan
table(overlap_supp_info$Guatemala, overlap_supp_info$India) #no overlap Guatemala/India
table(overlap_supp_info$Guatemala, overlap_supp_info$Pakistan) #no overlap Guatemala/Pakistan
table(overlap_supp_info$Pakistan, overlap_supp_info$India) #no overlap Pakistan/India

#Only Guatemala and DRC have overlaps in results, and even then it is only at time!