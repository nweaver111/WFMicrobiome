################## Data Prep ######################

library(here)
library("vegan")
library("cluster")
library("gplots")
library ("ggplot2")
library(RColorBrewer)
library(dplyr)
library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Load Daniel Frank's scripts:
source(here("Code", "D_Frank_Code",'dnf_transforms.R'))
source(here("Code", "D_Frank_Code",'dnf_matrix_corr.R'))
source(here("Code", "D_Frank_Code", 'dnf_misc.R'))
source(here("Code", "D_Frank_Code", 'dnf_multiple_tests.R'))

##### Data Load and Merge #####
#Load the OTU data
df_meta <- (read.metatable(here("Data", "WFSTL_metadata_21Aug2019.txt"),order=TRUE, lib.nms = "Lib"))
df_cts <- read.otutable(here("Data", "WFSTL_alltaxa_cts_16Feb2022.txt"),order = TRUE, rm.prev = FALSE, explicet=TRUE, otu.nms = "OTU_Name", transform = TRUE, meta = df_meta, drop.zeros = TRUE)
df_counts <- as.data.frame(df_cts$cts)
df_counts$Lib <- rownames(df_counts)

#merge meta and cts data:
df_merged <- (merge(df_meta, df_counts, by = "Lib"))

#Create a variable for supplement status:
df_merged$Supp_Status <- 0

for (i in 1:length(df_merged$Lib)) {
  if (df_merged$Timepoint[i] == "12Weeks" & df_merged$Arm[i] == 1) {
    df_merged$Supp_Status[i] <- 1
  }
  if (df_merged$Timepoint[i] == "34Weeks" & df_merged$Arm[i] != 3) {
    df_merged$Supp_Status[i] <- 1
  }
}

#Now load the outcome and diet data:
Moms <- read_excel_allsheets(filename = here("Data", "Women First_Fecal Microbiome_Final Data Set_21-Jul-2020.xlsx")) #renames columns
dem <- Moms$Microbiome_AB_WF01A_Mom_Baby_Me #2732 x 147
diet <- Moms$`Dietary Data` #6860 x 79
diet$Subject_ID

#merge dem information with counts:
df_merged_2 <- left_join(df_merged, dem, by = "AliquotLabel")

rm("Moms", "dem", "df_counts", "df_cts", "df_meta") #clear up some files we no longer need



##### Data Cleaning #####



#Diet data needs to be reduced to single entry per person instead of per meal
str(diet) #Everything is being read as character strings, need to workaround to get correct
# numbers.

write.csv(diet, file = here("Data", "tmp", "WFDiet.csv"), col.names = T)
diet <- read.csv(file = here("Data", "tmp","WFDiet.csv"), header = T, as.is = T, ) #6860 x 80, 

str(diet) #work around works :)

dietvars <- colnames(diet)[c(5:6,8,24:44, 48:57, 59:75, 78)] #Values to keep

#remove rows that contain exactly the same diet information
tmp <- unique(diet[,-c(1,14:19)]) #1608 x 73

#Average dietary data between the two days to get one value per woman
diet_dat <- tmp %>% group_by(Subject_ID) %>% summarize_at(dietvars, mean) #1339 x 53



##Addressing Jennifer's Email: 1 -- merge with 16s data.

# This is point 3 in Jennifer's email. Aliquot Label 113 is removed becuase it was reanalyzed as
# 155. So the databases don't line up this way but it shouldn't matter because the information
# I will be using from the Moms database is not reflected in the sample! 155 doesn't appear in the
# data set!

# Point 2 in Jennifer's email is not applicable here
# Point 3 in Jennifer's email is applied above
# Point 4 in Jennifer's email is about twin data. Only one observation to be concerned about is 615972D.

#Is 615972D in the data set?
dim(df_merged_2[df_merged_2$SubjectID == "615972D",])

#Yes, remove these rows!
df_merged_2 <- (df_merged_2[df_merged_2$SubjectID != "615972D",])

# Point 5 in Jennifer's email deals with DRC data only. Remove miscarriage observations:
df_merged_2 <- (df_merged_2[!(df_merged_2$SubjectID %in% c("236420H", "236421H","241720O","241721M","243020J","243021H")),])

# Point 6 in Jennifer's email is about antibiotic usage before stool sample collection which can be found in the ABX_LAST_MO column:
unique(df_merged_2$ABX_LAST_MO) #More NA's than 1 (1 means yes to antibiotic usage). What is the Z?
length(df_merged_2[df_merged_2$ABX_LAST_MO == "1", "SubjectID"]) #Only 88 observations on antibiotic. 

## Only use this variable for sensitivity analysis! ##

# Point 7 in Jennifer's email is about the dietary data which I addressed at the start of 
#  this file. Merge the data in now:

df_final <- merge(df_merged_2, diet_dat[,c(1,5:53)], by.x = "SubjectID", by.y = "Subject_ID")
dim(df_final) # 2732 x 1258.

# Point 8 doesn't need any edits. It is about education data.
# Point 9 in Jennifer's email is about extra variables she included in the data set. Will
#  ignore these for this analysis
# Point 10 in Jennifer's email is about supplement and inflammation data. We will want to
#   use these in our analysis as well!
# Point 11 in Jennifer's email is about the birth outcomes given in the data set. These
#  are not of use right now, so I will remove them from the data base when we decide which
#  variables we actually want to use moving forward.


#Load the mother data


#remove the excess files to clear space
rm("df_merged", "df_merged_2", "diet", "diet_dat", "tmp")


##### Data Prep #####

#Subset to only mothers at 12 and 34 weeks:
df_final <- df_final[df_final$Mother_Infant.x == "Mother",] #1708 x 1258
df_final <- df_final[df_final$Timepoint %in% c("12Weeks", "34Weeks"),] #1028 x 1258

#Remove any mother with Read count < 3700
df_final <- (df_final[!(df_final$Reads < 3700),]) #983 x 1258

#Keep OTUs that were seen in at least 5% of samples (across all sites)
otus <- colnames(df_final)[21:1062]
otu_presence <- apply(df_final[,otus], 2, function(x) sum(x != 0)/length(df_final[,1]))

sum(otu_presence > 0.05) #194 OTUs pass the threshold!

otus_keep <- names(which(otu_presence>0.05)) #Store the OTUs we want to keep

#Keep OTUs that had a relative abundance > 0.001 in at least 1 mother
otu_abundance <- as.data.frame(t(apply(df_final[,otus], 1, function(x) (x/sum(x))))) #get relative abundances
otu_abund_track <- apply(otu_abundance, 2, function(x) sum(x > 0.001)) #Determine how many samples have rel abundance above 0.001 for each OTU
otus_potential <- names(which(otu_abund_track != 0)) #Which taxa have 0 samples above 0.001?

otus <- intersect(otus_keep, otus_potential) #What is the intersection of the two filters?

#final OTUs list
length(otus)


##### Data Output #####
df_final <- df_final[,c(colnames(df_final)[1:20], otus, colnames(df_final)[1063:ncol(df_final)])] #983 x 374

write.csv(df_final, file = here("Data", "Correct_Data_2022_02_19.csv"), row.names = F)


#OTU list:
otus <- gsub(":", ".", otus) #remove colon issue
otus <- gsub("-", ".", otus) #remove dash issue
saveRDS(otus, file = here("Data", "otu_list_2022_02_19.rds"))
