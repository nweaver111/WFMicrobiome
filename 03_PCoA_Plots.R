############### PCoA Analysis #################

library(here)
library(ggplot2)
library(vegan)
library(ape)

##### Data Input #####

#Dataset:
df_mothers <- read.csv(here("Data", "Correct_Data_2022_02_19.csv"))
rownames(df_mothers) <- df_mothers$Lib

otus <- readRDS(here("Data", "otu_list_2022_02_19.rds"))

#Covariate Lists:
cov_list_bc <- readRDS(here("Data", "Covariate", "brayCurtis_covariates.rds"))
cov_list_jac <- readRDS(here("Data", "Covariate", "jaccard_covariates.rds"))


##### Oridnation Plots #####

# Calculate beta diversity values

#Need proportions, not counts:
rel_dat <- (t(apply(df_mothers[,otus], 1, function(x) x/sum(x))))

#Bray-Curtis metric
bc_reldist <- vegdist(rel_dat)

#Jaccard:
jac_reldist <- vegdist(rel_dat, method = "jaccard")


# Calculate PCoAs

pcoa_bc <- pcoa(bc_reldist, rn = rownames(df_mothers))
pcoa_jac <- pcoa(jac_reldist,rn = rownames(df_mothers))


# Visualize the PCoA Information

#Relative Eigenvalue plot:
df_bc <- data.frame(Relative_Eigenvalue = pcoa_bc$values[,"Rel_corr_eig"], PCoA = 1:max(dim(df_mothers)))

ggplot(df_bc, aes(x = PCoA, y = Relative_Eigenvalue)) + geom_point() + theme_bw() +
  labs(x = "Principal Component Axis (ordered)",
       y = "% of Variance Explained",
       title = "Percent of Variance Explained by PCoA\nBray-Curtis\nAll Sites and Timepoints Combined")
ggsave(filename = here("Plots", "Ordination", "Prop_Var_Explained_BC.jpg"),
       width = 8.77, height = 7.04)

#Cumulative Sum plot
cumsum_bc <- c()
for(i in 1:max(dim(df_mothers))){
  cumsum_bc <- append(cumsum_bc, sum(pcoa_bc$values[1:i,"Rel_corr_eig"]))
}

df_bc_cumsum <- data.frame(Cumulative_Sum = cumsum_bc, PCoA = 1:max(dim(df_mothers)))
ggplot(df_bc_cumsum, aes(x = PCoA, y = Cumulative_Sum)) + geom_point() + theme_bw()
ggsave(filename = here("Plots", "Ordination", "Cumulative_Sum_Var_Explained_BC.jpg"),
       width = 8.77, height = 7.04)

# There really isn't anything too helpful here...

# First two PCoAs:
df <- data.frame(PCoA.1 = pcoa_bc$vectors[,1], PCoA.2 = pcoa_bc$vectors[,2],
                 Lib = rownames(pcoa_bc$vectors))


df <- merge(df, df_mothers[,c("Lib","SubjectID", "Timepoint","SiteLocation.x","Supp_Status","Arm.x","BatchAdj",
                              "ShannonE_3700Seqs")],
            by = "Lib")

dim(df)
View(df)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = SiteLocation.x)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_bc$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_bc$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Bray-Curtis\nColored by Site")+
  theme(legend.position = "bottom")+ stat_ellipse()
ggsave(filename = here("Plots", "Ordination", "PCoA_all_Site_BC.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = SiteLocation.x, shape = Timepoint)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_bc$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_bc$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Bray-Curtis\nColored by Site and Timepoint")+
  theme(legend.position = "bottom")+ stat_ellipse(aes(lty = Timepoint))
ggsave(filename = here("Plots", "Ordination", "PCoA_all_SiteandTimepoint_BC.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = ShannonE_3700Seqs)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_bc$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_bc$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Bray-Curtis\nColored by Alpha Diversity (Shannon E)")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_ShannonE_All_Site_BC.jpg"),
       width = 8.77, height = 7.04)


#Jaccard:
#Relative Eigenvalue plot:
df_jac <- data.frame(Relative_Eigenvalue = pcoa_jac$values[,"Rel_corr_eig"], PCoA = 1:max(dim(df_mothers)))

ggplot(df_jac, aes(x = PCoA, y = Relative_Eigenvalue)) + geom_point() + theme_bw() +
  labs(x = "Principal Component Axis (ordered)",
       y = "% of Variance Explained",
       title = "Percent of Variance Explained by PCoA\nJaccard\nAll Sites and Timepoints Combined")
ggsave(filename = here("Plots", "Ordination", "Prop_Var_Explained_Jaccard.jpg"),
       width = 8.77, height = 7.04)

#Cumulative Sum plot
cumsum_jac <- c()
for(i in 1:max(dim(df_mothers))){
  cumsum_jac <- append(cumsum_jac, sum(pcoa_jac$values[1:i,"Rel_corr_eig"]))
}

df_jac_cumsum <- data.frame(Cumulative_Sum = cumsum_jac, PCoA = 1:max(dim(df_mothers)))
ggplot(df_jac_cumsum, aes(x = PCoA, y = Cumulative_Sum)) + geom_point() + theme_bw()
ggsave(filename = here("Plots", "Ordination", "Cumulative_Sum_Var_Explained_Jaccard.jpg"),
       width = 8.77, height = 7.04)

# There really isn't anything too helpful here...

# First two PCoAs:
df <- data.frame(PCoA.1 = pcoa_jac$vectors[,1], PCoA.2 = pcoa_jac$vectors[,2],
                 Lib = rownames(pcoa_jac$vectors))


df <- merge(df, df_mothers[,c("Lib","SubjectID", "Chao1_3700Seqs", "Sobs_3700Seqs", "ShannonE_3700Seqs", "ShannonH_3700Seqs",
                              "Timepoint","SiteLocation.x","Supp_Status","Arm.x","BatchAdj")],
            by = "Lib")

dim(df)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = SiteLocation.x)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Site")+
  theme(legend.position = "bottom")+ stat_ellipse()
ggsave(filename = here("Plots", "Ordination", "PCoA_all_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = SiteLocation.x, shape = Timepoint)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Site and Timepoint")+
  theme(legend.position = "bottom")+ stat_ellipse(aes(lty = Timepoint))
ggsave(filename = here("Plots", "Ordination", "PCoA_all_SiteandTimepoint_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = Chao1_3700Seqs)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Alpha Diversity (Chao1)")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_Chao1_All_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = Sobs_3700Seqs)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Alpha Diversity (Sobs)")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_Sobs_All_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = ShannonE_3700Seqs)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Alpha Diversity (Shannon E)")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_ShannonE_All_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = ShannonH_3700Seqs)) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Alpha Diversity (Shannon H)")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_ShannonH_All_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = as.factor(Supp_Status))) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Supplement Status")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_Supplement_all_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = as.factor(Arm.x))) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Arm")+
  theme(legend.position = "bottom")
ggsave(filename = here("Plots", "Ordination", "PCoA_Arm_all_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = as.factor(Timepoint))) + theme_bw() + geom_point() +
  labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
       y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
       title = "PCoA 1 vs PCoA 2 using Jaccard\nColored by Timepoint")+
  theme(legend.position = "bottom") + stat_ellipse()
ggsave(filename = here("Plots", "Ordination", "PCoA_Timepoint_all_Site_Jaccard.jpg"),
       width = 8.77, height = 7.04)

#Maybe one of the covariates??
unique(unlist(cov_list_bc))
df_mothers$SES_common <- as.numeric(df_mothers$SES_common)
df_mothers$BatchAdj <- as.factor(df_mothers$BatchAdj)
df_mothers$F1AA16 <- as.numeric(df_mothers$F1AA16)
df_mothers$impr_water_01 <- as.factor(df_mothers$impr_water_01)
df_mothers$deliv_date <- as.numeric(df_mothers$deliv_date)
df_mothers$Supp_Status <- as.factor(df_mothers$Supp_Status)
df_mothers$std_reads <- as.numeric(df_mothers$std_reads)
df_mothers$s1_comp <- as.numeric(df_mothers$s1_comp)
df_mothers$impr_cook_fuel_01 <- as.factor(df_mothers$impr_cook_fuel_01)

for(i in unique(unlist(cov_list_jac))){
  
  df <- data.frame(PCoA.1 = pcoa_jac$vectors[,1], PCoA.2 = pcoa_jac$vectors[,2],
                   Lib = rownames(pcoa_jac$vectors))
  
  df <- merge(df, df_mothers[,c("Lib","SubjectID", i, "Timepoint","SiteLocation.x","Supp_Status","Arm.x","BatchAdj")],
              by = "Lib")
  
  ggplot(df, aes(x = PCoA.1, y = PCoA.2, color = df[,i])) + theme_bw() + geom_point() +
    labs(x = paste("Principal Coordinate Axis 1 (",round(pcoa_jac$values[1,"Rel_corr_eig"],4)*100,"%)", sep =""),
         y = paste("Principal Coordinate Axis 2 (",round(pcoa_jac$values[2,"Rel_corr_eig"],4)*100,"%)", sep =""),
         title = paste("PCoA 1 vs PCoA 2 using Jaccard\nColored by ", i,sep=""))+
    theme(legend.position = "bottom")
  ggsave(filename = here("Plots", "Ordination", paste("PCoA_",i,"_all_Site_Jaccard.jpg",sep="")),
         width = 8.77, height = 7.04)
}

#Shannon Diversity is explaining the separation on PCoA1, but that is expected. Reviewer asked for us to do this,
# but in my opinion it is redundant. We want to visualize our predictors here, we just can't do it well!
