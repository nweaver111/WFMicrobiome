########################## Extra Numbers in Paper #############################

#### Libraries ####
library(here)
library(stringr)

#### Load Data ####
df <- read.csv(here("Data", "Correct_Data_2022_02_19.csv"))


## Distribution of Read depths ##
summary(df$Reads)
IQR(df$Reads)


## Covariates by Site ##
covars <- readRDS(here("Data", "Covariate", "brayCurtis_covariates.rds"))

covars$DRC
covars$Guatemala
covars$India
covars$Pakistan


#### Avg Relative Abundance across all 4 sites ####
otus <- readRDS(here("Data", "otu_list_2022_02_19.rds"))

df_rel <- df
df_rel[,otus] <- t(apply(df_rel[,otus], 1, function(x) x/sum(x)))

avg_rel_abund <- apply(df_rel[,otus], 2, mean)

df_rel <- as.data.frame(avg_rel_abund)
df_rel$OTU <- rownames(df_rel)
write.csv(df_rel, here("Results", "Avg_relative_abundance_all_sites.csv"))

#### Table 4 ####
firms <- colnames(df)[str_detect(colnames(df), "Firm")]
acti <- colnames(df)[str_detect(colnames(df), "Acti")]
bact <- colnames(df)[str_detect(colnames(df), "Bact")]
prot <- colnames(df)[str_detect(colnames(df), "Prot")]

table_4 <- data.frame(Subject = NA, Batch = NA, Site = NA, Time = NA, Firm = NA, Acti = NA, Bact = NA, Prot = NA, Other = NA)
taxa <- otus[1:9]

#Taxa we want to plot:
#taxa <- c("Bifidobacterium",
#          "Prevotella",
#          "Lachnospiraceae",
#          "Ruminococcaceae",
#          "Escherichia_Shigella",
#          "Streptococcus",
#          "Faecalibacterium",
#          "Pseudobutyrivibrio",
#          "Succinivibrio")

bar_df <- data.frame(Site = NA, Time = NA, Taxa = NA, Relative_Abundance = NA)

for(i in unique(df$SiteLocation.x)){
  for(j in c("12Weeks", "34Weeks")){
    rel_abund <- t(apply((df[df$SiteLocation.x==i & df$Timepoint == j,otus]), 1, function(x) x/sum(x)))
    mean_rel <- apply(rel_abund, 2, mean)
    table_4 <- rbind(table_4, cbind(Subject = df[df$SiteLocation.x==i & df$Timepoint==j,"SubjectID"],
                           Batch = df[df$SiteLocation.x==i & df$Timepoint==j,"BatchAdj"],
                           Site = i,
                           Time = j,
                           Firm = (apply(rel_abund[,firms], 1, sum)),
                           Acti = (apply(rel_abund[,acti],1,sum)),
                           Bact = (apply(rel_abund[,bact],1,sum)),
                           Prot = (apply(rel_abund[,prot],1,sum)),
                           Other = (apply(rel_abund[,!(names(rel_abund) %in% c(firms,acti,bact,prot))],1,sum))))
    for(k in taxa){
      tmp_val <- mean_rel[k]
      
      bar_df <- rbind(bar_df, cbind(Site = i, 
                                    Time = j,
                                    Taxa = k,
                                    Relative_Abundance = as.numeric(tmp_val)))
    }
    bar_df$Relative_Abundance <- as.numeric(bar_df$Relative_Abundance)
    bar_df <- rbind(bar_df, cbind(Site = i,
                                  Time = j,
                                  Taxa = "Other",
                                  Relative_Abundance = 1-sum(bar_df[bar_df$Site==i & bar_df$Time==j,"Relative_Abundance"], na.rm = T)))
  }
}
table_4 <- table_4[-1,]
View(table_4)

##Perform statistical tests on phylum level data
library(lmerTest)
library(emmeans)

table_4$Site <- as.factor(table_4$Site)
table_4$Batch <- as.factor(table_4$Batch)
table_4$Time <- as.factor(table_4$Time)

lm_noB_vals <- lm_b_vals <- re_noB_vals <- re_b_vals  <- vector('list')
mod_comps <- data.frame(Outcome = NA, Batch = NA, Model = NA, Interaction_Name = NA, Interaction_Effect = NA, Interaction_Pvalue = NA)
aov_comps <- data.frame(Outcome = NA, Batch = NA, Model = NA, Name = NA, F_value = NA, P_value = NA)
#Create interaction models:
for(i in c("Firm", "Acti", "Bact", "Prot")){
  table_4[,i] <- as.numeric(table_4[,i])
  lm_noB_vals[[i]] <- lmer(as.formula(paste(i, "~ Site+Time+(1|Subject)",sep="")),
                         data = table_4)
  lm_b_vals[[i]] <- lmer(as.formula(paste(i, "~ Site+Time + Batch +(1|Subject)",sep="")),
                         data = table_4) 
  re_noB_vals[[i]] <- lmer(as.formula(paste(i,"~Site*Time + (1|Subject)",sep="")),
                           data = table_4)
  re_b_vals[[i]] <- lmer(as.formula(paste(i,"~Site*Time + Batch + (1|Subject)",sep="")),
                           data = table_4)
  
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "lm",
                                      Interaction_Name = "Time34Weeks",
                                      Interaction_Effect = summary(lm_noB_vals[[i]])$coefficients["Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(lm_noB_vals[[i]])$coefficients["Time34Weeks","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "lm",
                                      Interaction_Name = 'SitePakistan',
                                      Interaction_Effect = summary(lm_noB_vals[[i]])$coefficients["SitePakistan","Estimate"],
                                      Interaction_Pvalue = summary(lm_noB_vals[[i]])$coefficients["SitePakistan","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "lm",
                                      Interaction_Name = 'SiteIndia',
                                      Interaction_Effect = summary(lm_noB_vals[[i]])$coefficients["SiteIndia","Estimate"],
                                      Interaction_Pvalue = summary(lm_noB_vals[[i]])$coefficients["SiteIndia","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "lm",
                                      Interaction_Name = 'SiteGuatemala',
                                      Interaction_Effect = summary(lm_noB_vals[[i]])$coefficients["SiteGuatemala","Estimate"],
                                      Interaction_Pvalue = summary(lm_noB_vals[[i]])$coefficients["SiteGuatemala","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "lm",
                                      Interaction_Name = "Time34Weeks",
                                      Interaction_Effect = summary(lm_b_vals[[i]])$coefficients["Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(lm_b_vals[[i]])$coefficients["Time34Weeks","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "lm",
                                      Interaction_Name = 'SitePakistan',
                                      Interaction_Effect = summary(lm_b_vals[[i]])$coefficients["SitePakistan","Estimate"],
                                      Interaction_Pvalue = summary(lm_b_vals[[i]])$coefficients["SitePakistan","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "lm",
                                      Interaction_Name = 'SiteIndia',
                                      Interaction_Effect = summary(lm_b_vals[[i]])$coefficients["SiteIndia","Estimate"],
                                      Interaction_Pvalue = summary(lm_b_vals[[i]])$coefficients["SiteIndia","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "lm",
                                      Interaction_Name = 'SiteGuatemala',
                                      Interaction_Effect = summary(lm_b_vals[[i]])$coefficients["SiteGuatemala","Estimate"],
                                      Interaction_Pvalue = summary(lm_b_vals[[i]])$coefficients["SiteGuatemala","Pr(>|t|)"]))
  
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "RE",
                                      Interaction_Name = "SiteGuatemala:Time34Weeks",
                                      Interaction_Effect = summary(re_noB_vals[[i]])$coefficients["SiteGuatemala:Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(re_noB_vals[[i]])$coefficients["SiteGuatemala:Time34Weeks","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "RE",
                                      Interaction_Name = 'SitePakistan:Time34Weeks',
                                      Interaction_Effect = summary(re_noB_vals[[i]])$coefficients["SitePakistan:Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(re_noB_vals[[i]])$coefficients["SitePakistan:Time34Weeks","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "No",
                                      Model = "RE",
                                      Interaction_Name = 'SiteIndia:Time34Weeks',
                                      Interaction_Effect = summary(re_noB_vals[[i]])$coefficients["SiteIndia:Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(re_noB_vals[[i]])$coefficients["SiteIndia:Time34Weeks","Pr(>|t|)"]))
  
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "RE",
                                      Interaction_Name = "SiteGuatemala:Time34Weeks",
                                      Interaction_Effect = summary(re_b_vals[[i]])$coefficients["SiteGuatemala:Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(re_b_vals[[i]])$coefficients["SiteGuatemala:Time34Weeks","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "RE",
                                      Interaction_Name = 'SitePakistan:Time34Weeks',
                                      Interaction_Effect = summary(re_b_vals[[i]])$coefficients["SitePakistan:Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(re_b_vals[[i]])$coefficients["SitePakistan:Time34Weeks","Pr(>|t|)"]))
  mod_comps <- rbind(mod_comps, cbind(Outcome = i,
                                      Batch = "Yes",
                                      Model = "RE",
                                      Interaction_Name = 'SiteIndia:Time34Weeks',
                                      Interaction_Effect = summary(re_b_vals[[i]])$coefficients["SiteIndia:Time34Weeks","Estimate"],
                                      Interaction_Pvalue = summary(re_b_vals[[i]])$coefficients["SiteIndia:Time34Weeks","Pr(>|t|)"]))

  #Comparison of ANOVAs:
  aov_comps <- rbind(aov_comps,cbind(Outcome = i,
                                     Batch = "Yes",
                                     Model = "Interaction",
                                     Name = "Interaction",
                                     F_value = anova(re_b_vals[[i]])["Site:Time","F value"],
                                     P_value = anova(re_b_vals[[i]])["Site:Time","Pr(>F)"]))
  
  aov_comps <- rbind(aov_comps,cbind(Outcome = i,
                                     Batch = "Yes",
                                     Model = "Additive",
                                     Name = "Time",
                                     F_value = anova(lm_b_vals[[i]])["Time","F value"],
                                     P_value = anova(lm_b_vals[[i]])["Time","Pr(>F)"]))
  
  aov_comps <- rbind(aov_comps,cbind(Outcome = i,
                                     Batch = "Yes",
                                     Model = "Additive",
                                     Name = "Site",
                                     F_value = anova(lm_b_vals[[i]])["Site","F value"],
                                     P_value = anova(lm_b_vals[[i]])["Site","Pr(>F)"]))
  
}

emmeans(re_b_vals[[i]],pairwise ~ Site*Time)
mod_comps <- mod_comps[-1,]
aov_comps <- aov_comps[-1,]
View(aov_comps)

#paired comparisons if variable is associated:
aov_comps$Compare <- ifelse(as.numeric(aov_comps$P_value) < 0.05, "Yes", "No")
write.csv(aov_comps, here("Results", "Table_4", "Anova_Comparisons.csv"))
tracker <- 1
emmeans_obj <- vector('list')
for(i in 1:nrow(aov_comps)){
  
  if(aov_comps[i,"Compare"] == "Yes"){
    outcome <- aov_comps[i,"Outcome"]
    eff_nam <- aov_comps[i,"Name"]
    if(aov_comps[i,"Model"] == "Additive"){
      model <- lm_b_vals[[outcome]]
      emmeans_obj[[tracker]] <- emmeans(model, as.formula(paste("pairwise ~ ",eff_nam,sep="")))
    }

    if(aov_comps[i,"Model"] == "Interaction"){
      model <- re_b_vals[[outcome]]
      emmeans_obj[[tracker]] <- emmeans(model,pairwise~Site*Time)
    }
                                     
    tracker <- tracker + 1
  }
  
}

saveRDS(emmeans_obj, here("Results", "Table_4","Paired_Comps","All_comparisons.rds"))



#### Figure 3 ####

bar_df <- bar_df[-1,]

bar_df$Relative_Abundance <- as.numeric(bar_df$Relative_Abundance)
library(ggplot2)
library(viridis)

bar_df$Taxa <- as.factor(bar_df$Taxa)


bar_df$Taxa2 <- reorder(bar_df$Taxa, bar_df[,"Relative_Abundance"])
bar_df$Taxa3 <- relevel(bar_df$Taxa2, "Other")
ggplot(bar_df, aes(x = Time, fill = Taxa3, y = Relative_Abundance)) + theme_bw() +
  geom_bar(position = "fill", stat = "identity") + facet_grid(cols = vars(Site)) + 
  theme(legend.position = "bottom") +
  scale_x_discrete(labels=c("12Weeks" = "12 Weeks", "34Weeks" = "34 Weeks"))+
  labs(x = "Time",
       y = "Relative Abundance", fill = "Taxa") + scale_fill_brewer(type = "qual", palette = "Paired")
ggsave(filename = here("Plots", "rel_abundance_bar_charts.png"),
       width = 10, height = 6)


#### Figure 1 ####
#Alpha Diversity plots:
sub_df <- df[,c("SubjectID", "Timepoint", "SiteLocation.x", "Supp_Status","Chao1_3700Seqs", 
                "ShannonE_3700Seqs", "ShannonH_3700Seqs", "Sobs_3700Seqs")]

library(tidyr)
alpha_df <- pivot_longer(sub_df, cols = c("Chao1_3700Seqs", "ShannonE_3700Seqs", "ShannonH_3700Seqs", "Sobs_3700Seqs"),
                         names_to = "Alpha_Diversity", values_to = "Index")

#Fix lables/names:
alpha_df[alpha_df$Alpha_Diversity == "ShannonE_3700Seqs", "Alpha_Diversity"] <- "Shannon H/Hmax"
alpha_df[alpha_df$Alpha_Diversity == "ShannonH_3700Seqs", "Alpha_Diversity"] <- "Shannon H"
alpha_df[alpha_df$Alpha_Diversity == "Chao1_3700Seqs", "Alpha_Diversity"] <- "Chao 1"
alph_df <- alpha_df[alpha_df$Alpha_Diversity %in% c("Shannon H/Hmax", "Shannon H", "Chao 1"),]
alph_df$Supp_Status <- ifelse(alph_df$Supp_Status == 0, "No", "Yes")
alph_df$Timepoint <- ifelse(alph_df$Timepoint == "12Weeks", "12 Weeks", "34 Weeks")
alph_df$Timepoint <- ifelse(alph_df$Timepoint == "12 Weeks", "12", "34")

dodge <- position_dodge(width = 0.8)

#regular boxplots:
ggplot(alph_df, aes(x=Timepoint, y = Index, fill = as.factor(Supp_Status))) + theme_bw() + geom_boxplot() +
  facet_grid(Alpha_Diversity~SiteLocation.x, scales = "free_y",
             switch = "y",
             labeller = as_labeller(c(`Chao 1` = "Chao 1 Index", `Shannon H` = "Shannon H Index", 
                                      `Shannon H/Hmax` = "Shannon H/Hmax Index", Guatemala = "Guatemala",
                                      Pakistan = "Pakistan", DRC = "DRC", India = "India"))) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(fill = "On Supplement?", y = "", x = "Time (Weeks)",
       title = "Alpha Diversity Comparisons")

#Violin plots:
ggplot(alph_df, aes(x=Timepoint, y = Index, fill = as.factor(Supp_Status))) + theme_bw() + geom_violin(alpha = 0.4, position = dodge) +
  facet_grid(Alpha_Diversity~SiteLocation.x, scales = "free_y",
             switch = "y",
             labeller = as_labeller(c(`Chao 1` = "Chao 1 Index", `Shannon H` = "Shannon H Index", 
                                      `Shannon H/Hmax` = "Shannon H/Hmax Index", Guatemala = "Guatemala",
                                      Pakistan = "Pakistan", DRC = "DRC", India = "India"))) + 
  geom_boxplot(alpha = 0.4, width = 0.2, position = dodge, outlier.colour = NA) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(fill = "On Supplement?", y = "", x = "Time (Weeks)",
       title = "Alpha Diversity Comparisons")

#box plots with scatterplot overlay:
ggplot(alph_df, aes(x=Timepoint, y = Index, fill = as.factor(Supp_Status))) + theme_bw() +
  facet_grid(Alpha_Diversity~SiteLocation.x, scales = "free_y", switch = "y",
             labeller = as_labeller(c(`Chao 1` = "Chao 1 Index", `Shannon H` = "Shannon H Index", 
                                      `Shannon H/Hmax` = "Shannon H/Hmax Index", Guatemala = "Guatemala",
                                      Pakistan = "Pakistan", DRC = "DRC", India = "India"))) + 
  geom_point(alpha = 0.4, position=position_jitterdodge()) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) + 
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(fill = "On Supplement?", y = "", x = "Time (Weeks)",
       title = "Alpha Diversity Comparisons")
ggsave(here("Plots", "Alpha_Diversity_Reviewer_Response.jpg"),
       width = 8,height = 7)
