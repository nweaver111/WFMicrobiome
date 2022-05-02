####### Alpha Diversity Re-Analysis #######

#Why: Need to adjust for ALL covairates by site#


#### Load Libraries ####
library(here)
library(ggplot2)
library(emmeans)
library(geepack)


#### Load Data ####
df <- read.csv(here("Data", "Correct_Data_2022_02_19.csv"), header = T)
covars <- readRDS(here("Data", "Covariate", "brayCurtis_covariates.rds"))

covars_small <- lapply(covars, function(x) gsub("BatchAdj", "std_reads",x))
covars_small <- lapply(covars_small, function(x) unique(x))

#### Data Checks ####
str(df) #all alpha diversity metrics are numeric.
df$Supp_Status <- as.factor(df$Supp_Status)

#### Model Prep ####
sites <- unique(df$SiteLocation.x)
alpha_metrics <- c("Sobs_3700Seqs", "Chao1_3700Seqs", "ShannonH_3700Seqs", "ShannonE_3700Seqs")

tst_obj <- vector('list')

figure_dat <- data.frame(Interaction = NA, Outcome = NA, Timepoint = NA, Site = NA, Supp_Status = NA, 
                         Mean = NA, LCL = NA, UCL = NA)

#### Interaction Models ####
interaction_res <- vector('list')
sensitivity_res <- vector('list')
for(i in sites){
  for(j in alpha_metrics){
    df2 <- df[df$SiteLocation.x == i,]
    df2 <- df2[order(df2$SubjectID, df2$Timepoint),] #Make sure data is ordered by grouping factor!
    
    #Grouping factor needs to be numeric, so do that:
    df2$id1 <- NA
    for(k in 1:length(unique(df2$SubjectID))){
      tmp <- unique(df2$SubjectID)[k]
      df2[df2$SubjectID == tmp, "id1"] <- k
    }
    
    tryCatch({
    interaction_res[[paste(i,"_",j,sep="")]] <- geeglm(as.formula(paste(j, 
                            "~ Timepoint*Supp_Status + ",
                            paste(covars[[i]], collapse = "+"),sep= "")),
           data = df2,
           id = id1)
    
    #sensitivity to outliers:
    mod_tmp <- summary(interaction_res[[paste(i,"_",j,sep="")]])
    sens_resid <- mod_tmp$deviance.resid
    sens_IQR <- quantile(sens_resid, 0.75) - quantile(sens_resid, 0.25)
    sens_cut <- sens_IQR*3
    sens_out <- which(sens_resid < (quantile(sens_resid,0.25) - sens_cut) | sens_resid > (quantile(sens_resid,0.75) + sens_cut))
    
    df_sens <- df2[!(rownames(df2) %in% sens_out),]
    
    sensitivity_res[[paste(i,"_",j,sep="")]] <- geeglm(as.formula(paste(j, 
                            "~ Timepoint*Supp_Status + ",
                            paste(covars[[i]], collapse = "+"),sep= "")),
           data = df_sens,
           id = id1)
    
    #Get adjusted means plots (with CIs) for Time and Supplement for each outcome:
    tst_obj[[paste(i,"_",j,"_interaction",sep="")]] <- emmeans(interaction_res[[paste(i,"_",j,sep="")]]
                       ,specs = c("Timepoint","Supp_Status"))
    
    tst <- summary(tst_obj[[paste(i,"_",j,"_interaction",sep="")]],
                   infer = TRUE,
                   level = 0.95)
    
    
    figure_dat <- rbind(figure_dat, cbind(Interaction = "Yes",
                                          Outcome = j, 
                                          Timepoint = "12Weeks", 
                                          Site = i, 
                                          Supp_Status = 0,
                                          Mean = tst[tst$Timepoint=="12Weeks" & tst$Supp_Status == 0, "emmean"], 
                                          LCL = tst[tst$Timepoint=="12Weeks"& tst$Supp_Status == 0,"asymp.LCL"],
                                          UCL = tst[tst$Timepoint=="12Weeks"& tst$Supp_Status == 0,"asymp.UCL"]))
    
    
    figure_dat <- rbind(figure_dat, cbind(Interaction = "Yes",
                                          Outcome = j, 
                                          Timepoint = "34Weeks", 
                                          Site = i, 
                                          Supp_Status = 0,
                                          Mean = tst[tst$Timepoint=="34Weeks"& tst$Supp_Status == 0, "emmean"], 
                                          LCL = tst[tst$Timepoint=="34Weeks"& tst$Supp_Status == 0,"asymp.LCL"],
                                          UCL = tst[tst$Timepoint=="34Weeks"& tst$Supp_Status == 0,"asymp.UCL"]))
    
    figure_dat <- rbind(figure_dat, cbind(Interaction = "Yes",
                                          Outcome = j, 
                                          Timepoint = "12Weeks", 
                                          Site = i,
                                          Supp_Status = 1,
                                          Mean = tst[tst$Timepoint =="12Weeks" & tst$Supp_Status == 1, "emmean"],
                                          LCL = tst[tst$Timepoint == "12Weeks" & tst$Supp_Status == 1, "asymp.LCL"],
                                          UCL = tst[tst$Timepoint == "12Weeks" & tst$Supp_Status == 1, "asymp.UCL"]))
    
    figure_dat <- rbind(figure_dat, cbind(Interaction = "Yes",
                                          Outcome = j, 
                                          Timepoint = "34Weeks", 
                                          Site = i, 
                                          Supp_Status = 1,
                                          Mean = tst[tst$Timepoint =="34Weeks" & tst$Supp_Status == 1, "emmean"],
                                          LCL = tst[tst$Timepoint == "34Weeks" & tst$Supp_Status == 1, "asymp.LCL"],
                                          UCL = tst[tst$Timepoint == "34Weeks" & tst$Supp_Status == 1, "asymp.UCL"]))
    
    
    }, error=function(e){
      #nothing
    })
  }
}

figure_dat <- figure_dat[-1,]


#### No Interaction Models ####
addition_res <- vector('list')
addition_sensitivity_res <- vector('list')
for(i in sites){
  for(j in alpha_metrics){
    
    tryCatch({
      
      df2 <- df[df$SiteLocation.x == i,]
      df2 <- df2[order(df2$SubjectID, df2$Timepoint),] #Make sure data is ordered by grouping factor!
      
      #Grouping factor needs to be numeric, so do that:
      df2$id1 <- NA
      for(k in 1:length(unique(df2$SubjectID))){
        tmp <- unique(df2$SubjectID)[k]
        df2[df2$SubjectID == tmp, "id1"] <- k
      }
      
      addition_res[[paste(i,"_",j,sep="")]] <- geeglm(as.formula(paste(j, 
                                                                          "~ Timepoint + ",
                                                                          paste(covars[[i]], collapse = "+"),sep= "")),
                                                         data = df2,
                                                         id = id1)
      
      
      #sensitivity to outliers:
      mod_tmp <- summary(addition_res[[paste(i,"_",j,sep="")]])
      sens_resid <- mod_tmp$deviance.resid
      sens_IQR <- quantile(sens_resid, 0.75) - quantile(sens_resid, 0.25)
      sens_cut <- sens_IQR*3
      sens_out <- which(sens_resid < (quantile(sens_resid,0.25) - sens_cut) | sens_resid > (quantile(sens_resid,0.75) + sens_cut))
      
      df_sens <- df2[!(rownames(df2) %in% sens_out),]
      
      addition_sensitivity_res[[paste(i,"_",j,sep="")]] <- geeglm(as.formula(paste(j, 
                                                                          "~ Timepoint + ",
                                                                          paste(covars[[i]], collapse = "+"),sep= "")),
                                                         data = df_sens,
                                                         id = id1)
      
      if(i == "India"){
        
        tst_obj[[paste(i,"_",j,sep="")]] <- emmeans(addition_res[[paste(i,"_",j,sep="")]]
                                                                   ,specs = c("Timepoint","Supp_Status"))
        
        tst <- summary(tst_obj[[paste(i,"_",j,sep="")]],
                       infer = TRUE,
                       level = 0.95)
        
        figure_dat <- rbind(figure_dat, cbind(Interaction = "No",
                                              Outcome = j, 
                                              Timepoint = "12Weeks", 
                                              Site = i, 
                                              Supp_Status = 0,
                                              Mean = tst[tst$Timepoint=="12Weeks" & tst$Supp_Status == 0, "emmean"], 
                                              LCL = tst[tst$Timepoint=="12Weeks"& tst$Supp_Status == 0,"asymp.LCL"],
                                              UCL = tst[tst$Timepoint=="12Weeks"& tst$Supp_Status == 0,"asymp.UCL"]))
        
        
        figure_dat <- rbind(figure_dat, cbind(Interaction = "No",
                                              Outcome = j, 
                                              Timepoint = "12Weeks", 
                                              Site = i,
                                              Supp_Status = 1,
                                              Mean = tst[tst$Timepoint =="12Weeks" & tst$Supp_Status == 1, "emmean"],
                                              LCL = tst[tst$Timepoint == "12Weeks" & tst$Supp_Status == 1, "asymp.LCL"],
                                              UCL = tst[tst$Timepoint == "12Weeks" & tst$Supp_Status == 1, "asymp.UCL"]))
        
        figure_dat <- rbind(figure_dat, cbind(Interaction = "No",
                                              Outcome = j, 
                                              Timepoint = "34Weeks", 
                                              Site = i, 
                                              Supp_Status = 1,
                                              Mean = tst[tst$Timepoint =="34Weeks" & tst$Supp_Status == 1, "emmean"],
                                              LCL = tst[tst$Timepoint == "34Weeks" & tst$Supp_Status == 1, "asymp.LCL"],
                                              UCL = tst[tst$Timepoint == "34Weeks" & tst$Supp_Status == 1, "asymp.UCL"]))
        
      }
    }, error=function(e){
      #nothing
    })
  }
}


#### Sensitivity Checks ####

#loop through all scenarios and compare full data to sensitivity data results and then export table
table_comp <- data.frame(Site = NA, Alpha = NA, Model = NA, Covariate = NA,
                         Estimate = NA, Estimate_p = NA,
                         Sens_Estimate = NA, Sens_Estimate_p = NA,
                         Conclusion = NA)
for(i in sites){
  for(j in alpha_metrics){
    #get coefficients, but India will cause an error for interaction so use tryCatch
    tryCatch({
      inter_coefs <- summary(interaction_res[[paste(i,"_",j,sep="")]])$coefficients
      inter_coefs_sens <- summary(sensitivity_res[[paste(i,"_",j,sep="")]])$coefficients
      
      #Subset to coefficient we need for each and add to table:
      table_comp <- rbind(table_comp,
                          cbind(Site = i, Alpha = j, Model = "Interaction", Covariate = "Interaction",
                                Estimate = inter_coefs["Timepoint34Weeks:Supp_Status1", "Estimate"],
                                Estimate_p = inter_coefs["Timepoint34Weeks:Supp_Status1","Pr(>|W|)"],
                                Sens_Estimate = inter_coefs_sens["Timepoint34Weeks:Supp_Status1", "Estimate"],
                                Sens_Estimate_p = inter_coefs_sens["Timepoint34Weeks:Supp_Status1", "Pr(>|W|)"],
                                Conclusion = NA))

      
    }, error = function(e){
      #Do nothing
    })
    
    add_coefs <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients
    add_coefs_sens <- summary(addition_sensitivity_res[[paste(i,"_",j,sep="")]])$coefficients
    
    table_comp <- rbind(table_comp,
                        cbind(Site = i, Alpha = j, Model = "Additive", Covariate = "Time",
                              Estimate = add_coefs["Timepoint34Weeks", "Estimate"],
                              Estimate_p = add_coefs["Timepoint34Weeks","Pr(>|W|)"],
                              Sens_Estimate = add_coefs_sens["Timepoint34Weeks", "Estimate"],
                              Sens_Estimate_p = add_coefs_sens["Timepoint34Weeks", "Pr(>|W|)"],
                              Conclusion = NA))
    
    table_comp <- rbind(table_comp,
                        cbind(Site = i, Alpha = j, Model = "Additive", Covariate = "Supplement Status",
                              Estimate = add_coefs["Supp_Status1", "Estimate"],
                              Estimate_p = add_coefs["Supp_Status1","Pr(>|W|)"],
                              Sens_Estimate = add_coefs_sens["Supp_Status1", "Estimate"],
                              Sens_Estimate_p = add_coefs_sens["Supp_Status1", "Pr(>|W|)"],
                              Conclusion = NA))
  }
}
View(table_comp)
table_comp <- table_comp[-1,]

write.csv(table_comp, here("Results", "Alpha", "Sensitivity_Check.csv"))
#### Plots ####

#prep:
figure_dat$Outcome <- ifelse(figure_dat$Outcome == "Sobs_3700Seqs",
                             "Sobs",
                             ifelse(figure_dat$Outcome == "Chao1_3700Seqs",
                                    "Chao 1",
                                    ifelse(figure_dat$Outcome == "ShannonE_3700Seqs",
                                           "Shannon H/Hmax", "Shannon H")))

str(figure_dat)
figure_dat[,c("Mean", "LCL", "UCL")] <- apply(figure_dat[,c("Mean", "LCL", "UCL")], 2, as.numeric)

dat_text <- expand.grid(
  Site = sites,
  Index = alpha_metrics
)
dat_text$Label <- NA
dat_text$Interaction <- NA
dat_text$Time <- NA
dat_text$Supplement <- NA
dat_text$Interaction_Estimate <- NA
dat_text$Time_Estimate <- NA
dat_text$Supplement_Estimate <- NA
dat_text$Interaction_Standard_Error <- NA
dat_text$Time_Standard_Error <- NA
dat_text$Supplement_Standard_Error <- NA


for(i in sites){
  for(j in alpha_metrics){
    if(i != "India"){
      dat_text[dat_text$Site == i & dat_text$Index == j, "Interaction_Estimate"] <- summary(interaction_res[[paste(i,"_",j,sep="")]])$coefficients["Timepoint34Weeks:Supp_Status1","Estimate"]
      dat_text[dat_text$Site == i & dat_text$Index == j, "Interaction_Standard_Error"] <- summary(interaction_res[[paste(i,"_",j,sep="")]])$coefficients["Timepoint34Weeks:Supp_Status1","Std.err"]
      dat_text[dat_text$Site == i & dat_text$Index == j, "Interaction"] <- summary(interaction_res[[paste(i,"_",j,sep="")]])$coefficients["Timepoint34Weeks:Supp_Status1","Pr(>|W|)"]
      
    }

    dat_text[dat_text$Site == i & dat_text$Index == j, "Time_Estimate"] <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients["Timepoint34Weeks","Estimate"]
    dat_text[dat_text$Site == i & dat_text$Index == j, "Time_Standard_Error"] <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients["Timepoint34Weeks","Std.err"]
    dat_text[dat_text$Site == i & dat_text$Index == j, "Time"] <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients["Timepoint34Weeks","Pr(>|W|)"]
    dat_text[dat_text$Site == i & dat_text$Index == j, "Supplement_Estimate"] <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients["Supp_Status1","Estimate"]
    dat_text[dat_text$Site == i & dat_text$Index == j, "Supplement_Standard_Error"] <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients["Supp_Status1","Std.err"]
    dat_text[dat_text$Site == i & dat_text$Index == j, "Supplement"] <- summary(addition_res[[paste(i,"_",j,sep="")]])$coefficients["Supp_Status1","Pr(>|W|)"]
    
  }
}

#Adjust p-values for multiple testing:
dat_text$Interaction_FDR <- p.adjust(dat_text$Interaction, 'fdr')
dat_text$Time_FDR <- p.adjust(dat_text$Time, 'fdr')
dat_text$Supplement_FDR <- p.adjust(dat_text$Supplement, 'fdr')

dat_text$Label <- paste("p_time :",round(dat_text$Time_FDR,4)
                        ,"\np_supplement: ", round(dat_text$Supplement_FDR,4)
                        ,sep="")

dat_text$Supp_Status <- as.factor(0)
dat_text$Timepoint <- "12Weeks"
dat_text$Mean <- 3
dat_text$Outcome <- dat_text$Index
dat_text$Outcome <- ifelse(dat_text$Outcome == "Sobs_3700Seqs",
                             "Sobs",
                             ifelse(dat_text$Outcome == "Chao1_3700Seqs",
                                    "Chao 1",
                                    ifelse(dat_text$Outcome == "ShannonE_3700Seqs",
                                           "Shannon H/Hmax", "Shannon H")))

pd = position_dodge(0.5)
ggplot(figure_dat[figure_dat$Outcome != "Sobs",], aes(x = Timepoint, y = Mean, color = Supp_Status)) + theme_bw() +
  geom_point(size = 5, shape = 18, position=pd) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL),width = .4, position=pd) +
  labs(x = "Time (Weeks)",
       y = "",
       color = "On Supplement?") +
  scale_color_discrete(labels=c("0" = "No", "1" = "Yes")) +
  scale_x_discrete(labels=c("12Weeks" = "12",
                            "34Weeks" = "34")) + 
  facet_grid(Outcome ~ Site, scales = "free_y", switch = "y",
             labeller = as_labeller(c(`Chao 1` = "Chao 1 Index", `Shannon H` = "Shannon H Index", 
                                      `Shannon H/Hmax` = "Shannon H/Hmax Index", Guatemala = "Guatemala",
                                      Pakistan = "Pakistan", DRC = "DRC", India = "India"))) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside") + 
  geom_text(
    color = "black",
    size    = 3.5,
    data    = dat_text[dat_text$Outcome != "Sobs",],
    mapping = aes(x = Inf, y = Inf, label = Label),
    vjust = 1.02,
    hjust = 1.02
  )



ggsave(here("Plots", "Alpha_Diversity_Figure_1.jpg"),
       width = 10,height = 7)
write.csv(dat_text, here("Results", "Alpha", "All_Alpha_Models_Updated.csv"))
