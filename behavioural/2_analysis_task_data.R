#####################################################################################
########################### script analysis task paper ##############################
#####################################################################################

options(scipen=999)

# empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
analysis_dir <- getwd()
fig_dir <- file.path(analysis_dir, "figures")
if (dir.exists(fig_dir) == F ) {
  dir.create(fig_dir)
}

filename_tables <- "tables_task_paper.xlsx"

# delete output files
ifelse(file.exists(filename_tables), file.remove(filename_tables), FALSE)

# helper functions and packages #
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/errorbars.R?raw=TRUE")
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/rbindcolumns.R?raw=TRUE")

library(psych)
library(dplyr)
library(lme4)
library(ggplot2)
library(gt)

library(metafor)


# define version 
versions <- c("pilot", "main", "fmri")
versions_redo <- c("Behavioural study", "Replication", "fMRI study")

# read in data sets
dfWide_all <- read.csv(paste0("wide_MagicBehavioural.csv"), stringsAsFactors = F)
dfLong_all <- read.csv(paste0("long_MagicBehavioural.csv"), stringsAsFactors = F)

# dataLong_all has two NAs in the recall data
dfLong_all$cuedRecallLenient <- ifelse(is.na(dfLong_all$cuedRecallLenient), 0, dfLong_all$cuedRecallLenient)
dfLong_all$cuedRecallStrict <- ifelse(is.na(dfLong_all$cuedRecallStrict), 0, dfLong_all$cuedRecallStrict)

#### CREATE DEMOGRAPHICS TABLE FOR PAPER (i.e. Table 1) ####
demogs <- data.frame() # 1. col: variable, 2. col = control, 3. col = experimental

col = 0
for (version in versions){
  
  i = 0 # row iterator
  col = col + 2 # col iterator
  varCol <- 1
  contCol <- col
  expCol <- contCol + 1
  
  # subject data
  dfWide <- subset(dfWide_all, dfWide_all$sample == version)
  
  # n per group
  N_per_group <- plyr::count(dfWide, vars = "group")
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Subjects per group"
  }
  demogs[i,contCol] <- paste0("n = ", N_per_group$freq[N_per_group$group == "cont"]) # control group
  demogs[i,expCol] <- paste0("n = ", N_per_group$freq[N_per_group$group == "exp"]) # incentive group
  
  # age
  psych::describe(dfWide$age)
  age <- psych::describeBy(dfWide[,"age"], group=dfWide$group)
  age <- as.data.frame(rbind(age$cont, age$exp))
  age$mot <- rep(c("cont","exp"), each = 1)
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Age"
  }
  demogs[i,contCol] <- paste(paste0(format(round(age$mean[1],2), nsmall = 2), " (", format(round(age$sd[1], 2), nsmall = 2), ")"), paste0("[",age$min[1],"; ", age$max[1],"]"), sep = '\n')  # control group
  demogs[i, expCol] <- paste(paste0(format(round(age$mean[2],2), nsmall = 2), " (", format(round(age$sd[2], 2), nsmall = 2), ")"), paste0("[",age$min[2],"; ", age$max[2],"]"), sep = '\n')  # incentive group
  
  # gender
  plyr::count(dfWide, vars =  "gender")
  gender <- plyr::count(dfWide, vars = c("gender","group"))
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Gender (% female)"
  }
  demogs[i,contCol] <- paste0(format(round(gender$freq[gender$gender == "Female" & gender$group == "cont"]/N_per_group$freq[N_per_group$group == "cont"]*100,2), nsmall = 2)) # control group
  demogs[i, expCol] <- paste0(format(round(gender$freq[gender$gender == "Female" & gender$group == "exp"]/N_per_group$freq[N_per_group$group == "exp"]*100,2), nsmall = 2)) # experimental group
  
  # ethnicity
  ethnicity <- plyr::count(dfWide, vars = c("ethnicity","group"))
  dfWide$ethnicity_BAME <- ifelse(dfWide$ethnicity == "White British" | dfWide$ethnicity == "Other White", "White", "BAME")
  BAME <- plyr::count(dfWide, vars = c("ethnicity_BAME","group"))
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Ethnicity (% BAME)"
  }
  demogs[i,contCol] <- paste0(format(round(BAME$freq[BAME$ethnicity_BAME == "BAME" & BAME$group == "cont"]/N_per_group$freq[N_per_group$group == "cont"]*100,2), nsmall = 2)) # control group
  demogs[i, expCol] <- paste0(format(round(BAME$freq[BAME$ethnicity_BAME == "BAME" & BAME$group == "exp"]/N_per_group$freq[N_per_group$group == "exp"]*100,2), nsmall = 2)) # experimental group
  
  # education
  dfWide$yearsOfEducation_num <- ifelse(dfWide$yearsOfEducation == "17 years", 17,
                                        ifelse(dfWide$yearsOfEducation == "20 (18.5 excluding time as phd student)", 20, 
                                               as.numeric(as.character(dfWide$yearsOfEducation))))
  
  psych::describe(dfWide$yearsOfEducation_num)
  edu <- psych::describeBy(dfWide[,"yearsOfEducation_num"], group=dfWide$group)
  edu <- as.data.frame(rbind(edu$cont, edu$exp))
  edu$mot <- rep(c("cont","exp"), each = 1)
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Years of Education"
  }
  demogs[i,contCol] <- paste(paste0(format(round(edu$mean[1],2), nsmall = 2), " (", format(round(edu$sd[1], 2), nsmall = 2), ")"), paste0("[",edu$min[1],"; ", edu$max[1],"]"), sep = '\n')  # control group
  demogs[i, expCol] <- paste(paste0(format(round(edu$mean[2],2), nsmall = 2), " (", format(round(edu$sd[2], 2), nsmall = 2), ")"), paste0("[",edu$min[2],"; ", edu$max[2],"]"), sep = '\n')  # incentive group
  
  # days between incidental encoding and memory test
  psych::describe(dfWide$daysBetweenExpAndMemory)
  days <- psych::describeBy(dfWide[,"daysBetweenExpAndMemory"], group=dfWide$group)
  days <- as.data.frame(rbind(days$cont, days$exp))
  days$mot <- rep(c("cont","exp"), each = 1)
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Days btw sessions"
  }
  demogs[i,contCol] <- paste(paste0(format(round(days$mean[1],2), nsmall = 2), " (", format(round(days$sd[1],2), nsmall = 2), ")"), paste0("[",format(round(days$min[1],2), nsmall = 2),"; ", format(round(days$max[1],2), nsmall = 2),"]"), sep = '\n')  # control group
  demogs[i, expCol] <- paste(paste0(format(round(days$mean[2],2), nsmall = 2), " (", format(round(days$sd[2],2), nsmall = 2), ")"), paste0("[",format(round(days$min[2],2), nsmall = 2),"; ", format(round(days$max[2],2), nsmall = 2),"]"), sep = '\n')  # incentive group
  
  # experience with magic tricks
  dfWide$magictrickExperience_num <- ifelse(dfWide$magictrickExperience == " Very frequently", 6,
                                            ifelse(dfWide$magictrickExperience == " Frequently", 5,
                                                   ifelse(dfWide$magictrickExperience == " Occasionally ", 4,
                                                          ifelse(dfWide$magictrickExperience == " Rarely ", 3,
                                                                 ifelse(dfWide$magictrickExperience == " Very rarely ", 2,
                                                                        ifelse(dfWide$magictrickExperience == " Never", 1, NA))))))
  
  magic <- psych::describeBy(dfWide[,"magictrickExperience_num"], group=dfWide$group)
  magic <- as.data.frame(rbind(magic$cont, magic$exp))
  magic$mot <- rep(c("cont","exp"), each = 1)
  
  # add to table
  i = i+1
  if (version == "pilot"){
    demogs[i,varCol] <- "Experience with magic"
  }
  demogs[i,contCol] <- paste(paste0(format(round(magic$mean[1],2), nsmall = 2), " (", format(round(magic$sd[1],2), nsmall = 2), ")"), paste0("[",format(round(magic$min[1],2), nsmall = 2),"; ", format(round(magic$max[1],2), nsmall = 2),"]"), sep = '\n')  # control group
  demogs[i, expCol] <- paste(paste0(format(round(magic$mean[2],2), nsmall = 2), " (", format(round(magic$sd[2],2), nsmall = 2), ")"), paste0("[",format(round(magic$min[2],2), nsmall = 2),"; ", format(round(magic$max[2],2), nsmall = 2),"]"), sep = '\n')  # incentive group
}

# add two header rows
demogs_write <- rbind(c("", "Behavioural study", "", "Replication", "", "fMRI study", ""), 
                      #                      c("", "Control group", "Incentive group", "Control group", "Incentive group", "Control group", "Incentive group"),
                      c("", paste("Control", "group", sep = "\n"), paste("Incentive", "group", sep = "\n"), paste("Control", "group", sep = "\n"), paste("Incentive", "group", sep = "\n"), paste("Control", "group", sep = "\n"), paste("Incentive", "group", sep = "\n")),
                      demogs)

# remove all objects no longer needed
rm(age, BAME, days, edu, ethnicity, gender, N_per_group, magic)

# save demogs
xlsx::write.xlsx(demogs_write, file=filename_tables, sheetName = "Table_1", append = T, row.names = F, col.names = F)

### IMI scores

#### Compute t-tests and effect sizes for between-group differences in IMI questionnaire scores #### 
scales <-c("IMI_intrinsicMotivation", "IMI_taskEngagement", "IMI_interest", "IMI_boredom", "IMI_effort", "IMI_pressure")
scales_new <-c("Intrinsic motivation", "Task engagement", "Interest", "Boredom", "Effort", "Pressure")

imi <- data.frame() # 1. col: variable, 2. col = control, 3. col = experimental, 4. col = group comparison, 5. col = cohen's d
imi_2 <- data.frame() # 1. col: descriptor, col 2 - 7 = imi scales
j = 0
for (v in seq_along(versions)){
  
  i = 0 # row iterator
  
  # subset data
  dfWide <- subset(dfWide_all, dfWide_all$sample == versions[v])
  
  # determine columns
  if (versions[v] == "pilot"){
    varCol <- 1
    contCol <- varCol + 1
  } else {
    contCol <- cohenCol + 1
  }
  expCol <- contCol + 1
  testCol <- expCol + 1
  cohenCol <- testCol + 1
  
  for(scale in 1:length(scales)) {
    print(scales[scale])
    
    # compute mean and SD for within each group
    group_stat <- psych::describeBy(dfWide[,scales[scale]], group=dfWide$group)
    group_stat <- as.data.frame(rbind(group_stat$cont, group_stat$exp))
    group_stat$mot <- rep(c("cont","exp"), each = 1)
    
    # check whether data is normally distributed within each group
    cont <- shapiro.test(dfWide[dfWide$group == "cont",scales[scale]])
    exp <- shapiro.test(dfWide[dfWide$group == "exp",scales[scale]])
    
    # checl whether shapiro wilk is significant in either of the groups
    #if (cont$p.value >= 0.05 && exp$p.value >= 0.05){
    # compute two sample t-test
    groupcomp <- t.test(dfWide[,scales[scale]]~dfWide$group)
    name_test <- "Two Sample t-Test"
    test_output <- paste("t(",format(round(groupcomp$parameter, 3), nsmall = 3), ") = ", format(round(groupcomp$statistic, 3), nsmall = 3), ", p = ", format(round(groupcomp$p.value, 3), nsmall = 3))
    test_output <- paste(paste0("t(",format(round(groupcomp$parameter, 3), nsmall = 3), ") = ", format(round(groupcomp$statistic, 3), nsmall = 3), ","), paste0("p = ", format(round(groupcomp$p.value, 3), nsmall = 3)), sep = "\n")
    #} else {
    #  # compute wilcox
    #  groupcomp <- wilcox.test(dfWide[,scales[scale]]~dfWide$group, exact = FALSE) 
    #  name_test <- "Wilcoxon Rank Sum Test"
    #  test_output <- paste0("W = ", round(groupcomp$statistic, 3), nsmall = 3), ", p = ", round(groupcomp$p.value, 3), nsmall = 3))
    #}
    rm(cont, exp)
    
    # compute effect size
    if (group_stat$mean[1] != group_stat$mean[2]) {
      data <- dfWide[,c("group", scales[scale])]
      psych::cohen.d(data, "group")
      d <- psych::cohen.d(data, "group")
      cohen_output <- paste0(format(round(d$cohen.d[2], 2), nsmall = 2), " [", format(round(d$cohen.d[1], 2), nsmall = 2), "; ", format(round(d$cohen.d[3], 2), nsmall = 2),"]")
    }
    
    
    # put it all together into a table (imi)
    i = i+1
    if (versions[v] == "pilot"){
      imi[i,varCol] <- paste(scales_new[scale])
    }
    imi[i,contCol] <- paste0(format(round(group_stat$mean[1],2), nsmall = 2), " (", format(round(group_stat$sd[1], 2), nsmall = 2), ") [",format(round(group_stat$min[1], 2), nsmall = 2),"; ", format(round(group_stat$max[1], 2), nsmall = 2),"]") # control
    imi[i,expCol] <- paste0(format(round(group_stat$mean[2],2), nsmall = 2), " (", format(round(group_stat$sd[2], 2), nsmall = 2), ") [",format(round(group_stat$min[2], 2), nsmall = 2),"; ", format(round(group_stat$max[2], 2), nsmall = 2),"]")
    imi[i,testCol] <- paste(test_output)
    imi[i,cohenCol] <- paste(cohen_output)
    
    # put it all together into a table (imi_2)
    imi_2[j+1,varCol] <- paste0(versions_redo[v])
    imi_2[j+2,varCol] <- "Control group"
    imi_2[j+3,varCol] <- "Incentive group"
    imi_2[j+4,varCol] <- "Group comparison"
    imi_2[j+5,varCol] <- "Cohen's d"
    scaleCol <- scale + 1
    imi_2[j+1,scaleCol] <- ""
    imi_2[j+2,scaleCol] <- paste0(format(round(group_stat$mean[1],2), nsmall = 2), " (", format(round(group_stat$sd[1], 2), nsmall = 2), ") [",format(round(group_stat$min[1], 2), nsmall = 2),"; ", format(round(group_stat$max[1], 2), nsmall = 2),"]")
    #imi_2[j+2,scaleCol] <- paste(paste0(format(round(group_stat$mean[1],2), nsmall = 2), " (", format(round(group_stat$sd[1], 2), nsmall = 2), ")"), paste0("[",format(round(group_stat$min[1], 2), nsmall = 2),"; ", format(round(group_stat$max[1], 2), nsmall = 2),"]"), sep = '\n')
    
    imi_2[j+3,scaleCol] <- paste0(format(round(group_stat$mean[2],2), nsmall = 2), " (", format(round(group_stat$sd[2], 2), nsmall = 2), ") [",format(round(group_stat$min[2], 2), nsmall = 2),"; ", format(round(group_stat$max[2], 2), nsmall = 2),"]")
    #imi_2[j+3,scaleCol] <- paste(paste0(format(round(group_stat$mean[2],2), nsmall = 2), " (", format(round(group_stat$sd[2], 2), nsmall = 2), ")"), paste0("[",format(round(group_stat$min[2], 2), nsmall = 2),"; ", format(round(group_stat$max[2], 2), nsmall = 2),"]"), sep = '\n')
    
    imi_2[j+4,scaleCol] <- paste(test_output)
    imi_2[j+5,scaleCol] <- paste(cohen_output)
    
  }
  j <- j + 5
}

#names(imi) <- c("", "Control group", "Incentive group", "group comparison", "Cohen's d [CI_l; CI_u]", "Control group", "Incentive group", "group comparison", "Cohen's d [CI_l; CI_u]", "Control group", "Incentive group", "Group comparison", "Cohen's d [CI_l; CI_u]")
names(imi_2) <- c("", scales_new)

# add two header rows
imi_write <- rbind(c("", "Behavioural study", "", "", "", "Replication", "", "", "", "fMRI study", "", "", ""), 
                   c("", "Control group", "Incentive group", "Group comparison", "Cohen's d", "Control group", "Incentive group", "Group comparison", "Cohen's d", "Control group", "Incentive group", "Group comparison", "Cohen's d"),
                   imi)

# save imi
#xlsx::write.xlsx(imi_write, file=filename_tables, sheetName = "Table_S1", append = T, row.names = F, col.names = F) # note: row.names contain variables
xlsx::write.xlsx(imi_2, file=filename_tables, sheetName = "Table_S1", append = T, row.names = F) # note: row.names contain variables


#### Behaviour during the task ####

dependentVariables <- c("responseCuriosity") # note: this is the raw value, not mean-centered
varNames <- c("Curiosity ratings")

for (v in seq_along(versions)){
  
  # subset data
  dfLong <- subset(dfLong_all, dfLong_all$sample == versions[v])
  
  #  make sure that ID and stimID are factors    
  dfLong$ID <- as.factor(as.character(dfLong$ID)) 
  dfLong$stimID <- as.factor(as.character(dfLong$stimID))
  
  for (DV in 1:length(dependentVariables)){
    
    # compute LME
    LME <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + (1|ID) + (1|stimID), data = dfLong)
    sum_LME <- summary(LME)
    print(paste(versions[v], dependentVariables[DV]))
    LME_coef <- as.data.frame(round(sum_LME$coefficients, 3))
    
    LME_coef$param <- c("Intercept", "Incentive effect")
    
    LME_coef <- rbind(c("", "", "", versions_redo[v]), LME_coef)
    
    #row.names(LME_coef) <- c(paste("LME",version, "Intercept"),
    #                         paste("LME",version, "group"))
    
    if(versions[v] == "pilot"){
      LME_cur <- LME_coef
    } else {
      LME_cur <- rbind(LME_cur, LME_coef)
    }
    rm(LME_coef, sum_LME)
    
  }
  
}

# re-order column
LME_cur <- LME_cur[, c("param", "Estimate", "Std. Error", "t value")]
names(LME_cur) <- c("", "Estimate", "SE", "t value")

# save imi
xlsx::write.xlsx(LME_cur, file=filename_tables, sheetName = "Table_S2", append = T, row.names = F)


#### Effects of curiosity, incentive and their interaction on memory ####

dv_recog <- c("recognition", "recognitionConfLevel_above_1", "recognitionConfLevel_above_2", "recognitionConfLevel_above_3", 
              "recognitionConfLevel_above_4", "recognitionConfLevel_above_5", 
              "cuedRecallStrict")

dv_recog_new <- c("Recog [Conf > 0]", "Recog [Conf > 1]", "Recog [Conf > 2]", "Recog [Conf > 3]", "Recog [Conf > 4]", "Recog [Conf > 5]", 
                  "Recall")

main_dvs <- c("Recognition", "High confidence recognition", 
              "Cued recall")
effects <- c("Monetary incentive", "Curiosity", "Interaction")
effects <- c("Curiosity", "Monetary incentive", "Interaction")

models <- c("full", "reduced")

# go into output directory for figures
setwd(fig_dir)

for (model in models){
  
  # create empty data frames
  output_temp <- data.frame()
  
  for (version in versions){
    
    # subset data
    dfLong <- subset(dfLong_all, dfLong_all$sample == version)
    
    #  make sure that ID and stimID are factors    
    dfLong$ID <- as.factor(as.character(dfLong$ID)) 
    dfLong$stimID <- as.factor(as.character(dfLong$stimID)) 
    
    # compute the main model of interest and save the output in a csv file
    for(DV in seq_along(dv_recog)){
      
      DVV <- length(dv_recog)+DV
      DVVV <- 2*length(dv_recog)+DV
      
      # compute LME model
      if (model == "full"){
        gLME <- glmer(dfLong[, dv_recog[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dfLong)
      } else { 
        gLME <- glmer(dfLong[, dv_recog[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dfLong)
      }
      
      # extract coefficients: effect estimate
      curiosityEffect <- fixef(gLME)[3] # curiosity effect
      incentiveEffect <- fixef(gLME)[2] # incentive effect
      interactionEffect <- fixef(gLME)[4] # curiosity effect
      
      # extract coefficients: standard error
      SEcuriosityEffect <- sqrt(vcov(gLME)[3, 3]) # curiosity SE
      SEincentiveEffect <- sqrt(vcov(gLME)[2, 2]) # incentive SE
      SEinteractionEffect <- sqrt(vcov(gLME)[4, 4]) # interaction SE
      
      # compute 95% CI for each effect
      curiosityLower <- curiosityEffect - 1.96*SEcuriosityEffect
      curiosityUpper <- curiosityEffect + 1.96*SEcuriosityEffect
      incentiveLower <- incentiveEffect - 1.96*SEincentiveEffect
      incentiveUpper <- incentiveEffect + 1.96*SEincentiveEffect
      interactionLower <- interactionEffect - 1.96*SEinteractionEffect
      interactionUpper <- interactionEffect + 1.96*SEinteractionEffect

      # insert into table #
      
      # data collection
      output_temp[DV,1] <- version
      output_temp[DVV,1] <- version
      output_temp[DVVV,1] <- version

      # effect
      output_temp[DV,2] <- "Curiosity" #effect
      output_temp[DVV,2] <- "Monetary incentive" #effect
      output_temp[DVVV,2] <- "Interaction" #effect
      
      # dependent variavble
      output_temp[DV,3] <- paste(dv_recog[DV]) # dependent var - curiosity effect
      output_temp[DVV,3] <- paste(dv_recog[DV]) # dependent var - incentive effect
      output_temp[DVVV,3] <- paste(dv_recog[DV]) # # dependent var - interaction effect
      
      # effect size
      output_temp[DV,4] <- curiosityEffect # beta
      output_temp[DVV,4] <- incentiveEffect # beta
      output_temp[DVVV,4] <- interactionEffect # beta
      
      # se effect
      output_temp[DV,5] <- SEcuriosityEffect # se
      output_temp[DVV,5] <- SEincentiveEffect # se
      output_temp[DVVV,5] <- SEinteractionEffect # se
      
      # lower ci
      output_temp[DV,6] <- curiosityLower # ci
      output_temp[DVV,6] <- incentiveLower # ci
      output_temp[DVVV,6] <- interactionLower # ci
      
      # upper ci
      output_temp[DV,7] <- curiosityUpper # ci
      output_temp[DVV,7] <- incentiveUpper # ci
      output_temp[DVVV,7] <- interactionUpper # ci

      # model fit
      output_temp[DV,8] <- isSingular(gLME)
      output_temp[DVV,8] <- isSingular(gLME)
      output_temp[DVVV,8] <- isSingular(gLME)
      
    }
    
    if (version == "pilot"){
      output <- output_temp
    }
    else {
      output <- rbind.all.columns(output, output_temp)
    }
  }
  
  # rename columns
  names(output) <- c("sample", "effect", "memory_level", "effect_estimate", "se_estimate", "ci_lower", "ci_upper", "singular_fit")
  
  # META ANALYSIS FIXED EFFECTS
  for(DV in seq_along(dv_recog)){
    
    # compute meta
    meta_c <- rma(yi = output$effect_estimate[output$effect=="Curiosity" & output$memory_level==dv_recog[DV]], sei =  output$se_estimate[output$effect=="Curiosity" & output$memory_level==dv_recog[DV]], method="FE")

    meta_r <- rma(yi = output$effect_estimate[output$effect=="Monetary incentive" & output$memory_level==dv_recog[DV]], sei =  output$se_estimate[output$effect=="Monetary incentive" & output$memory_level==dv_recog[DV]], method="FE")
    
    meta_i <- rma(yi = output$effect_estimate[output$effect=="Interaction" & output$memory_level==dv_recog[DV]], sei =  output$se_estimate[output$effect=="Interaction" & output$memory_level==dv_recog[DV]], method="FE")
    
    # collect information of meta analysis in one table --> results table
    meta_temp <- data.frame()
    meta_temp[1:length(effects), 1] <- paste0(dv_recog_new[DV]) # variable name dependent variable
    
    for (e in seq_along(effects)){
      
      # define correct output
      if(effects[e] == "Monetary incentive"){
        meta <- meta_r
      } else if(effects[e] == "Curiosity"){
        meta <- meta_c
      } else {
        meta <- meta_i
      }
      
      # add information for each effect into table
      meta_temp[e, 2] <- paste0(effects[e]) # effect
      meta_temp[e, 3] <- paste0(format(round(meta$b, 3), nsmall = 3), " (", format(round(meta$se, 3), nsmall = 3),")")  # b (se)
      meta_temp[e, 4] <- paste0(format(round(exp(meta$b), 2), nsmall = 2), " [", format(round(exp(meta$ci.lb), 2), nsmall = 2),"; ", format(round(exp(meta$ci.ub), 2), nsmall = 2),"]") # OR [lower ci; upper ci]
      
      meta_temp[e, 5] <- format(round(meta$zval, 3), nsmall = 3) # z
      meta_temp[e, 6] <- format(round(meta$pval, 3), nsmall = 3) # p
      meta_temp[e, 6] <- ifelse(meta_temp[e, 6] == "0.000", "< 0.001", meta_temp[e, 6])
      #meta_temp[e, 8] <- paste0("Q(df = ", meta$k -1,") = ", round(meta$QE, digits = 3), ", p = ", round(meta$QEp, digits = 3)) # q statistics
      #meta_temp[e, 9] <- paste0(round(meta$I2, digits = 3), "%") # I^2 (total heterogeneity / total variability)
      
      rm(meta)
      #names(meta_temp) <- c("dependent variable", "Effect", "b", "(SE)", "z", "p", "95%-CI", "Test for Heterogeneity", "I^2")
      names(meta_temp) <- c("dependent variable", "Effect", "b (SE)", " OR [95%-CI]", "z value", "p value")
      
    }
    
    # combine for all DVs
    if(DV == 1){
      meta_out <- meta_temp
    } else {
      meta_out <- rbind(meta_out, meta_temp)
    }
    rm(meta_temp)
    
    # add results from meta analyses to output table --> plotting
    output[(1+3*length(dv_recog)*length(versions)):(3*length(dv_recog)*(length(versions)+1)),1] <- "FE_meta" # sample
    output[(1+3*length(dv_recog)*length(versions)):(2*length(dv_recog)*(length(versions)+2)),2] <- "Curiosity" # effect
    output[(1+2*length(dv_recog)*(length(versions)+2)):(2*length(dv_recog)*(length(versions)+2)+length(dv_recog)),2] <- "Monetary incentive" # effect
    output[(1+2*length(dv_recog)*(length(versions)+2)+length(dv_recog)):(3*length(dv_recog)*(length(versions)+1)),2] <- "Interaction" # effect
    
    output[(DV+3*length(dv_recog)*length(versions)),3] <- paste(dv_recog[DV]) # memory level
    output[(DV+2*length(dv_recog)*(length(versions)+2)),3] <- paste(dv_recog[DV]) # memory level
    output[(DV+2*length(dv_recog)*(length(versions)+2)+length(dv_recog)),3] <- paste(dv_recog[DV]) # memory level
    
    output[(DV+3*length(dv_recog)*length(versions)),4] <- meta_c$beta # meta effect curiosity
    output[(DV+3*length(dv_recog)*length(versions)),5] <- meta_c$se # meta se curiosity
    
    output[(DV+3*length(dv_recog)*length(versions)),6] <- meta_c$ci.lb # ci curiosity
    output[(DV+3*length(dv_recog)*length(versions)),7] <- meta_c$ci.ub # ci curiosity
    
    output[(DV+2*length(dv_recog)*(length(versions)+2)),4] <- meta_r$beta # meta effect incentive
    output[(DV+2*length(dv_recog)*(length(versions)+2)),5] <- meta_r$se # meta se incentive
    
    output[(DV+2*length(dv_recog)*(length(versions)+2)),6] <- meta_r$ci.lb # ci incentive
    output[(DV+2*length(dv_recog)*(length(versions)+2)),7] <- meta_r$ci.ub # ci incentive
    
    output[(DV+2*length(dv_recog)*(length(versions)+2)+length(dv_recog)),4] <- meta_i$beta # meta effect interaction
    output[(DV+2*length(dv_recog)*(length(versions)+2)+length(dv_recog)),5] <- meta_i$se # meta se interaction
    
    output[(DV+2*length(dv_recog)*(length(versions)+2)+length(dv_recog)),6] <- meta_i$ci.lb # ci interaction
    output[(DV+2*length(dv_recog)*(length(versions)+2)+length(dv_recog)),7] <- meta_i$ci.ub # ci interaction
    
    rm(meta_c, meta_r, meta_i)
    
    # change memory level
    output$memory_level[output$memory_level == dv_recog[DV]] <- paste0(dv_recog_new[DV])
    
  }
  
  # create table with the results of the invidual gLMEs for supplementary material
  lme_sample <- data.frame()
  multi <- 0
  
  for(dv in seq_along(dv_recog_new)){
    
    # add effects
    for (e in seq_along(effects)){
      
      # add dv
      lme_sample[(multi*length(effects))+e,1] <- paste(dv_recog_new[dv])
      
      # add name of effect
      lme_sample[(multi*length(effects))+e,2] <- paste(effects[e])
      
      for (v in seq_along(versions)){
        
        version <- versions[v]
        
        # extract effect size from each version
        effect <-  output$effect_estimate[output$effect==effects[e] & output$memory_level==dv_recog_new[dv] & output$sample==version]
        se <-  output$se_estimate[output$effect==effects[e] & output$memory_level==dv_recog_new[dv] & output$sample==version]
        OR <- format(round(exp(effect), 2), nsmall = 2) # compute OR
        ci_l <- output$ci_lower[output$effect==effects[e] & output$memory_level==dv_recog_new[dv] & output$sample==version]
        ci_l <- format(round(exp(ci_l), 2), nsmall = 2) # overwrite after transforming to OR
        ci_u <- output$ci_upper[output$effect==effects[e] & output$memory_level==dv_recog_new[dv] & output$sample==version]
        ci_u <- format(round(exp(ci_u), 2), nsmall = 2) # overwrite after transforming to OR
        sing_fit <- output$singular_fit[output$effect==effects[e] & output$memory_level==dv_recog_new[dv] & output$sample==version]
        
        # add to table
        if (sing_fit == T) {
          lme_sample[(multi*length(effects))+e,v*length(effects)] <- paste0("* ",format(round(effect,3), nsmall = 3), " (",format(round(se,3), nsmall = 3),")") # beta (se)
        } else {
          lme_sample[(multi*length(effects))+e,v*length(effects)] <- paste0("  ",format(round(effect,3), nsmall = 3), " (",format(round(se,3), nsmall = 3),")") # beta (se)
        }
        lme_sample[(multi*length(effects))+e,v*length(effects)+1] <- OR # OR
        lme_sample[(multi*length(effects))+e,v*length(effects)+2] <- paste0("[",ci_l,"; ", ci_u,"]" ) # [95%-CI]
        
      }
      
    }
    
    multi <- multi + 1
  }
  
  
  # make a more compromised table with the results of the meta analysis
  for (dv in seq_along(dv_recog_new)){
    
    if (dv == 1){
      lme_sample_new <- rbind(c(paste0(dv_recog_new[dv]), paste0(dv_recog_new[dv]), "", "", "", "", "", "", "", "", ""), lme_sample[lme_sample$V1 == paste0(dv_recog_new[dv]),])
    } else {
      lme_sample_new <- rbind(lme_sample_new, c(paste0(dv_recog_new[dv]), paste0(dv_recog_new[dv]), "", "", "", "", "", "", "", "", ""), lme_sample[lme_sample$V1 == paste0(dv_recog_new[dv]),])
    }
  }
  
  # remove dv column
  lme_sample_new$V1 <- NULL
  
  # add two header rows
  lme_sample_new <- rbind(c("", "Behavioural study", "", "", "Replication", "", "", "fMRI study", "", ""), 
                     c("Fixed Effect", "b (SE)", "OR", "95%-CI", "b (SE)", "OR", "95%-CI", "b (SE)", "OR", "95%-CI"),
                     lme_sample_new)
  
  # CREATE TABLE SHOWING INTEGRATED RESULTS MAIN DVS
  
  # subset results
  meta_main_dvs <- subset(meta_out, meta_out$`dependent variable` == "Recog [Conf > 0]" | meta_out$`dependent variable` == "Recog [Conf > 3]" | meta_out$`dependent variable` == "Recall" )
  # change variable names
  meta_main_dvs$`dependent variable` <- ifelse(meta_main_dvs$`dependent variable` == "Recog [Conf > 0]", main_dvs[1],
                                               ifelse(meta_main_dvs$`dependent variable` == "Recog [Conf > 3]", main_dvs[2],
                                                      ifelse(meta_main_dvs$`dependent variable` == "Recall", main_dvs[3],NA)))
  # convert effects to factors to re-order data frame
  meta_main_dvs$Effect <- factor(meta_main_dvs$Effect, levels = effects)
  meta_main_dvs <- meta_main_dvs[order(meta_main_dvs$Effect),]
  
  # add column spanners and remove effect column
  meta_main <- rbind(c(paste0(effects[1]), "", "", "", "", ""), meta_main_dvs[meta_main_dvs$Effect == paste0(effects[1]),],
                     c(paste0(effects[2]), "", "", "", "", ""), meta_main_dvs[meta_main_dvs$Effect == paste0(effects[2]),],
                     c(paste0(effects[3]), "", "", "", "", ""), meta_main_dvs[meta_main_dvs$Effect == paste0(effects[3]),])
  meta_main$Effect <- NULL
  names(meta_main)[names(meta_main) == "dependent variable"] <- ""
  
  # make a more compromised table with the results of the meta analysis
  for (dv in seq_along(dv_recog_new)){
    
    if (dv == 1){
      meta_table <- rbind(c(paste0(dv_recog_new[dv]), paste0(dv_recog_new[dv]), "", "", "", ""), meta_out[meta_out$`dependent variable` == paste0(dv_recog_new[dv]),])
    } else {
    meta_table <- rbind(meta_table, c(paste0(dv_recog_new[dv]), paste0(dv_recog_new[dv]), "", "", "", ""), meta_out[meta_out$`dependent variable` == paste0(dv_recog_new[dv]),])
    }
  }
  
  # subset results
  meta_recog <- subset(meta_table, meta_table$`dependent variable` != "Recall")
  meta_recog$`dependent variable` <- NULL
  
  # remove dependent variable column
  meta_table$`dependent variable` <- NULL
  
  # CREATE BAR GRAPH PLOTTING EFFECTS ON RECOGNITION AS A FUNCTION OF CONFIDENCE

  # prepare data frame for plotting
  output_temp <- subset(output, output$memory_level == "Recog [Conf > 0]" | output$memory_level == "Recog [Conf > 1]" | output$memory_level == "Recog [Conf > 2]" | output$memory_level == "Recog [Conf > 3]" | output$memory_level == "Recog [Conf > 4]" | output$memory_level == "Recog [Conf > 5]" )
  output_temp$sample <- ifelse(output_temp$sample == "fmri", "fMRI study",
                               ifelse(output_temp$sample == "pilot", "Behavioural study",
                                      ifelse(output_temp$sample == "main", "Replication",
                                             ifelse(output_temp$sample == "FE_meta", "FE Meta-analysis", NA))))
  output_temp$sample <- factor(output_temp$sample, levels = c("FE Meta-analysis", "Behavioural study", "Replication", "fMRI study"))

  labels <- c("Conf > 0","Conf > 1","Conf > 2","Conf > 3","Conf > 4","Conf > 5")

  output_temp$cutoff <- rep(0:5)


  # full graph: all four point estimates
  graph_all <- ggplot(output_temp, aes(x=cutoff, y=effect_estimate, colour = sample))+
    geom_hline(yintercept=0, linetype="dashed", size = 0.2, color="grey") +
    stat_smooth(method ='lm', se = F, size = 0.5, position=position_dodge(width=0.5)) +
    geom_point( size = 2, aes(shape = sample, colour = sample), position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.1, position=position_dodge(width=0.5)) +
    theme_classic() +
    #theme(strip.text.y = element_text(size = 8)) +
    scale_x_continuous(breaks = c(0:5), labels = labels) +
    scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "Odds Ratio")) +
    facet_grid(factor(effect, levels=c('Curiosity','Monetary incentive','Interaction')) ~ ., scale = "free", switch = "y") +
    labs(x="Behavioural recognition memory performance", y=paste(model, "gLME FE beta estimate"), color = "Data collection", title = paste("Fixed effects as a function of confidence cut-off")) +
    #theme(axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), title=element_text(size =14, face="bold"), legend.title = element_text(size=12), legend.text = element_text(size =12), strip.text = element_text(size =12)) +
    theme(legend.position="bottom") +
    guides(colour = guide_legend("Data collection"), shape = guide_legend("Data collection"), mapping = guide_legend("Data collection")) +
    scale_colour_grey()
  print(graph_all)
  ggsave(paste0("gradual_effects_all_", model, ".jpg"), width = 20, height = 12, units = "cm")

  # graph with meta analysis only
  meta <- subset(output_temp, output_temp$sample == "FE Meta-analysis")

  graph_meta <- ggplot(meta, aes(x=cutoff, y=effect_estimate, colour = sample))+
    geom_hline(yintercept=0, linetype="dashed", size = 0.2, color="grey") +
    stat_smooth(method ='lm', se = T, size = 0.25, alpha = 0.2) +
    geom_point( size = 2, aes(shape = sample, colour = sample),position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.1, position=position_dodge(width=0.5)) +
    facet_grid(factor(effect, levels=c('Curiosity','Monetary incentive','Interaction')) ~ ., scale = "free", switch = "y") +
    theme_classic() +
    #theme(strip.text.y = element_text(size = 8)) +
    scale_x_continuous(breaks = c(0:5), labels = labels) +
    scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "Odds Ratio")) +
    #labs(x="Behavioural recognition memory performance", y=paste("b estimate (FE meta-analysis integration of", model, "gLME)"), color = "Data collection", title = paste("Fixed effects as a function of confidence cut-off")) +
    labs(x="Behavioural recognition memory performance", y=paste("b estimate (FE meta-analysis integration)"), color = "Data collection") +
    #theme(axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), title=element_text(size =14, face="bold"), legend.title = element_text(size=12), legend.text = element_text(size =12), strip.text = element_text(size =12)) +
    theme(legend.position="none") +
    #guides(colour = guide_legend("Data collection"), shape = guide_legend("Data collection"), mapping = guide_legend("Data collection")) +
    scale_colour_grey()
  print(graph_meta)
  ggsave(paste0("gradual_effects_meta_", model, ".jpg"), width = 15, height = 11, units = "cm")

  # graph with studies only
  studies <- subset(output_temp, output_temp$sample != "FE Meta-analysis")
  graph_studies <- ggplot(studies, aes(x=cutoff, y=effect_estimate, colour = sample, fill = sample))+
    geom_hline(yintercept=0, linetype="dashed", size = 0.2, color="grey") +
    stat_smooth(method ='lm', se = T, size = 0.5, alpha = 0.2, position=position_dodge(width=0.2)) +
    geom_point( size = 2, aes(shape = sample, colour = sample),position=position_dodge(width=0.2)) +
    geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.1, position=position_dodge(width=0.2)) +
    facet_grid(factor(effect, levels=c('Curiosity','Monetary incentive','Interaction')) ~ ., scale = "free", switch = "y") +
    theme_classic() +
    #theme(strip.text.y = element_text(size = 8)) +
    scale_x_continuous(breaks = c(0:5), labels = labels) +
    scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "Odds Ratio")) +
    #labs(x="Behavioural recognition memory performance", y=paste(model, "gLME FE beta estimate"), color = "Data collection", title = paste("Fixed effects as a function of confidence cut-off")) +
    labs(x="Behavioural recognition memory performance", y=paste("b estimate (gLME model -", model, "RE structure)"), color = "Data collection") +
    #theme(axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), title=element_text(size =14, face="bold"), legend.title = element_text(size=12), legend.text = element_text(size =12), strip.text = element_text(size =12)) +
    theme(legend.position="bottom") +
    guides(fill = guide_legend("Data collection"), colour = guide_legend("Data collection"), shape = guide_legend("Data collection"), mapping = guide_legend("Data collection"))
  print(graph_studies)
  ggsave(paste0("gradual_effects_studies_", model, ".jpg"), width = 15, height = 12, units = "cm")


  # COMPUTE RELATIONSHIP BETWEEN EFFECT SIZE AND CONFIDENCE THRESHOLD
  cutoff <- c(0:5) # confidence threshold
  
  # linear model curiosity
  curiosityBeta <- output_temp$effect_estimate[output_temp$sample == "FE Meta-analysis" & output_temp$effect == "Curiosity"]
  lm_curiosityBeta <- summary(lm(curiosityBeta ~ cutoff))
  # model coefficients
  coef_curiosity <- as.data.frame(lm_curiosityBeta$coefficients)
  coef_curiosity$Effect <- "Curiosity"
  coef_curiosity <- merge(coef_curiosity, as.data.frame(confint(lm(curiosityBeta ~ cutoff))), by = 0)
  coef_curiosity[,2:5] <- format(round(coef_curiosity[,2:5], digits = 3), nsmall = 3) # round values
  coef_curiosity[,7:8] <- format(round(coef_curiosity[,7:8], digits = 3), nsmall = 3) # round values
  # overall model fit
  fstat <- format(round(lm_curiosityBeta$fstatistic["value"], 3), nsmall = 3)
  df1 <- round(lm_curiosityBeta$fstatistic["numdf"])
  df2 <- round(lm_curiosityBeta$fstatistic["dendf"])
  pval <- format(round(pf(lm_curiosityBeta$fstatistic[1], lm_curiosityBeta$fstatistic[2], lm_curiosityBeta$fstatistic[3],lower.tail=F), 3), nsmall = 3)
  r2 <- format(round(lm_curiosityBeta$r.squared, 3), nsmall = 3)
  fit_curiosityBeta <- paste(paste0("F(", df1, ",", df2, ") = ", fstat), paste0("p = ", pval), paste0( "R2 = ", r2), sep = '\n')
  fit_curiosityBeta <- paste0( "R2 = ", r2)
  #coef_curiosity$Fit <- ""
  #coef_curiosity <- rbind(coef_curiosity, c("", "", "", "", "", "Curiosity", "", "",  fit_curiosityBeta))

  # linear model incentive
  incentiveBeta <- output_temp$effect_estimate[output_temp$sample == "FE Meta-analysis" & output_temp$effect == "Monetary incentive"]
  lm_incentiveBeta <- summary(lm(incentiveBeta ~ cutoff))
  # model coefficients
  coef_incentive <- as.data.frame(lm_incentiveBeta$coefficients)
  coef_incentive$Effect <- "Monetary incentive"
  coef_incentive <- merge(coef_incentive, as.data.frame(confint(lm(incentiveBeta ~ cutoff))), by = 0)
  coef_incentive[,2:5] <- format(round(coef_incentive[,2:5], digits = 3), nsmall = 3) # round values
  coef_incentive[,7:8] <- format(round(coef_incentive[,7:8], digits = 3), nsmall = 3) # round values
  # overall model fit
  fstat <- format(round(lm_incentiveBeta$fstatistic["value"], 3), nsmall = 3)
  df1 <- round(lm_incentiveBeta$fstatistic["numdf"])
  df2 <- round(lm_incentiveBeta$fstatistic["dendf"])
  pval <- format(round(pf(lm_incentiveBeta$fstatistic[1], lm_incentiveBeta$fstatistic[2], lm_incentiveBeta$fstatistic[3],lower.tail=F), 3), nsmall = 3)
  r2 <- format(round(lm_incentiveBeta$r.squared, 3), nsmall = 3)
  fit_incentiveBeta <- paste(paste0("F(", df1, ",", df2, ") = ", fstat), paste0("p = ", pval), paste0( "R2 = ", r2), sep = '\n')
  fit_incentiveBeta <- paste0( "R2 = ", r2)
  #coef_incentive$Fit <- ""
  #coef_incentive <- rbind(coef_incentive, c("", "", "", "", "", "Monetary incentive", "", "",  fit_incentiveBeta))
  
  # linear model interaction
  interactionBeta <- output_temp$effect_estimate[output_temp$sample == "FE Meta-analysis" & output_temp$effect == "Interaction"]
  lm_interactionBeta <- summary(lm(interactionBeta ~ cutoff))
  # model coefficients
  coef_interaction <- as.data.frame(lm_interactionBeta$coefficients)
  coef_interaction$Effect <- "Interaction"
  coef_interaction <- merge(coef_interaction, as.data.frame(confint(lm(interactionBeta ~ cutoff))), by = 0)
  coef_interaction[,2:5] <- format(round(coef_interaction[,2:5], digits = 3), nsmall = 3) # round values
  coef_interaction[,7:8] <- format(round(coef_interaction[,7:8], digits = 3), nsmall = 3) # round values
  # overall model fit
  fstat <- format(round(lm_interactionBeta$fstatistic["value"], 3), nsmall = 3)
  df1 <- round(lm_interactionBeta$fstatistic["numdf"])
  df2 <- round(lm_interactionBeta$fstatistic["dendf"])
  pval <- format(round(pf(lm_interactionBeta$fstatistic[1], lm_interactionBeta$fstatistic[2], lm_interactionBeta$fstatistic[3],lower.tail=F), 3), nsmall = 3)
  r2 <- format(round(lm_interactionBeta$r.squared, 3), nsmall = 3)
  fit_interactionBeta <- paste(paste0("F(", df1, ",", df2, ") = ", fstat), paste0("p = ", pval), paste0( "R2 = ", r2), sep = '\n')
  fit_interactionBeta <- paste0( "R2 = ", r2)
  #coef_interaction$Fit <- ""
  #coef_interaction <- rbind(coef_interaction, c("", "", "", "", "", "Interaction", "", "",  fit_interactionBeta))

  # combine all three effects
  lm_beta <- rbind(coef_curiosity, coef_incentive, coef_interaction)
  lm_beta$Parameter <- c("Intercept", "Slope")
  row.names(lm_beta) <- NULL
  lm_beta <- lm_beta[,c("Effect", "Parameter", "Estimate", "Std. Error", "2.5 %", "97.5 %",  "t value", "Pr(>|t|)")]

  # transform values to CI
  lm_beta$`2.5 %` <- paste0("[", lm_beta$`2.5 %`,"; ", lm_beta$`97.5 %`, "]")
  lm_beta$`97.5 %` <- NULL
  names(lm_beta)[names(lm_beta) == "2.5 %"] <- "95%-CI"
  names(lm_beta)[names(lm_beta) == "Estimate"] <- "b (SE)"
  names(lm_beta)[names(lm_beta) == "Std. Error"] <- "SE"
  names(lm_beta)[names(lm_beta) == "Pr(>|t|)"] <- "p value"

  lm_beta$`b (SE)` <- paste0(lm_beta$`b (SE)`, " (", lm_beta$SE, ")")
  lm_beta$SE <- NULL
  lm_beta$Fit <- ""
  
  # add column spanners and model fit and remove effect column
  lm_beta <- rbind(
    # curiosity 
    c(paste0(effects[1]), paste0(effects[1]), "", "", "", "", ""), 
    lm_beta[lm_beta$Effect == paste0(effects[1]),],
    c(paste0(effects[1]), "", "", "", "", "", fit_curiosityBeta),
    # incentive
    c(paste0(effects[2]), paste0(effects[2]), "", "", "", "", ""), 
    lm_beta[lm_beta$Effect == paste0(effects[2]),],
    c(paste0(effects[2]), "", "", "", "", "", fit_incentiveBeta),
    # interaction               
    c(paste0(effects[3]), paste0(effects[3]), "", "", "", "", ""), 
    lm_beta[lm_beta$Effect == paste0(effects[3]),],
    c(paste0(effects[3]), "", "", "", "", "", fit_interactionBeta)
    )
  lm_beta$Effect <- NULL
  lm_beta$id_col <- c(1:dim(lm_beta)[1])

  # remove not needed objects
  rm(incentiveBeta, lm_incentiveBeta, coef_incentive, curiosityBeta, lm_curiosityBeta, coef_curiosity, interactionBeta, lm_interactionBeta, coef_interaction)
  
  # save output
  if (model == "full"){
    output_full <- output
    meta_full <- meta_table
    beta_full <- lm_beta
    main_dv_full <- meta_main
    recog_full <- meta_recog
    lme_sample_full <- lme_sample_new
  } else { 
    output_reduced <- output
    meta_reduced <- meta_table
    beta_reduced <- lm_beta
    main_dv_reduced <- meta_main
    recog_reduced <- meta_recog
    lme_sample_reduced <- lme_sample_new
    
  }
  #rm(output, meta_out, lm_beta, meta_main_dvs, meta_main, meta_recog, lme_sample, lme_sample_new)
  
}

# lm_beta tables to create table_s6
beta_both <- merge(beta_full, beta_reduced, by = c("id_col"))
beta_both$Parameter.y <- NULL
beta_both$id_col <- NULL
beta_both <- rbind(c("", "Full gLME model", "", "", "", "", "Reduced gLME model", "", "", "", ""),
                   c("", "B (SE)", "95%-CI", "t value", "p value", "Fit", "B (SE)", "95%-CI", "t value", "p value", "Fit"),
                   beta_both)

# go back into analysis directory
setwd(analysis_dir)

# save all tables
xlsx::write.xlsx(main_dv_full, file=filename_tables, sheetName = "Table_2", append = T, row.names = F) 
#xlsx::write.xlsx(beta_full, file=filename_tables, sheetName = "Table_3", append = T, row.names = F) 
xlsx::write.xlsx(recog_full, file=filename_tables, sheetName = "Table_S3", append = T, row.names = F) 
xlsx::write.xlsx(lme_sample_full, file=filename_tables, sheetName = "Table_S4", append = T, row.names = F, col.names = F) 

xlsx::write.xlsx(meta_reduced, file=filename_tables, sheetName = "Table_S5", append = T, row.names = F) # all FE meta results for all 7 DVs
xlsx::write.xlsx(beta_both, file=filename_tables, sheetName = "Table_S6", append = T, row.names = F, col.names = F) # lm betas for both models
xlsx::write.xlsx(lme_sample_reduced, file=filename_tables, sheetName = "Table_S7", append = T, row.names = F, col.names = F)  # equivalent to table S4




# print header to markdown
cat("  \n### Comparison between full and reduced model")
cat("  \n The output of both moodels was compared by determining whether the numbers were positive (or negative) for both models.")
cat("  \n For this, the following formula was used: sum((estimate from full model > 0) != (estimate from reduced model > 0))")
cat("  \n This comparison was done for the beta estimate, as well as the lower and upper bound of the 95% confidence interval.")
cat("  \n ")
cat("  \n\n",paste("The beta estimate differed between full and reduced model in", sum((output_full$effect_estimate > 0) != (output_reduced$effect_estimate > 0)), "occasion(s)."))
cat("  \n\n",paste("The lower bound of the confidence interval differed between full and reduced model in", sum((output_full$ci_upper > 0) != (output_reduced$ci_upper > 0)), "occasion(s)."))
cat("  \n\n",paste("The lower bound of the confidence interval differed between full and reduced model in", sum((output_full$ci_lower > 0) != (output_reduced$ci_lower > 0)), "occasion(s)."))




###### rain cloud plots #####

library(ggplot2)

# RECOGNITION

data <- dfLong_all %>%
  group_by(Username, group, recognition) %>% 
  summarise(mean_curiosity = mean(responseCuriosity, na.rm = TRUE))

data$cat <- ifelse(data$group == "exp" & data$recognition == 1, "incentive R",
                   ifelse(data$group == "exp" & data$recognition == 0, "incentive F",
                          ifelse(data$group == "cont" & data$recognition == 1, "cont R",
                                 ifelse(data$group == "cont" & data$recognition == 0, "cont F",NA))))


library(ggplot2)
ggplot(data, aes(x = as.factor(cat), y = mean_curiosity)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")


# HIGH CONF RECOGNITION

data <- dfLong_all %>%
  group_by(Username, group, recognitionConfLevel_4_5_6) %>% 
  summarise(mean_curiosity = mean(responseCuriosity, na.rm = TRUE))

data$cat <- ifelse(data$group == "exp" & data$recognitionConfLevel_4_5_6 == 1, "incentive R",
                   ifelse(data$group == "exp" & data$recognitionConfLevel_4_5_6 == 0, "incentive F",
                          ifelse(data$group == "cont" & data$recognitionConfLevel_4_5_6 == 1, "cont R",
                                 ifelse(data$group == "cont" & data$recognitionConfLevel_4_5_6 == 0, "cont F",NA))))

ggplot(data, aes(x = as.factor(cat), y = mean_curiosity)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")




# CUED RECALL

data <- dfLong_all %>%
  group_by(Username, group, cuedRecallStrict) %>% 
  summarise(mean_curiosity = mean(responseCuriosity, na.rm = TRUE))

data$cat <- ifelse(data$group == "exp" & data$cuedRecallStrict == 1, "incentive R",
                   ifelse(data$group == "exp" & data$cuedRecallStrict == 0, "incentive F",
                          ifelse(data$group == "cont" & data$cuedRecallStrict == 1, "cont R",
                                 ifelse(data$group == "cont" & data$cuedRecallStrict == 0, "cont F",NA))))


#data <- subset(dfLong_all, sample = "main")

library(ggplot2)
ggplot(data, aes(x = as.factor(cat), y = mean_curiosity)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")



# HIGHEST CONFIDENCE ONLY

data <- dfLong_all %>%
  group_by(Username, group, recognitionConfLevel_6) %>% 
  summarise(mean_curiosity = mean(responseCuriosity, na.rm = TRUE))

data$cat <- ifelse(data$group == "exp" & data$recognitionConfLevel_6 == 1, "incentive R",
                   ifelse(data$group == "exp" & data$recognitionConfLevel_6 == 0, "incentive F",
                          ifelse(data$group == "cont" & data$recognitionConfLevel_6 == 1, "cont R",
                                 ifelse(data$group == "cont" & data$recognitionConfLevel_6 == 0, "cont F",NA))))

ggplot(data, aes(x = as.factor(cat), y = mean_curiosity)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

