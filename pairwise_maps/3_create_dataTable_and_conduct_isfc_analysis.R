#####################################################################################
########################### script analysis task paper ##############################
#####################################################################################


# this script creates the dataTable files needed for AFNI's 3dISC command to specify the LME-CRE to analyse

rm(list = ls())
options(scipen=999)

library(dplyr)
library(lme4)

# define directories
anal_dir <- getwd() # current directory where the project file is saved
ISFC_dir <- file.path(anal_dir, "ISFC") # directory including the extracted ROI time series
trials_dir <- file.path(anal_dir, "trials_files") # directory including files with onset, duration, and behavioural measures of each trial

# define root directory and go back to anal_dir
setwd('..')
root_dir <- getwd()
setwd(anal_dir)

# define output directoty for tables
behav_dir <- file.path(root_dir, "behavioural")

# determine file name of output excel file
filename_tables <- "tables_task_paper.xlsx"

# delete any dataTable files
file.remove(list.files(pattern = "dataTable_IS"))

# determine memory labels
# memoryLevels <- c("cuedRecallStrict", 
#                      "recognition", "recognitionConfLevel_4_5_6")
# memoryLabels <- c("cuedRecallStrict",
#                      "allConf", "highConf")
memoryLevels <- c("recognitionConfLevel_4_5_6")
memoryLabels <- c("highConf")

# define version 
version <- "MAGMOT"
version_official <- "fmri"
dataset_name <- "MMC"


### downlaod behavioural datasets from OSF ###
MMC_project <- osfr::osf_retrieve_node("eyzwb") 

osfr::osf_ls_files(MMC_project, pattern = paste0(dataset_name, "_demographics.csv")) %>%
  osfr::osf_download(conflicts = "overwrite")
osfr::osf_ls_files(MMC_project, pattern = paste0(dataset_name, "_experimental_data.csv")) %>%
  osfr::osf_download(conflicts = "overwrite")

# read in behavioural data
MAGMOT <- read.csv(paste0(dataset_name, "_demographics.csv"), stringsAsFactors = F)
dfLong <- read.csv(paste0(dataset_name, "_experimental_data.csv"), stringsAsFactors = F)

# # delete behavioural data
# file.remove(paste0(dataset_name, "_demographics.csv"))
# file.remove(paste0(dataset_name, "_experimental_data.csv"))

### compute memory performance scores ###

# effect-code group
dfLong$group_ec <- ifelse(dfLong$group == "exp", 1, -1) # effect code group

# center curiosity:  within cluster (CWC, i.e., group-mean centering)
subjects <- MAGMOT$BIDS
dfLong$curiosity_cwc <- NA
for (s in seq_along(subjects)) {
  dfLong[dfLong$BIDS == subjects[s], "curiosity_cwc"] <- dfLong[dfLong$BIDS == subjects[s], "responseCuriosity"] - mean(dfLong[dfLong$BIDS == subjects[s], "responseCuriosity"], na.rm = T)
}

#  make sure that BIDS and stimID are factors    
dfLong$BIDS <- as.factor(as.character(dfLong$BIDS)) 
dfLong$stimID <- as.factor(as.character(dfLong$stimID))

for (mem in 1:length(memoryLevels)) {
  
  # compute full gLME model: random intercept for stimulus ID and participant, random slope for curiosity
  print(paste("full model",memoryLevels[mem]))
  
  # full Random Effects structure
  gLME_full <- glmer(dfLong[, memoryLevels[mem]] ~ group_ec*curiosity_cwc + (1+curiosity_cwc|BIDS) + (1|stimID) , family = "binomial"(link = 'logit'), data = dfLong)
  temp <- coef(gLME_full)$BIDS # put coefficients into temp object
  temp$BIDS <- row.names(temp) # add column with BIDS
  row.names(temp) <- NULL # remove row names
  temp <- temp[, c("BIDS", "curiosity_cwc")] # reduce to BIDS column and curiosity betas
  names(temp)[names(temp) == "curiosity_cwc"] <- paste0("curBetaFull_", memoryLabels[mem]) # rename variable
  MAGMOT <- merge(MAGMOT, temp, by = "BIDS") # merge with dfWide
  rm(temp)
  MAGMOT[[paste0("curBetaFull_c_", memoryLabels[mem])]] <-  scale(MAGMOT[[paste0("curBetaFull_", memoryLabels[mem])]], center = T, scale = F) # center values
  
  # compute reduced gLME model: random intercept for participant, random slope for curiosity
  print(paste("reduced model",memoryLevels[mem]))
  
  # reduced randmom effects structure
  gLME_red <- glmer(dfLong[, memoryLevels[mem]] ~ group_ec*curiosity_cwc + (1+curiosity_cwc|BIDS), family = "binomial"(link = 'logit'), data = dfLong)
  temp <- coef(gLME_red)$BIDS # put coefficients into temp object
  temp$BIDS <- row.names(temp) # add column with BIDS
  row.names(temp) <- NULL # remove row names
  temp <- temp[, c("BIDS", "curiosity_cwc")] # reduce to BIDS column and curiosity betas
  names(temp)[names(temp) == "curiosity_cwc"] <- paste0("curBetaRed_", memoryLabels[mem]) # rename variable
  MAGMOT <- merge(MAGMOT, temp, by = "BIDS") # merge with dfWide
  rm(temp)
  MAGMOT[[paste0("curBetaRed_c_", memoryLabels[mem])]] <-  scale(MAGMOT[[paste0("curBetaRed_", memoryLabels[mem])]], center = T, scale = F) # center values
  
}


# get a list with all subject BIDS strings
subjectsCorr <- subjects
subjectsToCorrelate <- subjectsCorr[-1]

N <- 0.5*length(subjectsCorr)*length(subjectsToCorrelate)

dataTable_ISFC <- data.frame()
x <- 0

vol <- 1:594 # volumes in concatenated time series
ROI_data <- as.data.frame(vol)


if ( identical(subjectsToCorrelate, character(0)) == F){
  
  # loop over subject s
  for(s in seq_along(subjectsCorr)){
    
    # determine subject s1
    s1 <- subjectsCorr[s]
    
    # behavioural data #
    
    # determine behavioural timecourse files
    file_s1 <- list.files(path = trials_dir, pattern = paste0(s1, "_task-magictrickwatching_trials.tsv"))
    events_s1 <- read.delim(file = file.path(trials_dir, file_s1), header = T, sep="\t", na = "n/a")
    rm(file_s1)
    
    # pick relevant columns
    events_s1 <- events_s1[,c("stim_file", "responseCuriosity", 
                            "trial_type_cuedRecallStrict", 
                            "trial_type_allConf", "trial_type_highConf",  
                            "responseConfidence")]
    
    # mean-center curiosity and confidence
    events_s1$responseCuriosity <- events_s1$responseCuriosity - mean(events_s1$responseCuriosity, na.rm = T) #mean center curiosity
    events_s1$responseConfidence <- events_s1$responseConfidence - mean(events_s1$responseConfidence, na.rm = T) #mean center confidence
    
    # change column names so that they include the subject ID
    names(events_s1) <- paste0(names(events_s1), "_", s1)
    names(events_s1)[names(events_s1)== paste0("stim_file_", s1)] <- "stim_file"
    
    # ROI data #
    
    # determine ROI timeseries
    HPC_s1 <- read.delim(file = file.path(ISFC_dir, list.files(path = ISFC_dir, pattern = paste0(s1, "_task-magictrickwatching_maskave_aHPC.txt"))), header = F) # average aHPC
    HPC_s1 <- HPC_s1$V1
    
    VTA_s1 <- read.delim(file = file.path(ISFC_dir, list.files(path = ISFC_dir, pattern = paste0(s1, "_task-magictrickwatching_maskave_VTA.txt"))), header = F) # average VTA
    VTA_s1 <- VTA_s1$V1
    
    # add to data frame
    if (s == 1) {
      
      ROI_data[[paste0("aHPC_", s1)]] <- HPC_s1
      ROI_data[[paste0( "VTA_", s1)]] <- VTA_s1
      
      }
    
    # update subject list
    subjectsToCorrelate <-subjectsCorr[-c(1:s)]

    # loop over subject s+1
    for(ss in seq_along(subjectsToCorrelate)){
      
      # determine subject s2
      s2 <- subjectsToCorrelate[ss]
      
      # behavioural data #
      
      # determine behavioural timecourse files
      file_s2 <- list.files(path = trials_dir, pattern = paste0(s2, "_task-magictrickwatching_trials.tsv"))
      events_s2 <- read.delim(file = file.path(trials_dir, file_s2), header = T, sep="\t", na = "n/a")
      rm(file_s2)
      
      # pick relevant columns
      events_s2 <- events_s2[,c("stim_file", "responseCuriosity", 
                                "trial_type_cuedRecallStrict", 
                                "trial_type_allConf", "trial_type_highConf",  
                                "responseConfidence")]
      
      # mean-center curiosity and confidence
      events_s2$responseCuriosity <- events_s2$responseCuriosity - mean(events_s2$responseCuriosity, na.rm = T) #mean center curiosity
      events_s2$responseConfidence <- events_s2$responseConfidence - mean(events_s2$responseConfidence, na.rm = T) #mean center confidence
      
      # change column names so that they include the subject ID
      names(events_s2) <- paste0(names(events_s2), "_", s2)
      names(events_s2)[names(events_s2)== paste0("stim_file_", s2)] <- "stim_file"
      
      # merge both subjects
      # using stim_file (e.g., 01_H35_combined_small.mp4) de-shuffles the trials and arranges the data in ascending order
      events <- merge(events_s1, events_s2, by = c("stim_file")) 
      
      for (mem in 1:length(memoryLevels)) {
        # behavioural parcellation: define match in memory performance
        events[[paste0("behavParcel_", memoryLabels[mem])]] <- ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",s1)]] == "remembered" & events[[paste0("trial_type_", memoryLabels[mem], "_", s2)]] == "remembered", "bothRemembered",
                                                                      ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",s1)]] == "forgotten" & events[[paste0("trial_type_", memoryLabels[mem], "_", s2)]] == "forgotten", "bothForgotten",
                                                                             "differentResponses"))
        # create an effect coded memory variable
        events[[paste0("trial_type_", memoryLabels[mem],"_",s1, "_effectCoded")]] <- ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",s1)]] == "remembered", 1, -1)
        events[[paste0("trial_type_", memoryLabels[mem],"_",s2, "_effectCoded")]] <- ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",s2)]] == "remembered", 1, -1)
      }
      
      # create effect coded group variable
      events[[paste0("reward_",s1, "_effectCoded")]] <- ifelse(grepl("cont", s1), -1, 1)
      events[[paste0("reward_",s2, "_effectCoded")]] <- ifelse(grepl("cont", s2), -1, 1)
      
      # ROI data #
      
      # determine ROI timeseries
      HPC_s2 <- read.delim(file = file.path(ISFC_dir, list.files(path = ISFC_dir, pattern = paste0(s2, "_task-magictrickwatching_maskave_aHPC.txt"))), header = F) # average aHPC
      HPC_s2 <- HPC_s2$V1
      
      VTA_s2 <- read.delim(file = file.path(ISFC_dir, list.files(path = ISFC_dir, pattern = paste0(s2, "_task-magictrickwatching_maskave_VTA.txt"))), header = F) # average VTA
      VTA_s2 <- VTA_s2$V1
      
      # add data to data frame and change column names
      if (s == 1 ) {
        
        ROI_data[[paste0("aHPC_", s2)]] <- HPC_s2
        ROI_data[[paste0( "VTA_", s2)]] <- VTA_s2

      }
      
      # compute aHPC-midbrain ISFC for each combination
      ISFC1 <- cor(HPC_s1, VTA_s2)
      ISFC2 <- cor(HPC_s2, VTA_s1)

      # compute mean ISFC for the pair after fisher's z transformation
      ISFC <- mean(atanh(ISFC1), atanh(ISFC2))

      ### fill in information to dataTable (group coding) ###
      x <- x+1
      
      ### fill in information to dataTable_ISFC (group coding) ###
      dataTable_ISFC[x,1] <- s1
      dataTable_ISFC[x,2] <- s2
      # we adopt deviation coding for the two groups by replacing two groups G1 (cont) and G2 (exp) with 0.5 and -0.5. 
      # Then add up the two values for each row (each subject pair), resulting in three possible values of 1, -1 and 0.
      dataTable_ISFC[x,3] <- ifelse(grepl("cont", s1) & grepl("cont", s2), 1, 
                                    ifelse(grepl("cont", s1) & grepl("exp", s2), 0, 
                                           ifelse(grepl("exp", s1) & grepl("exp", s2), -1,NA )))
      
      #curiosity and curiosity interaction
      dataTable_ISFC[x,4] <- cor(events[, paste0("responseCuriosity_", s1)], events[, paste0("responseCuriosity_", s2)] )
      dataTable_ISFC[x,5] <- dataTable_ISFC[x,4] * dataTable_ISFC[x,3] 
      
      # confidence and confidence interaction
      dataTable_ISFC[x,6] <- cor(events[, paste0("responseConfidence_", s1)], events[, paste0("responseConfidence_", s2)] )
      dataTable_ISFC[x,7] <- dataTable_ISFC[x,6] * dataTable_ISFC[x,3] 
      
      addCol <- 0
      for (mem in 1:length(memoryLevels)) {
        # memory and reward-memory interaction
        memoCol <- 8 + (6*addCol) # add more columns
        dataTable_ISFC[x,memoCol] <- cor(events[, paste0("trial_type_", memoryLabels[mem],"_",s1, "_effectCoded")], events[, paste0("trial_type_", memoryLabels[mem],"_",s2, "_effectCoded")])
        dataTable_ISFC[x,memoCol+1] <-  dataTable_ISFC[x,memoCol] * dataTable_ISFC[x,3] # third col has group information
        
        # curiosity beta and curiosity beta interaction
        betaCol <- memoCol + 2
        dataTable_ISFC[x,betaCol] <- mean(c(MAGMOT[[paste0("curBetaFull_c_", memoryLabels[mem])]][MAGMOT$BIDS == s1], MAGMOT[[paste0("curBetaFull_c_", memoryLabels[mem])]][MAGMOT$BIDS == s2])) # AnnaK formulation
        dataTable_ISFC[x,betaCol+1] <- dataTable_ISFC[x,betaCol] * dataTable_ISFC[x,3] # interaction
        
        # curiosity only beta and curiosity only beta interaction
        betaOnlyCol <- betaCol + 2
        dataTable_ISFC[x,betaOnlyCol] <- mean(c(MAGMOT[[paste0("curBetaRed_c_", memoryLabels[mem])]][MAGMOT$BIDS == s1], MAGMOT[[paste0("curBetaRed_c_", memoryLabels[mem])]][MAGMOT$BIDS == s2]))  # AnnaK formulation
        dataTable_ISFC[x,betaOnlyCol+1] <- dataTable_ISFC[x,betaOnlyCol] * dataTable_ISFC[x,3] # interaction
        
        addCol <- addCol + 1 # updates multiplier
      }
      
      # add file name and back slashs
      nextCol <- betaOnlyCol + 2 # looks at current number of columns in object and adds 1
      
      dataTable_ISFC[x,nextCol] <- ISFC
      # dataTable_ISFC[x,nextCol+1] <- ISFCfc
      
      
      dataTable_ISFC[x,nextCol+1] <- paste0("ISFC_",s1,s2,"_task-magictrickwatching_mask-ROI_z.nii.gz")
      if(x < N) {
        dataTable_ISFC[x,nextCol+2] <- '\\' #add back slash at end of the row
      }
      
      # remove data not needed anymore
      rm(HPC_s2, VTA_s2, ISFC1, ISFC2, ISFC)
      
    }
    
    # remove data not needed anymore
    rm(HPC_s1, VTA_s1)
    
    
  }
  
  # create vector with column names for dataTable
  for (mem in 1:length(memoryLevels)) {
    if (mem == 1){
      ISC_table_names <- c(paste0("corr_", memoryLabels[mem]))
    } else {
      ISC_table_names <-  c(ISC_table_names, paste0("corr_", memoryLabels[mem], collapse = ", "))
    }
    ISC_table_names <-  c(ISC_table_names, paste0("grCorr_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("curBetaFull_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("grCurBetaFull_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("curBetaRed_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("grCurBetaRed_", memoryLabels[mem], collapse = ", "))
  }
  
  # add column names to ISC table dummy
  names(dataTable_ISFC) <- c("Subj1", "Subj2", "grp", "corr_curiosity", "grCorr_curiosity", "corr_confidence", "grCorr_confidence",
                             ISC_table_names, "HPC_VTA_ISFC", "InputFile", "\\")
  
  
  ### grand-mean center values ###
  
  # create copy of dataTable
  dataTable_ISFC_c <- dataTable_ISFC
  
  # curiosity and curiosity interaction: first Fisher'z transform correlation coefficients
  dataTable_ISFC_c[,4] <- scale(atanh(dataTable_ISFC_c[,4]), center = T, scale = F) # grand mean center values
  dataTable_ISFC_c[,5] <- dataTable_ISFC_c[,4] * dataTable_ISFC_c[,3] # interaction (computed with centered vals)
  # confidence and confidence interaction: first Fisher'z transform correlation coefficients
  dataTable_ISFC_c[,6] <- scale(atanh(dataTable_ISFC_c[,6]), center = T, scale = F) # grand mean center values
  dataTable_ISFC_c[,7] <- dataTable_ISFC_c[,6] * dataTable_ISFC_c[,3] # interaction (computed with centered vals)
  
  addCol <- 0
  for (mem in 1:length(memoryLevels)) {
    # memory and memory interaction: first Fisher'z transform correlation coefficients
    memoCol <- 8 + (6*addCol) # add more columns
    dataTable_ISFC_c[,memoCol] <- scale(atanh(dataTable_ISFC_c[,memoCol]), center = T, scale = F) # grand mean center values
    dataTable_ISFC_c[,memoCol+1] <-  dataTable_ISFC_c[,memoCol] * dataTable_ISFC_c[,3] # third col has group information
    
    # curiosity beta and curiosity beta interaction
    betaCol <- memoCol + 2
    dataTable_ISFC_c[,betaCol] <- scale(dataTable_ISFC_c[,betaCol], center = T, scale = F) # grand mean center values
    dataTable_ISFC_c[,betaCol+1] <- dataTable_ISFC_c[,betaCol] * dataTable_ISFC_c[,3] # interaction (computed with centered vals)
    
    # curiosity only beta and curiosity only beta interaction
    betaOnlyCol <- betaCol + 2
    dataTable_ISFC_c[,betaOnlyCol] <- scale(dataTable_ISFC_c[,betaOnlyCol], center = T, scale = F) # grand mean center values
    dataTable_ISFC_c[,betaOnlyCol+1] <- dataTable_ISFC_c[,betaOnlyCol] * dataTable_ISFC_c[,3] # interaction (computed with centered vals)
    
    addCol <- addCol + 1 # updates multiplier
  }
  
  
  ### create data table with unique effects of curiosity, memory and their interaction ###
  dataTable_ISFC_unique <- dataTable_ISFC[, c("Subj1", "Subj2", "grp")] # copy some columns from dataTable_ISFC
  
  addCol2 <- 0
  for (mem in 1:length(memoryLevels)) {
    # unique curiosity
    curCol <- 4 + (8*addCol2) # add more columns
    uniqueCuriosity <- lm(atanh(dataTable_ISFC$corr_curiosity) ~ atanh(dataTable_ISFC[[paste0("corr_", memoryLabels[mem])]]))
    dataTable_ISFC_unique[,curCol] <- residuals(uniqueCuriosity) # residuals from model
    dataTable_ISFC_unique[,curCol] <- scale(dataTable_ISFC_unique[,curCol], center = T, scale = F) # grand mean center values
    dataTable_ISFC_unique[,curCol+1] <- dataTable_ISFC_unique[,curCol] * dataTable_ISFC_unique[,3] # interaction (computed with centered vals)
    
    # unique memory
    memCol <- curCol + 2 # add more columns
    uniqueMemory <- lm(atanh(dataTable_ISFC[[paste0("corr_", memoryLabels[mem])]]) ~ atanh(dataTable_ISFC$corr_curiosity))
    dataTable_ISFC_unique[,memCol] <- residuals(uniqueMemory) # residuals from model
    dataTable_ISFC_unique[,memCol] <- scale(dataTable_ISFC_unique[,memCol], center = T, scale = F) # grand mean center values
    dataTable_ISFC_unique[,memCol+1] <- dataTable_ISFC_unique[,memCol] * dataTable_ISFC_unique[,3] # interaction (computed with centered vals)
    
    # curiosity beta and curiosity beta interaction
    betaCol <- memCol + 2
    dataTable_ISFC_unique[,betaCol] <- dataTable_ISFC_c[,paste0("curBetaFull_", memoryLabels[mem])] # use centered values
    dataTable_ISFC_unique[,betaCol+1] <- dataTable_ISFC_c[,paste0("grCurBetaFull_", memoryLabels[mem])] # use centered values
    
    # curiosity only beta and curiosity only beta interaction
    betaOnlyCol <- betaCol + 2
    dataTable_ISFC_unique[,betaOnlyCol] <- dataTable_ISFC_c[,paste0("curBetaRed_", memoryLabels[mem])] # use centered values
    dataTable_ISFC_unique[,betaOnlyCol+1] <- dataTable_ISFC_c[,paste0("grCurBetaRed_", memoryLabels[mem])] # use centered values
    
    addCol2 <- addCol2 + 1 # updates multiplier
  }
  
  # add file name and back slashs
  nextCol2 <- betaOnlyCol + 2 # looks at current number of columns in object and adds 1
  
  #dataTable_ISFC_unique[,nextCol2] <- dataTable_ISFC[, "HPC_VTA_ISFC"]
  # dataTable_ISFC_unique[,nextCol2+1] <- dataTable_ISFC[, "HPC_fcVTA_ISFC"]
  
  dataTable_ISFC_unique[,nextCol2] <- dataTable_ISFC[, "InputFile"]
  dataTable_ISFC_unique[,nextCol2+1] <- '\\' #add back slash at end of the row
  dataTable_ISFC_unique[nrow(dataTable_ISFC_unique),ncol(dataTable_ISFC_unique)] <- NA # no // for last column
  
  # create vector with column names for dataTable
  for (mem in 1:length(memoryLevels)) {
    if (mem == 1){
      ISC_table_names <- c(paste0("uniqueCur_", memoryLabels[mem]))
    } else {
      ISC_table_names <-  c(ISC_table_names, paste0("uniqueCur_", memoryLabels[mem], collapse = ", "))
    }
    ISC_table_names <-  c(ISC_table_names, paste0("grUniqueCur_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("uniqueMem_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("grUniqueMem_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("curBetaFull_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("grCurBetaFull_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("curBetaRed_", memoryLabels[mem], collapse = ", "))
    ISC_table_names <-  c(ISC_table_names, paste0("grCurBetaRed_", memoryLabels[mem], collapse = ", "))
  }
  
  # add column names to dataTable_ISFC_unique
  names(dataTable_ISFC_unique) <- c("Subj1", "Subj2", "grp", #"uniqueConf", "grpUniqueConf",
                                       ISC_table_names,
                                    "InputFile", "\\")
}

# create dataTable_ISFC for each seed roi #
dataTable_ISFC_VTA <- dataTable_ISFC_unique
dataTable_ISFC_VTA$InputFile <- gsub("ROI", "VTA", dataTable_ISFC_VTA$InputFile)

dataTable_ISFC_HPC <- dataTable_ISFC_unique
dataTable_ISFC_HPC$InputFile <- gsub("ROI", "aHPC", dataTable_ISFC_HPC$InputFile)

dataTable_ISC <- dataTable_ISFC_unique
dataTable_ISC$InputFile <- gsub("_mask-ROI", "", dataTable_ISC$InputFile)
dataTable_ISC$InputFile <- gsub("ISFC", "ISC", dataTable_ISC$InputFile)

# save files #
write.table(dataTable_ISFC_VTA, file="dataTable_ISFC_magictrickwatching_VTA.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
write.table(dataTable_ISFC_HPC, file="dataTable_ISFC_magictrickwatching_aHPC.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
write.table(dataTable_ISC, file="dataTable_ISC_magictrickwatching.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")


### correlation of CMLE between both models
cor(MAGMOT$curBetaFull_c_highConf, MAGMOT$curBetaRed_c_highConf)
cor(dataTable_ISC$curBetaFull_highConf, dataTable_ISC$curBetaRed_highConf)

#### compute LME-CRE for aHPC-midbrain-ISFC #### 

# add HPC_VTA_ISFC values to data frames
dataTable_ISFC_unique$HPC_VTA_ISFC <- dataTable_ISFC$HPC_VTA_ISFC

# effect of incentive
lmer_incentive <- lme4::lmer(HPC_VTA_ISFC ~ grp+(1|Subj1)+(1|Subj2), data = dataTable_ISFC)
summary(lmer_incentive)

sum_LME <- summary(lmer_incentive)
LME_coef <- as.data.frame(round(sum_LME$coefficients, 3))
LME_coef$param <- c("Intercept", "Incentive effect")
LME_coef <- rbind(c("", "", "", "Incentive Effect"), LME_coef)
LME_incentive <- LME_coef[, c("param", "Estimate", "Std. Error", "t value")]
names(LME_incentive) <- c("", "Estimate", "SE", "t value")

# effect of curiosity
lmer_cur <- lme4::lmer(HPC_VTA_ISFC ~ grp+uniqueCur_highConf+grUniqueCur_highConf+(1|Subj1)+(1|Subj2), data = dataTable_ISFC_unique)
summary(lmer_cur)

# effect of memory
lmer_mem <- lme4::lmer(HPC_VTA_ISFC ~ grp+uniqueMem_highConf+grUniqueMem_highConf+(1|Subj1)+(1|Subj2), data = dataTable_ISFC_unique)
summary(lmer_mem)

# effect of CMLE (full)
lmer_clme_full <- lme4::lmer(HPC_VTA_ISFC ~ grp+curBetaFull_highConf+grCurBetaFull_highConf+(1|Subj1)+(1|Subj2), data = dataTable_ISFC_unique)
summary(lmer_clme_full)

#### put everything in one model ####

# main effects
lmer_main <- lme4::lmer(HPC_VTA_ISFC ~ grp+corr_curiosity+corr_highConf+curBetaFull_highConf+(1|Subj1)+(1|Subj2), data = dataTable_ISFC_c)
summary(lmer_main)

sum_LME <- summary(lmer_main)
LME_coef <- as.data.frame(round(sum_LME$coefficients, 3))
LME_coef$param <- c("Intercept", "Incentive effect", "Curiosity effect", "Memory effect", "CMLE effect")
LME_coef <- rbind(c("", "", "", "Main Effects"), LME_coef)
LME_main <- LME_coef[, c("param", "Estimate", "Std. Error", "t value")]
names(LME_main) <- c("", "Estimate", "SE", "t value")

# interaction effects
lmer_interaction <- lme4::lmer(HPC_VTA_ISFC ~ grp+
                                     corr_curiosity+grCorr_curiosity+
                                     corr_highConf+grCorr_highConf+
                                     curBetaFull_highConf+grCurBetaFull_highConf+
                                     (1|Subj1)+(1|Subj2), data = dataTable_ISFC_c)
summary(lmer_interaction)

sum_LME <- summary(lmer_interaction)
LME_coef <- as.data.frame(round(sum_LME$coefficients, 3))
LME_coef$param <- c("Intercept", "Incentive effect", "Curiosity effect", "Curiosity incentive interaction", "Memory effect", "Memory incentive interaction", "CMLE effect", "CMLE incentive interaction")
LME_coef <- rbind(c("", "", "", "Main and Interaction Effects"), LME_coef)
LME_int <- LME_coef[, c("param", "Estimate", "Std. Error", "t value")]
names(LME_int) <- c("", "Estimate", "SE", "t value")

# combine all aHPC-midbrain-isfc tables
LME_isfc <- rbind(LME_incentive, LME_main, LME_int)
setwd(anal_dir)
xlsx::write.xlsx(LME_isfc, file=file.path(behav_dir, filename_tables), sheetName = "Table_S11", append = T, row.names = F) 
