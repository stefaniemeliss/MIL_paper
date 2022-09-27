#####################################################################################
######################## script pre-processing task paper ###########################
#####################################################################################

## this script creates the dataset files pooling data across all 3 experiments ##

# empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
analysis_dir <- getwd()

# helper functions and packages #
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/errorbars.R?raw=TRUE")
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/rbindcolumns.R?raw=TRUE")

library(lme4)
library(psych)
library(ggplot2)
library(dplyr)

# define version 
versions <- c("pilot", "main", "fmri")

# define OSF directory
osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "data") # looks at all files and directories in the project and defines the match with "data"

for (version_official in versions){
  
  # downlaod data sets from OSF
  sub_dir <- osfr::osf_mkdir(target_dir, path = paste0(version_official)) # add folder in OSF data dir
  
  osfr::osf_ls_files(sub_dir, pattern = "MagicBehavioural") %>%
    osfr::osf_download(conflicts = "overwrite")
  
  # read in datasets
  wide <- read.csv(paste0("wide_MagicBehavioural_", version_official, ".csv"), stringsAsFactors = F)
  long <- read.csv(paste0("long_MagicBehavioural_", version_official, ".csv"), stringsAsFactors = F)
  
  # remove dataset files after reading in 
  file.remove(paste0("wide_MagicBehavioural_", version_official, ".csv"))
  file.remove(paste0("long_MagicBehavioural_", version_official, ".csv"))
                          
  # add study column
  wide$sample <- version_official
  long$sample <- version_official
  
  # rename columns to match 
  if (version_official == "fmri"){
    # long format
    names(long)[names(long)=="BIDS"] <- "Username"  
    names(long)[names(long)=="trial"] <- "trialMain"  
    # wide format
    names(wide)[names(wide)=="BIDS"] <- "Username"  
    names(wide)[names(wide)=="endExperiment_UTC"] <- "endMain_UTC"  
    names(wide)[names(wide)=="startExperiment_UTC"] <- "startMain_UTC"  
    
    

  }
  
  # reorder column
  wide_2 <- wide[, c(dim(wide)[2], 1:dim(wide)[2]-1)]
  long_2 <- long[, c(dim(long)[2], 1:dim(long)[2]-1)]
  

  # combine dataframes
  if (version_official == "pilot"){
    dfWide <- wide_2
    dfLong <- long_2
  } else {
    dfWide <- rbind.all.columns(dfWide, wide_2)
    dfLong <- rbind.all.columns(dfLong, long_2)
  }
  
  rm(wide, wide_2, long, long_2)
}

# remove columns
dfLong <- subset(dfLong, select = -c(motivation, orderNumber, block,acq, triggerTaskBlockRaw,startBlock,endBlock,
                            vidFileName,whichVid,tTrialStart,tTrialEnd, durationTrial, fixationInitialDuration, pre_stim_rest, timingCorrection, 
                            displayVidOnset, displayVidOffset, displayVidDuration, fixationPostVidOnset, fixationPostVidDuration, jitterVideo_trial, displayAnswerOnset, displayAnswerDuration, timeoutAnswer, timestampAnswer,timestampAnswerWhite, 
                            fixationPostAnswerOnset,fixationPostAnswerDuration, betweenRatingFixation, displayCuriosityOnset, displayCuriosityDuration, timeoutCuriosity, timestampCuriosity,timestampCuriosityWhite, startValueCuriosity, clicksCuriosity,fixationPostCuriosityOnset,fixationPostCuriosityDuration,	
                            jitterRating_trial,mockOffset,cueImage,momentOfSurprise_1,momentOfSurprise_2,momentOfSurprise_3,momentOfSurprise_4,momentOfSurprise_5,momentOfSurprise_6,momentOfSurprise_7,additionalMarker_momentOfSurprise_1,additionalMarker_momentOfSurprise_2,fMRI,
                            startExperiment_UTC, endExperiment_UTC,endPractice_UTC,meanCuriosity_MAGMOT,meanCuriosityStandardised_MAGMOT,mediansplitCuriosity_MAGMOT,avgVidDur_MAGMOT,curiosity,curiosity_RT,curiosity_tooSlow,mediansplitCuriosityWithinSubject_updated,curiosityGroupMeanCentered_updated,rewardByCuriosity_updated,answer,answer_RT,answer_tooSlow,responseEstimate,rtEstimate,displayBlankDuration) )

dfWide <- subset(dfWide, select = -c(completeData,post1,post2,post3,post4,post5,post6,post7,post8,post9,post10,post11,post12,post13,post14,post15,post16,post17,post18,post19,post20,post21,post22,post23,post24,post25,post26,
                                     post1_score,post2_score,post3_score,post4_score,post5_score,post6_score,post7_score,post8_score,post9_score,post10_score,post11_score,post12_score,post13_score,post14_score,post15_score,post16_score,post17_score,post18_score,post19_score,post20_score,post21_score,post22_score,
                                     answer_tooSlow,curiosity_tooSlow,endPractice_UTC,durScanningSession,durInSecs_firstBlock,durInMins_firstBlock,lmBeta_cuedRecallStrict,lmBeta_cuedRecallLenient,lmBeta_allConf,lmBeta_highConf,lmBeta_aboveAvgConf,lmBeta_rememberedStrictAboveAvg,lmBeta_rememberedLenientAboveAvg,lmBeta_rememberedStrictHigh,lmBeta_rememberedLenientHigh,durInSecs_secondBlock,durInMins_secondBlock,durInSecs_thirdBlock,durInMins_thirdBlock,durInSecs,durInMins,curBetaRed_cuedRecallStrict,curBetaRed_c_cuedRecallStrict,curBetaRed_cuedRecallLenient,curBetaRed_c_cuedRecallLenient,curBetaRed_allConf,curBetaRed_c_allConf,curBetaRed_highConf,curBetaRed_c_highConf,curBetaRed_aboveAvgConf,curBetaRed_c_aboveAvgConf,curBetaRed_rememberedStrictAboveAvg,curBetaRed_c_rememberedStrictAboveAvg,curBetaRed_rememberedLenientAboveAvg,curBetaRed_c_rememberedLenientAboveAvg,curBetaRed_rememberedStrictHigh,curBetaRed_c_rememberedStrictHigh,curBetaRed_rememberedLenientHigh,curBetaRed_c_rememberedLenientHigh,RSFC_VTAHPC_run1,RSFC_VTAHPC_run2,RSFC_VTAHPC_run1_z,RSFC_VTAHPC_run2_z,RSFC_VTAHPC_diff,RSFC_VTAHPC_run1_spearman,RSFC_VTAHPC_run2_spearman,RSFC_VTAHPC_run1_z_spearman,RSFC_VTAHPC_run2_z_spearman,RSFC_VTAHPC_diff_spearman,
                                     BISBAS_inhibition,BISBAS_rewardresponsiveness,BISBAS_drive,BISBAS_funseeking,NeedForCognition,FearOfFailure,ApproachTemperament,AvoidanceTemperament,TraitCuriosity,corsiSpan,NBacks,nonNBacks,nback_hits,nback_misses_inclTooSlow,nback_misses_exclTooSlow,nback_correctrejections,nback_falsealarms_inclTooSlow,nback_falsealarms_exclTooSlow,nback_hitrate,nback_falsealarmrate,nback_accurary,postFile,startPost,endPost,durPost,StateCuriosity,
                                     alcohol,alcoholAmount,eyetracking,fieldmap,preLearningRest,taskBlock1,taskBlock2,taskBlock3,postLearningRest,T1w,eyetrackingData,taskData,questionnaireData,motivation, fMRI, memoryFile, corsiFile,nbackFile,preFile,startPre,endPre,durPre, durMemory, durMainExp, DOB,sex))

# save files
write.csv(dfLong, "long_MagicBehavioural.csv", row.names = F)
write.csv(dfWide, "wide_MagicBehavioural.csv", row.names = F)

