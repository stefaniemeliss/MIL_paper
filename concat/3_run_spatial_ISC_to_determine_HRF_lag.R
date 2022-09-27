################################################################################
# compute spatial ISC / Intersubject pattern correlation (ISPC)
################################################################################

.rs.restartR()

rm(list = ls())
gc()

library(ggplot2)
library(dplyr)

# define necessary directories
anal_dir <- getwd() # current directory where the project file is saved
fig_dir <- file.path(anal_dir, "figures")
if (dir.exists(fig_dir) == F ) {
  dir.create(fig_dir)
}

concat_dir <- file.path(anal_dir, "init_concat_V2")

file_list <- list.files(concat_dir, pattern = "V2_subject")

# prepare objects to save data in
MZ_list <- c()

# create s iterator
s <- 0

# compute intersubject pattern correlation for each subject with the mean of all other subjects #

for (file in file_list) {

  # define subject
  s <- s + 1
  subject <- gsub("V2_subject_", "", file)
  subject <- gsub(".txt", "", subject)
  print(subject)

  # measure time
  start_time <- Sys.time()

  # read in data for each subject
  df_sub <- read.table(file.path(concat_dir, paste0("V2_subject_", subject, ".txt" )))
  df_mean <- read.table(file.path(concat_dir, paste0("V2_sample_mean_", subject, ".txt" )))

  # create output data frame
  print("computing correlation")
  df_out <- data.frame()

  for (t in 1:dim(df_sub)[2]){

    for (v in 1:dim(df_mean)[2]){

      # compute correlation between
      df_out[t,v] <- cor(df_sub[, t], df_mean[, v])

    }

  }

  # save data #
  print("saving values")

  # transform df to matrix and ***fisher's z transform*** the correlation values
  M <- as.matrix(df_out)
  Z <- atanh(M)

  # add row names
  dimnames(M)[[1]] <- paste0("t", dimnames(M)[[1]])
  dimnames(Z)[[1]] <- paste0("t", dimnames(Z)[[1]])

  # save data frame with distinct name
  assign(paste0(subject, "_cor_raw"), M)
  assign(paste0(subject, "_cor_z"), Z)

  # combine data from all subjects in one list
  MZ_list[[(length(MZ_list) + 1)]] <- Z

  # remove what is not needed anymore
  rm(df_out, M, Z, df_sub, df_mean)

  # measure time
  end_time <- Sys.time()
  print(end_time - start_time)

}

# compute sample mean and transform values back to correlations #
mean_MZ_all <- apply(simplify2array(MZ_list), 1:2, mean, na.rm = T)
mean_cor_all <- tanh(mean_MZ_all)
rm(MZ_list)

# save correlation matrices #
write.csv(mean_MZ_all, file.path(concat_dir, "mean_cor_z.csv"))
write.csv(mean_cor_all, file.path(concat_dir, "mean_cor_raw.csv"))

# remove that is not needed anymore data #
rm(list=ls(pattern = "sub-"))
gc()

# save correlation matrix graph
png(file.path(fig_dir,"mean_corr_z.png"), 20000, 20000)
ggcorrplot::ggcorrplot(mean_MZ_all)
dev.off()

png(file.path(fig_dir,"mean_corr_raw.png"), 20000, 20000)
ggcorrplot::ggcorrplot(mean_cor_all)
dev.off()

# extract diagonal of matrix (Fisher's z transformed correlation coefficients)
diag <- diag(mean_MZ_all)

# extract upper off diagonal (Fisher's z transformed correlation coefficients)
upper1_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) - 1)]
upper2_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) - 2)]

# extract lower off diagonal (Fisher's z transformed correlation coefficients)
lower1_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) + 1)]
lower2_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) + 2)]


##### create a data frame with resulta of spatial ISC #####
sISC <- data.frame(vol = 1:length(diag))

sISC$diag <- diag # reliability of spatial response in V2 across subjects
sISC$upper1_diag <- c(NA, upper1_diag) # t(i) and V(i + 1)
sISC$upper2_diag <- c(NA, NA, upper2_diag) # t(i) and V(i + 2)
sISC$lower1_diag <- c(NA, lower1_diag) # V(t) and t(i + 1)
sISC$lower2_diag <- c(NA, NA, lower2_diag) # V(t) and t(i + 2)

sISC$upper1_diag <- c(upper1_diag, NA) # t(i) and V(i + 1)
sISC$upper2_diag <- c(upper2_diag, NA, NA) # t(i) and V(i + 2)
sISC$lower1_diag <- c(lower1_diag, NA) # V(t) and t(i + 1)
sISC$lower2_diag <- c(lower2_diag, NA, NA) # V(t) and t(i + 2)

# add labels to DF
dataset_name <- "MMC"

### downlaod behavioural datasets from OSF ###
MMC_project <- osfr::osf_retrieve_node("eyzwb") 

osfr::osf_ls_files(MMC_project, pattern = paste0(dataset_name, "_experimental_data.csv")) %>%
  osfr::osf_download(conflicts = "overwrite")

dataLong <- read.csv("MMC_experimental_data.csv")
TR = 2

# transform long data into wide data
displayVidDuration <- reshape::cast(dataLong, ID~vidFileName,value="displayVidDuration")

# compute mean for the display duration of each magic trick
mean_displayVidDuration <- as.data.frame(mapply(mean, displayVidDuration[-1]))
names(mean_displayVidDuration) <- "mean_displayVidDuration"
mean_displayVidDuration$stim_file <- row.names(mean_displayVidDuration)

# calculate mean_displayVidDuration without mock
mean_displayVidDuration$mean_displayVidDuration <- mean_displayVidDuration$mean_displayVidDuration - 6 # subtract mock
mean_displayVidDuration$mean_displayVidDuration_TR <- ceiling(round(mean_displayVidDuration$mean_displayVidDuration, digits = 0)/TR)

# duration of magic tricks in volumes (already ordered)
dur_vol <- mean_displayVidDuration$mean_displayVidDuration_TR

# remove data
rm(dataLong, displayVidDuration, mean_displayVidDuration)

# create vector with labels
for (t in 1:length(dur_vol)) {
  if (t == 1){
    labels <- c(rep("pre_vid_fix", 1), rep("mock", 3), rep("vid", dur_vol[t]), rep("post_vid_fix", 6))
    labels2 <- c(rep("pre_vid_fix", 1), "mock_v1","mock_v2","mock_v3", "vid_v1", rep("vid", dur_vol[t]-2), "vid_vn", "post_vid_fix", rep("fix_rate", 5))
    labels_detail <- c(rep("pre_vid_fix", 1), rep(paste0("mock_", t), 3), rep(paste0("vid_", t), dur_vol[t]), rep("post_vid_fix", 6))
  } else {
    labels <-  c(labels, rep("pre_vid_fix", 1), rep("mock", 3), rep("vid", dur_vol[t]), rep("post_vid_fix", 6))
    labels2 <- c(labels2, rep("pre_vid_fix", 1), "mock_v1","mock_v2","mock_v3", "vid_v1", rep("vid", dur_vol[t]-2), "vid_vn", "post_vid_fix", rep("fix_rate", 5))
    labels_detail <- c(labels_detail, rep("pre_vid_fix", 1), rep(paste0("mock_", t), 3), rep(paste0("vid_", t), dur_vol[t]), rep("post_vid_fix", 6))
  }

}

# add labels to df
sISC$label_detail <- labels_detail
sISC$label_lag0 <- labels
sISC$label2_lag0 <- labels2

# add labels assuming different lags
sISC$label_lag1 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-1]))
sISC$label_lag1[1] <- NA
sISC$label2_lag1 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-1]))
sISC$label2_lag1[1] <- NA
sISC$label_lag2 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-2]))
sISC$label_lag2[1:2] <- NA
sISC$label2_lag2 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-2]))
sISC$label2_lag2[1:2] <- NA
sISC$label_lag3 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-3]))
sISC$label_lag3[1:3] <- NA
sISC$label2_lag3 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-3]))
sISC$label2_lag3[1:3] <- NA
sISC$label_lag4 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-4]))
sISC$label_lag4[1:4] <- NA
sISC$label2_lag4 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-4]))
sISC$label2_lag4[1:4] <- NA
sISC$label_lag5 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-5]))
sISC$label_lag5[1:5] <- NA
sISC$label2_lag5 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-5]))
sISC$label2_lag5[1:5] <- NA
sISC$label_lag6 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-6]))
sISC$label_lag6[1:6] <- NA
sISC$label2_lag6 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-6]))
sISC$label2_lag6[1:6] <- NA


# compute average correlation for each segment and lag
index <- c("diag", "upper1_diag", "upper2_diag", "lower1_diag", "lower2_diag")

lag <- c("label2_lag0", "label2_lag1", "label2_lag2", "label2_lag3", "label2_lag4", "label2_lag5", "label2_lag6")

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


for (l in 1:length(lag)) {
  for (i in 1:length(index)) {

    # compute mean
    temp <- tapply(sISC[, index[i]], sISC[, lag[l]], mean, na.rm=TRUE)
    temp <- as.data.frame(t(temp))

    temp$lag <- lag[l]
    temp$index <- index[i]

    # combine
    if(i == 1 && l == 1){
      mean <- temp
    } else {
      mean <- rbind(mean, temp)
    }

    rm(temp)

    # compute standard error
    temp <- tapply(sISC[, index[i]], sISC[, lag[l]], stderr)
    temp <- as.data.frame(t(temp))

    temp$lag <- lag[l]
    temp$index <- index[i]

    # combine
    if(i == 1 && l == 1){
      se <- temp
    } else {
      se <- rbind(se, temp)
    }

    rm(temp)

  }
}

# change to long format and transform z values to correlation values
mean_long <- reshape2::melt(mean, id.vars = c("lag", "index"))
names(mean_long)[names(mean_long) == "value"] <- "mean_corr_z"
mean_long$mean_corr <- tanh(mean_long$mean_corr_z)
se_long <- reshape2::melt(se, id.vars = c("lag", "index"))
names(se_long)[names(se_long) == "value"] <- "se_corr_z"
se_long$se_corr <- tanh(se_long$se_corr_z)

# merge mean and SE
out_long <- merge(mean_long, se_long, by = c("lag", "index", "variable" ))

out_long$label <- ifelse(out_long$variable == "mock_v1" | out_long$variable == "mock_v2" | out_long$variable == "mock_v3", "mock",
                         ifelse(out_long$variable == "vid_v1" | out_long$variable == "vid" | out_long$variable == "vid_vn", "vid",
                                ifelse(out_long$variable == "post_vid_fix" | out_long$variable == "fix_rate", "post_vid", "pre_vid")))

out_long$variable <- factor(out_long$variable, levels=c("pre_vid_fix", "mock_v1","mock_v2","mock_v3", "vid_v1", "vid", "vid_vn",  "post_vid_fix","fix_rate" ))
out_long$label <- factor(out_long$label, levels=c("pre_vid", "mock", "vid", "post_vid" ))


out_long$index <- ifelse(out_long$index == "lower2_diag", "off-diagonal [-2]",
                         ifelse(out_long$index == "lower1_diag", "off-diagonal [-1]",
                                ifelse(out_long$index == "diag", "diagonal",
                                       ifelse(out_long$index == "upper1_diag", "off-diagonal [+1]",
                                              ifelse(out_long$index == "upper2_diag", "off-diagonal [+2]",
                                                     NA)))))

out_long$index <- factor(out_long$index, levels = c("off-diagonal [-2]", "off-diagonal [-1]", "diagonal", "off-diagonal [+1]", "off-diagonal [+2]"))
#out_long$index <- factor(out_long$index, levels = c("lower2_diag", "lower1_diag", "diag", "upper1_diag", "upper2_diag"))

out_long$lag <- gsub("label2_lag", "HRF lag = ", out_long$lag)

# plot mean and se
p <- ggplot(data = out_long, aes(x = variable, y = mean_corr, col = label))
p <- p  + geom_hline(yintercept=0, col = "lightgrey", linetype="dashed") + geom_point(stat="identity", size = 1) +
  facet_grid(index ~ lag, scales = "free") +
  theme_bw() +
  geom_errorbar(aes(ymin=mean_corr-se_corr, ymax=mean_corr+se_corr), width=.35, position=position_dodge(.9), color = "darkgray") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95)) +
  theme(legend.position = "none") +
  xlab("Volume") + ylab("mean correlation") + labs(title = "Mean correlation averaged over events colour-coded for category")

print(p)

ggsave(file.path(fig_dir,"average_corr.jpg"), width = 30, height = 30, units = "cm")


### plot diagonal and off diagonal for each lag

sISC_long <- reshape2::melt(sISC, id.vars = c("vol", "label_detail",
                                              "label_lag0", "label2_lag0", "label_lag1", "label2_lag1",
                                              "label_lag2", "label2_lag2", "label_lag3", "label2_lag3",
                                              "label_lag4", "label2_lag4", "label_lag5", "label2_lag5",
                                              "label_lag6", "label2_lag6"))
sISC_long$value <- tanh(sISC_long$value) # back tranform z values into raw correlation values
names(sISC_long)[names(sISC_long) == "value"] <- "cor"


sISC_long$variable <- ifelse(sISC_long$variable == "lower2_diag", "off-diagonal [-2]",
                             ifelse(sISC_long$variable == "lower1_diag", "off-diagonal [-1]",
                                    ifelse(sISC_long$variable == "diag", "diagonal",
                                           ifelse(sISC_long$variable == "upper1_diag", "off-diagonal [+1]",
                                                  ifelse(sISC_long$variable == "upper2_diag", "off-diagonal [+2]",
                                                         NA)))))

sISC_long$variable <- factor(sISC_long$variable, levels = c("off-diagonal [-2]", "off-diagonal [-1]", "diagonal", "off-diagonal [+1]", "off-diagonal [+2]"))



#lags <- c("label_lag0", "label_lag1", "label_lag2", "label_lag3",
#          "label_lag4", "label_lag5", "label_lag6")

lags <- c("label2_lag0", "label2_lag1", "label2_lag2", "label2_lag3",
          "label2_lag4", "label2_lag5", "label2_lag6")

# go into directory where to save the figures
setwd(fig_dir)

for (lag in lags) {

  #sISC_long[,lag] <- factor(sISC_long[, lag], levels=c("pre_vid_fix", "mock", "vid","post_vid_fix" ))

  # declare factor
  sISC_long[,lag] <- factor(sISC_long[,lag], levels=c("pre_vid_fix", "mock_v1","mock_v2","mock_v3", "vid_v1", "vid", "vid_vn",  "post_vid_fix","fix_rate" ))

  # create ggplot
  p <- ggplot(data = sISC_long, aes(x = vol, y = cor, col = get(lag)))
  p <- p + geom_hline(yintercept=0, col = "lightgrey", linetype="dashed") +
    geom_line(col = "darkgrey") + geom_point(size = 1) +
    theme_bw() +
    facet_grid(variable ~ . , scales = "free") +
    scale_x_continuous(breaks = seq(0, 960, by=50), labels = as.character(seq(0, 960, by=50)), expand = expansion(add = 10)) +
    theme(legend.position = "bottom") + guides(col=guide_legend(nrow=1,byrow=TRUE)) +
    labs(title = paste0("spatial ISC in V2 colour-coded for ", gsub("label2_lag", "HRF lag = ", lag)), col = "") + xlab("Volume") + ylab("Correlation")
  print(p)
  # save ggplot
  ggsave(paste0(lag, ".jpg"), width = 50, height = 30, units = "cm")

}

# go back into main folder
setwd(anal_dir)
