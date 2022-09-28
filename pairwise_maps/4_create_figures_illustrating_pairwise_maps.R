####################################################################################
########################### script figures task paper ##############################
####################################################################################

###### setups ######

library(ggplot2)
library(wesanderson)
library(plot.matrix)

# empty work space, load libraries and functions
rm(list=ls())

# define directories
anal_dir <- getwd() # current directory where the project file is saved
fig_dir <- file.path(anal_dir, "figures")
if (dir.exists(fig_dir) == F ) {
  dir.create(fig_dir)
}
trials_dir <- file.path(anal_dir, "trials_files") # directory including files with onset, duration, and behavioural measures of each trial

# define root directory and go back to anal_dir
setwd('..')
root_dir <- getwd()
setwd(anal_dir)

# define input directoty for concat files
concat_dir <- file.path(root_dir, "concat", "init_concat_V2")


###### figure illustrating volumes selected in initial concatenation ######

timeline <- data.frame(name = c("Fixation pre video", "Mock video (v1)", "Mock video (v2)", "Mock video (v3)",     
                                "Magic trick video (v1)","Magic trick video", "Magic trick video (vn)",
                                "Fixation post video (v1)", "Fixation post video (v2)", "Fixation/Rating", "Rating"),
                       start = c(0, 1, 2, 3,
                                 4, 5, 20, 
                                 21, 22, 23, 28),
                       end = c(1, 2, 3, 4,
                               5, 20, 21,
                               22, 23, 28, 29),
                       Category = c("Pre video", "Mock video", "Mock video", "Mock video", 
                                    "Magic trick", "Magic trick", "Magic trick", 
                                    "Post video", "Post video", "Post video", "Post video"),
                       stringsAsFactors = FALSE)
timeline$Category <- factor(timeline$Category, levels=c("Pre video", "Mock video", "Magic trick", "Post video"))
timeline$mid <- (timeline$start + timeline$end) / 2
volumes <- data.frame(volume = c(0:29),
                      y = rep(0, 30))

# create graph
p <- ggplot(volumes, aes(x=volume,y=y)) +
  geom_rect(data=timeline, aes(NULL,NULL,xmin=start,xmax=end,fill=Category),
            ymin=-Inf,ymax=Inf, size=0.5, alpha = 0.75, col = "white") +
  scale_fill_manual(values=c("#FF0000", "#00A08A", "#F2AD00",  "#5BBCD6")) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0.5, 28.5, by=1), labels = rep("", 29)) + scale_y_discrete() +
  geom_text(data=timeline,aes(x=mid,y=0,label=name, fontface = "bold"), size=3, angle = 90) +
  ylab("") + xlab("Volumes") +
  theme(legend.position = "bottom") +
  annotate(
    geom = "rect", 
    # coordinates for the rectangle
    xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
    color = "black", alpha = 0, size = 1
  ) +
  annotate(
    geom = "rect", 
    # coordinates for the rectangle
    xmin = 21, xmax = 22, ymin = -Inf, ymax = Inf,
    color = "black", alpha = 0, size = 1
  )
  #scale_fill_manual(values=wes_palette(n=4, name="Darjeeling1"))

# print graph
print(p)
# save graph
# ggsave(file.path(fig_dir,"volume_category.jpg"), width = 15, height = 5, units = "cm")

###### figure illustrating initial and final concatenation ######

# read in data to extract mean duration of magic tricks and convert that to volumes #

# read in long data
dataLong <- read.csv("MMC_experimental_data.csv")

# center curiosity:  within cluster (CWC, i.e., group-mean centering)
subjects <- unique(dataLong$BIDS)
dataLong$curiosity_cwc <- NA
for (s in seq_along(subjects)) {
  dataLong[dataLong$BIDS == subjects[s], "curiosity_cwc"] <- dataLong[dataLong$BIDS == subjects[s], "responseCuriosity"] - mean(dataLong[dataLong$BIDS == subjects[s], "responseCuriosity"], na.rm = T)
}


# determine TR
TR = 2

# transform long data into wide data
displayVidDuration <- reshape::cast(dataLong, ID~vidFileName,value="displayVidDuration")

# compute mean for the display duration of each magic trick
mean_displayVidDuration <- as.data.frame(mapply(mean, displayVidDuration[-1]))
names(mean_displayVidDuration) <- "mean_displayVidDuration"
mean_displayVidDuration$stim_file <- row.names(mean_displayVidDuration)
row.names(mean_displayVidDuration) <- NULL

# calculate mean_displayVidDuration without mock
mean_displayVidDuration$mean_displayVidDuration <- mean_displayVidDuration$mean_displayVidDuration - 6 # subtract mock
mean_displayVidDuration$mean_displayVidDuration_TR <- ceiling(round(mean_displayVidDuration$mean_displayVidDuration, digits = 0)/TR)

# duration of magic tricks in volumes (already ordered)
dur_vol <- mean_displayVidDuration$mean_displayVidDuration_TR

# for each magic trick, add their start and end volume in the concatenated timeseries
for (r in 1:nrow(mean_displayVidDuration)) {
  
  if(r == 1){
    
    # initial concat
    mean_displayVidDuration$start_init[r] <- 1
    mean_displayVidDuration$end_init[r] <- mean_displayVidDuration$start_init[r] + mean_displayVidDuration$mean_displayVidDuration_TR[r] + 10 - 1
    
    # final concat
    mean_displayVidDuration$start_vol[r] <- 1
    mean_displayVidDuration$end_vol[r] <- mean_displayVidDuration$start_vol[r] + mean_displayVidDuration$mean_displayVidDuration_TR[r] - 1
    
  } else {
    
    # initial concat
    mean_displayVidDuration$start_init[r] <- mean_displayVidDuration$end_init[r-1] + 1
    mean_displayVidDuration$end_init[r] <- mean_displayVidDuration$start_init[r] + mean_displayVidDuration$mean_displayVidDuration_TR[r] + 10 - 1
    
    # final concat
    mean_displayVidDuration$start_vol[r] <- mean_displayVidDuration$end_vol[r-1] + 1
    mean_displayVidDuration$end_vol[r] <- mean_displayVidDuration$start_vol[r] + mean_displayVidDuration$mean_displayVidDuration_TR[r] - 1
  }
}

# add a short identifier for each magic trick (Trick xx)
mean_displayVidDuration$trick <- c(1:36)
mean_displayVidDuration$trick <- formatC(mean_displayVidDuration$trick, width = 2, format = "d", flag = "0")
mean_displayVidDuration$trick <- paste("Trick", mean_displayVidDuration$trick)

# compute mid point for use as label
mean_displayVidDuration$mid <- (mean_displayVidDuration$start_vol + mean_displayVidDuration$end_vol) / 2
mean_displayVidDuration$mid_init <- (mean_displayVidDuration$start_init + mean_displayVidDuration$end_init) / 2

# read in time series data of each subject #

# subject 1
sub1 <- read.table(file.path(anal_dir,"subj1.txt"))
sub1 <- t(sub1)

# subject 2
sub2 <- read.table(file.path(anal_dir,"subj2.txt"))
sub2 <- t(sub2)

# combine into one df adding a volume column
vol <- 1:1140
whole <- data.frame(vol, sub1)
whole$sub2 <- c(sub2, NA)

# read in onset data of each subject #

# subject 1
onset1 <- read.table(file.path(trials_dir,"sub-control001_task-magictrickwatching_trials.tsv"), header = T)
onset1 <- onset1[,c("stim_file", "mock")]
names(onset1)[names(onset1) == "mock"] <- "mock1"

#subject 2
onset2 <- read.table(file.path(trials_dir,"sub-control002_task-magictrickwatching_trials.tsv"), header = T)
onset2 <- onset2[,c("stim_file", "mock")]
names(onset2)[names(onset2) == "mock"] <- "mock2"

# combine data from both subjects
df_rect <- merge(mean_displayVidDuration, onset1, by = "stim_file")
df_rect <- merge(df_rect, onset2, by = "stim_file")
names(df_rect)[names(df_rect) == "mean_displayVidDuration_TR"] <- "dur_vol"

# compute start and end plus mid point for subject 1
df_rect$start_vol1 <- ceiling(round(df_rect$mock1, digits = 0)/TR) + 1
df_rect$end_vol1 <- df_rect$start_vol1 + df_rect$dur_vol - 1
df_rect$mid1 <- (df_rect$start_vol1 + df_rect$end_vol1) / 2

# also compute start and end for initial concat
df_rect$start_init1 <- df_rect$start_vol1 - 4
df_rect$end_init1 <- df_rect$end_vol1 + 6 
df_rect$mid_init1 <- (df_rect$start_init1 + df_rect$end_init1) / 2

# compute start and end plus mid point for subject 2
df_rect$start_vol2 <- ceiling(round(df_rect$mock2, digits = 0)/TR) + 1
df_rect$end_vol2 <- df_rect$start_vol2 + df_rect$dur_vol - 1
df_rect$mid2 <- (df_rect$start_vol2 + df_rect$end_vol2) / 2

# also compute start and end for initial concat
df_rect$start_init2 <- df_rect$start_vol2 - 4
df_rect$end_init2 <- df_rect$end_vol2 + 6 
df_rect$mid_init2 <- (df_rect$start_init2 + df_rect$end_init2) / 2

# create whole timeseries plots for initial concatenation #

# create plot for subject 1
p1_init <- ggplot(whole, aes(x=vol,y=sub1)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_init1,xmax=end_init1,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid_init1,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  #ylab("Subject A") + theme(axis.title.x=element_blank()) + 
  theme(axis.title=element_blank()) + 
  theme(legend.position = "none") + coord_cartesian(ylim = c(-1,2))
# print plot
print(p1_init)
# save plot
# ggsave(file.path(fig_dir,"timeseries_full_s1.jpg"), width = 15, height = 7, units = "cm")

# create plot for subject 2
p2_init <- ggplot(whole, aes(x=vol,y=sub2)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_init2,xmax=end_init2,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid_init2,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  #ylab("Subject B") + theme(axis.title.x=element_blank()) + 
  theme(axis.title=element_blank()) + 
  theme(legend.position = "none") + coord_cartesian(ylim = c(-1,2)) 
# print plot
print(p2_init)
# save plot
# ggsave(file.path(fig_dir,"timeseries_full_s2.jpg"), width = 15, height = 7, units = "cm")

# create whole timeseries plots #

# create plot for subject 1
p1 <- ggplot(whole, aes(x=vol,y=sub1)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_vol1,xmax=end_vol1,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() + 
  geom_line(col = "#FF2600") +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid1,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  ylab("Subject A") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "bottom") + coord_cartesian(ylim = c(-1.6, 2.6))
# print plot
print(p1)
# save plot
# ggsave(file.path(fig_dir,"timeseries_full_s1.jpg"), width = 15, height = 7, units = "cm")

# create plot for subject 2
p2 <- ggplot(whole, aes(x=vol,y=sub2)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_vol2,xmax=end_vol2,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() + 
  geom_line(col = "#0433FF") +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid2,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 2") +
  ylab("Subject B") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "bottom") + coord_cartesian(ylim = c(-1.6, 2.6))
# print plot
print(p2)
# save plot
# ggsave(file.path(fig_dir,"timeseries_full_s2.jpg"), width = 15, height = 7, units = "cm")


# create labels for volumes for inital and final concatenated timeseries #

# add label to timeseries for subject 1
whole$label_init1 <- NA
whole$label1 <- NA
for (v in 1:nrow(whole)) {
  for (l in 1:nrow(df_rect)) {
    if(whole$vol[v] == df_rect$start_init1[l]){
      while (whole$vol[v] <= df_rect$end_init1[l]) {
        whole$label_init1[v] <- df_rect$trick[l] # add label_init
        if(whole$vol[v] >= df_rect$start_vol1[l] & whole$vol[v] <= df_rect$end_vol1[l]){
          whole$label1[v] <- df_rect$trick[l] # add label
        }
        v <- v + 1
      }
    }
  }
} 

# add label to timeseries for subject 2
whole$label_init2 <- NA
whole$label2 <- NA
for (v in 1:nrow(whole)) {
  for (l in 1:nrow(df_rect)) {
    if(whole$vol[v] == df_rect$start_init2[l]){
      while (whole$vol[v] <= df_rect$end_init2[l]) {
        whole$label_init2[v] <- df_rect$trick[l] # add label_init
        if(whole$vol[v] >= df_rect$start_vol2[l] & whole$vol[v] <= df_rect$end_vol2[l]){
          whole$label2[v] <- df_rect$trick[l] # add label
        }
        v <- v + 1
      }
    }
  }
} 

# create df with initial concatenated volumes only #

# create time series data frames with concatened volumes only subject 1
s1 <- whole[complete.cases(whole$label_init1),]

# create time series data frames with concatened volumes only subject 2
s2 <- whole[complete.cases(whole$label_init2),]

# reorder volumes by magic trick and create new volume counter subject 1
s1 <- s1[with(s1, order(label_init1, vol)),]
s1$vol_new <- c(1:dim(s1)[1])

# reorder volumes by magic trick and create new volume counter subject 2
s2 <- s2[with(s2, order(label_init2, vol)),]
s2$vol_new <- c(1:dim(s2)[1])

# drop columns not needed and merge data from both subjects
s1 <- s1[,c("vol_new", "sub1", "label_init1")]
s2 <- s2[,c("vol_new", "sub2", "label_init2")]
temp <- merge(s1, s2, by = "vol_new")
concat <- data.frame(vol_new = vol)
concat <- temp
rm(temp)

# create plots of initial concatenated time series #

# create plot for subject 1
p1_init_new <- ggplot(concat, aes(x=vol_new,y=sub1)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_init,xmax=end_init,col = "grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid_init,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Concatenated time series Subject 1") +
  #ylab("Subject A") +  theme(axis.title.x=element_blank()) + 
  theme(axis.title=element_blank()) + 
  theme(legend.position = "bottom") + coord_cartesian(ylim = c(-1.6, 2.6))
# print plot
print(p1_init_new)
# save plot
# ggsave(file.path(fig_dir,"timeseries_concat_s1.jpg"), width = 15, height = 7, units = "cm")

# create plot for subject 2
p2_init_new <- ggplot(concat, aes(x=vol_new,y=sub2)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_init,xmax=end_init,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid_init,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Concatenated time series Subject 2") +
  #ylab("Subject B") + theme(axis.title.x=element_blank()) + 
  theme(axis.title=element_blank()) + 
  theme(legend.position = "bottom") + coord_cartesian(ylim = c(-1.6, 2.6))
# print plot
print(p2_init_new)
# save plot
# ggsave(file.path(fig_dir,"timeseries_concat_s2.jpg"), width = 15, height = 7, units = "cm")

# create df with final concatenated volumes only #

# create time series data frames with concatened volumes only subject 1
s1 <- whole[complete.cases(whole$label1),]

# create time series data frames with concatened volumes only subject 2
s2 <- whole[complete.cases(whole$label2),]

# reorder volumes by magic trick and create new volume counter subject 1
s1 <- s1[with(s1, order(label1, vol)),]
s1$vol_new <- c(1:dim(s1)[1])

# reorder volumes by magic trick and create new volume counter subject 2
s2 <- s2[with(s2, order(label2, vol)),]
s2$vol_new <- c(1:dim(s2)[1])

# drop columns not needed and merge data from both subjects
s1 <- s1[,c("vol_new", "sub1", "label1")]
s2 <- s2[,c("vol_new", "sub2", "label2")]
temp <- merge(s1, s2, by = "vol_new")
concat <- data.frame(vol_new = vol)
concat <- temp
rm(temp)

# create plots of final concatenated time series #

# create plot for subject 1
p1_new <- ggplot(concat, aes(x=vol_new,y=sub1)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_vol,xmax=end_vol,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() + 
  geom_line(col = "#FF2600") +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Concatenated time series Subject 1") +
  ylab("Subject A") +  theme(axis.title.x=element_blank()) + 
  theme(legend.position = "bottom") + coord_cartesian(ylim = c(-1.6, 2.6))
# print plot
print(p1_new)
# save plot
# ggsave(file.path(fig_dir,"timeseries_concat_s1.jpg"), width = 15, height = 7, units = "cm")

# create plot for subject 2
p2_new <- ggplot(concat, aes(x=vol_new,y=sub2)) +
  geom_rect(data=df_rect, aes(NULL,NULL,xmin=start_vol,xmax=end_vol,col="grey"),
            ymin=-Inf,ymax=Inf, colour=NA, size=0.5, alpha = 0.5) +
  theme_classic() + 
  geom_line(col = "#0433FF") +
  scale_x_continuous(breaks = seq(0, 1100, by=100), labels = seq(0, 1100, by=100)) + scale_y_discrete() +
  geom_text(data=df_rect,aes(x=mid,y=1.5,label=trick), size=3, angle = 90) +
  #ylab("") + xlab("") + labs(title = "Concatenated time series Subject 2") +
  ylab("Subject B") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "bottom") + coord_cartesian(ylim = c(-1.6, 2.6))
# print plot
print(p2_new)
# save plot
# ggsave(file.path(fig_dir,"timeseries_concat_s2.jpg"), width = 15, height = 7, units = "cm")

# combine all plots #
ggpubr::ggarrange(
                  ggpubr::ggarrange(p1_init, p1_init_new,                  
                  nrow = 2),
                  p,
                  ncol = 2,
                  labels = c("", "")   # Labels of the scatter plot
) 
# save plot: Figure S2A and S2B
ggsave(file.path(fig_dir,"timeseries_init.jpg"), width = 25, height = 8, units = "cm")


# combine all plots #
ggpubr::ggarrange(p1,p2, p1_new, p2_new,                  
                  nrow = 2, ncol = 2, 
                  labels = c("", "", "", "")   # Labels of the scatter plot
) 
# save plot: Figure 2A
ggsave(file.path(fig_dir,"timeseries_concat.jpg"), width = 30, height = 8, units = "cm")

###### figure illustrating spatial ISC ######

# load in data from spatial ISC #
mean_MZ_all <- read.csv(file.path(concat_dir, "mean_cor_z.csv"), row.names = 1)
mean_cor_all <- read.csv(file.path(concat_dir, "mean_cor_raw.csv"), row.names = 1)

mean_MZ_all <- as.matrix(mean_MZ_all)
mean_cor_all <- as.matrix(mean_cor_all)

# create matrix with random numbers between 0 and 1 
matrix <- data.frame(t0_a = runif(9),
                     t1_a = runif(9),
                     t2_a = runif(9),
                     tn_a = runif(9),
                     t0_mean = runif(9),
                     t1_mean = runif(9),
                     t2_mean = runif(9),
                     tn_mean = runif(9))

# plot matrix
png(file.path(fig_dir,"t0_a.png"), 500, 500)
plot(matrix(matrix$t0_a, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"t0_m.png"), 500, 500)
plot(matrix(matrix$t0_mean, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"t1_a.png"), 500, 500)
plot(matrix(matrix$t1_a, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"t1_m.png"), 500, 500)
plot(matrix(matrix$t1_mean, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"t2_a.png"), 500, 500)
plot(matrix(matrix$t2_a, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"t2_m.png"), 500, 500)
plot(matrix(matrix$t2_mean, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"tn_a.png"), 500, 500)
plot(matrix(matrix$tn_a, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# plot matrix
png(file.path(fig_dir,"tn_m.png"), 500, 500)
plot(matrix(matrix$tn_mean, ncol=3), key=NULL, axis.col=NULL, axis.row=NULL, xlab='', ylab='', breaks = 20)
dev.off()

# compute correlation matrix (based on the mean correlation matrix computed in the spatial_ISC.R script) +
cor <- mean_cor_all[3:6, 3:6]
colnames(cor) <- c("Mean t(i)", "Mean t(i + 1)", "Mean t(i + 2)", "Mean t(n)")
row.names(cor) <- c("Subject t(i)", "Subject t(i + 1)", "Subject t(i + 2)", "Subject t(n)")

ggcorrplot::ggcorrplot(cor, outline.col = "white", lab_col = "black", pch.col = "black",  tl.col = "black", ggtheme = ggplot2::theme_bw,
                       colors = c("darkred", "white", "lightblue"),
                       lab = T, show.legend = F)
ggsave(file.path(fig_dir,"spatialISC.jpg"), height = 10, width = 10, units = "cm")

# create correlation matrix #

# create smaller version of cor matrix
small_m <- mean_cor_all[1:100, 1:100]
ggcorrplot::ggcorrplot(small_m, outline.col = "white", lab_col = "black", pch.col = "black",  tl.col = "black", 
                       ggtheme = ggplot2::theme_classic(),
                       #ggtheme = ggplot2::theme_bw(),
                       colors = c("darkred", "white", "lightblue"),
                       #lab = T, lab_size = 2)
                       lab = F, lab_size = 2)
ggsave(file.path(fig_dir,"cor_matrix.jpg"), height = 12, width = 12, units = "cm")


###### figure illustrating label lag ######

# extract upper off diagonal
upper1_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) - 1)]
upper2_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) - 2)]

# create a data frame #
sISC <- data.frame(vol = 1:(length(upper1_diag)+1))

sISC$upper1_diag <- c(upper1_diag, NA) # t(i) and V(i + 1)
sISC$upper2_diag <- c(upper2_diag, NA, NA) # t(i) and V(i + 2)

# create vector with labels #
for (t in 1:length(dur_vol)) {
  if (t == 1){
    labels <- c(rep("pre_vid_fix", 1), rep("mock", 3), rep("vid", dur_vol[t]), rep("post_vid_fix", 6))
    labels2 <- c(rep("pre_vid_fix", 1), "mock_v1","mock_v2","mock_v3", "vid_v1", rep("vid", dur_vol[t]-2), "vid_vn", "post_vid_fix_v1", "post_vid_fix_v2", rep("fix_rate", 3), "rate")   
    labels_detail <- c(rep("pre_vid_fix", 1), rep(paste0("mock_", t), 3), rep(paste0("vid_", t), dur_vol[t]), rep("post_vid_fix", 6))   
  } else {
    labels <-  c(labels, rep("pre_vid_fix", 1), rep("mock", 3), rep("vid", dur_vol[t]), rep("post_vid_fix", 6)) 
    labels2 <- c(labels2, rep("pre_vid_fix", 1), "mock_v1","mock_v2","mock_v3", "vid_v1", rep("vid", dur_vol[t]-2), "vid_vn", "post_vid_fix_v1", "post_vid_fix_v2", rep("fix_rate", 3), "rate")   
    labels_detail <- c(labels_detail, rep("pre_vid_fix", 1), rep(paste0("mock_", t), 3), rep(paste0("vid_", t), dur_vol[t]), rep("post_vid_fix", 6))   
  }
}

# add labels to df
sISC$label_detail <- labels_detail 
sISC$label_lag0 <- labels
sISC$label2_lag0 <- labels2

# add labels assuming different lags
sISC$label_lag4 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-4]))
sISC$label_lag4[1:4] <- NA
sISC$label2_lag4 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label2_lag0[x-4]))
sISC$label2_lag4[1:4] <- NA

# reduce sISC to 100 vol
sISC_long <- sISC[301:400,]

# create data in long format
sISC_long <- reshape2::melt(sISC_long, id.vars = c("vol", "label_detail", 
                                              "label_lag0", "label2_lag0",  
                                              "label_lag4", "label2_lag4"))
sISC_long$value <- tanh(sISC_long$value) # back tranform z values
names(sISC_long)[names(sISC_long) == "value"] <- "cor"

sISC_long$variable <- ifelse(sISC_long$variable == "upper1_diag", "k = 1", 
                                                  ifelse(sISC_long$variable == "upper2_diag", "k = 2", 
                                                         NA))

sISC_long$variable <- factor(sISC_long$variable, levels = c("k = 1", "k = 2"))

# declare factor
sISC_long$label2_lag4 <- factor(sISC_long$label2_lag4, levels=c("pre_vid_fix", "mock_v1","mock_v2","mock_v3", "vid_v1", "vid", "vid_vn",  "post_vid_fix_v1",  "post_vid_fix_v2", "fix_rate", "rate"))
sISC_long$label_lag4 <- factor(sISC_long$label_lag4, levels=c("pre_vid_fix", "mock", "vid",   "post_vid_fix"))

# create ggplot not colour coded
sISC_long$header <- "Extracted values"
sISC_long$vol <- c(1001: 1100)
p_grey <- ggplot(data = sISC_long, aes(x = vol, y = cor))
p_grey <- p_grey + geom_hline(yintercept=0, col = "lightgrey", linetype="dashed") + 
  geom_line(col = "darkgrey") + geom_point(size = 1) +
  theme_bw() + 
  facet_grid(variable ~ header) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95, colour = "white")) +
  scale_x_continuous(breaks = c(1001:1100), labels = seq(1001, 1100, by = 1), expand = expansion(add = 2)) + 
  xlab("Volume") + ylab("Correlation")
# show plot
print(p_grey)
# save plot: Figure S2D
ggsave(file.path(fig_dir,"off_diag_grey.jpg"), width = 10, height = 15, units = "cm")

# create ggplot colour-coded
sISC_long$header <- "Labelled values"
p_col <- ggplot(data = sISC_long, aes(x = vol, y = cor, col = label_lag4))
p_col <- p_col + geom_hline(yintercept=0, col = "lightgrey", linetype="dashed") + 
  geom_line(col = "darkgrey") + geom_point(size = 1, aes(col = label_lag4)) +
  theme_bw() + 
  facet_grid(variable ~ header) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95, colour = "white")) +
  scale_x_continuous(breaks = c(1001:1100), labels = seq(1001, 1100, by = 1), expand = expansion(add = 2)) + 
  xlab("Volume") + ylab("Correlation") +
  scale_color_manual(values=c("#FF0000", "#00A08A", "#F2AD00",  "#5BBCD6")) +
  theme(legend.position = "none")
# show plot
print(p_col)
# save plot: Figure S2D
ggsave(file.path(fig_dir,"off_diag_col.jpg"), width = 10, height = 15, units = "cm")


###### figure showing the results of spatial ISC based labeling ###### 

# extract upper off diagonal
upper1_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) - 1)]
upper2_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) - 2)]

# extract lower off diagonal
lower1_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) + 1)]
lower2_diag <- mean_MZ_all[row(mean_MZ_all) == (col(mean_MZ_all) + 2)]

# create a data frame #
sISC <- data.frame(vol = 1:(length(upper1_diag)+1))

# add off diagonals to data frame and compute change from i+1 to i+2
sISC$upper1_diag <- c(upper1_diag, NA) # t(i) and V(i + 1)
sISC$upper2_diag <- c(upper2_diag, NA, NA) # t(i) and V(i + 2)
sISC$lower1_diag <- c(lower1_diag, NA) # V(t) and t(i + 1)
sISC$lower2_diag <- c(lower2_diag, NA, NA) # V(t) and t(i + 2)

# create vector with labels #
for (t in 1:length(dur_vol)) {
  if (t == 1){
    labels <- c(rep("pre_vid_fix", 1), "mock_v1","mock_v2","mock_v3", "vid_v1", rep("vid", dur_vol[t]-2), "vid_vn", "post_vid_fix_v1", "post_vid_fix_v2", rep("fix_rate", 3), "rate")   
    labels_detail <- c(rep("pre_vid_fix", 1), rep(paste0("mock_", t), 3), rep(paste0("vid_", t), dur_vol[t]), rep("post_vid_fix", 6))   
  } else {
    labels <- c(labels, rep("pre_vid_fix", 1), "mock_v1","mock_v2","mock_v3", "vid_v1", rep("vid", dur_vol[t]-2), "vid_vn", "post_vid_fix_v1", "post_vid_fix_v2", rep("fix_rate", 3), "rate")   
    labels_detail <- c(labels_detail, rep("pre_vid_fix", 1), rep(paste0("mock_", t), 3), rep(paste0("vid_", t), dur_vol[t]), rep("post_vid_fix", 6))   
  }
}

# add labels to df
sISC$label_detail <- labels_detail 
sISC$label_lag0 <- labels

# add labels assuming different lags
sISC$label_lag1 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-1]))
sISC$label_lag1[1] <- NA
sISC$label_lag2 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-2]))
sISC$label_lag2[1:2] <- NA
sISC$label_lag3 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-3]))
sISC$label_lag3[1:3] <- NA
sISC$label_lag4 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-4]))
sISC$label_lag4[1:4] <- NA
sISC$label_lag5 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-5]))
sISC$label_lag5[1:5] <- NA
sISC$label_lag6 <-  as.character(sapply(1:nrow(sISC), function(x) sISC$label_lag0[x-6]))
sISC$label_lag6[1:6] <- NA

# compute average correlation for each segment and lag
index <- c("upper1_diag", "upper2_diag", "lower1_diag", "lower2_diag")

lag <- c("label_lag0", "label_lag1", "label_lag2", "label_lag3", "label_lag4", "label_lag5", "label_lag6")

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

# add label
out_long$label <- ifelse(out_long$variable == "mock_v1" | out_long$variable == "mock_v2" | out_long$variable == "mock_v3", "mock",
                         ifelse(out_long$variable == "vid_v1" | out_long$variable == "vid" | out_long$variable == "vid_vn", "vid",
                                ifelse(out_long$variable == "post_vid_fix_v1" | out_long$variable == "post_vid_fix_v2" | out_long$variable == "fix_rate" | out_long$variable == "rate", "post_vid", "pre_vid")))
out_long$label <- factor(out_long$label, levels=c("pre_vid", "mock", "vid", "post_vid" ))

# add better namings for variables and make them factor
out_long$variable <- factor(out_long$variable, levels=c("pre_vid_fix", "mock_v1","mock_v2","mock_v3", "vid_v1", "vid", "vid_vn",  "post_vid_fix_v1", "post_vid_fix_v2",  "fix_rate", "rate" ))

out_long$lag <- gsub("label2_lag", "HRF lag = ", out_long$lag)


# subset data frame
out_long_key <- subset(out_long, variable == "pre_vid_fix" | variable == "post_vid_fix_v1")

# rename variable
out_long_key$lag <- gsub("label_lag", "HRF lag = ", out_long_key$lag)
out_long_key$variable <- ifelse(out_long_key$variable == "pre_vid_fix", "Pre", 
                              ifelse(out_long_key$variable == "post_vid_fix_v1", "Post", NA))

out_long_key$variable <- factor(out_long_key$variable, levels=c("Pre", "Post"))

out_long_key$diag <- ifelse(grepl("upper", out_long_key$index), "Lower diagonal",
                            ifelse(grepl("lower", out_long_key$index), "Upper diagonal",
                                   NA))
out_long_key$k <- ifelse(grepl("1", out_long_key$index), "1",
                            ifelse(grepl("2", out_long_key$index), "2",
                                   NA))


# plot mean and se
#  scale_fill_manual(values=c("#FF0000", "#00A08A", "#F2AD00",  "#5BBCD6"))

p_key <- ggplot(data = out_long_key, aes(x = variable, y = mean_corr, col = label, shape = k))
p_key <- p_key  + geom_hline(yintercept=0, col = "lightgrey", linetype="dashed") + 
  geom_errorbar(aes(ymin=mean_corr-se_corr, ymax=mean_corr+se_corr), width=.35, position=position_dodge(.9), color = "black") + 
  geom_point(aes(shape = k, color = label), stat="identity", size = 2, position=position_dodge(.9), alpha = 0.75) + 
  #scale_shape_manual(values=c(1, 2)) +
  scale_color_manual(values = c("#FF0000", "#5BBCD6"), guide = F) +
  facet_grid(diag ~ lag) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95)) +
  #theme(legend.position = "none") +
  xlab("Fixation volume") + ylab("Mean correlation") #+ labs(title = "Mean correlation during critical fixation volumes")
print(p_key)
# save file: Figure S2E
ggsave(file.path(fig_dir,"pre_post_fix.jpg"), width = 20, height = 15, units = "cm")

cor(out_long_key$mean_corr_z[out_long_key$diag == "Upper diagonal"], out_long_key$mean_corr_z[out_long_key$diag == "Lower diagonal"])

###### figure illustrating ISC map ######

# plot concatenated time series subject 1
subj_a_i <- ggplot(data = concat, aes(x = vol_new, y = sub1))
subj_a_i + geom_line(col = "#FF2600", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a_i.jpg"), width = 4, height = 1, units = "cm")


# plot concatenated time series subject 1
subj_b_i <- ggplot(data = concat, aes(x = vol_new, y = sub2))
subj_b_i + geom_line(col = "#0433FF", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_b_i.jpg"), width = 4, height = 1, units = "cm")

# plot concatenated time series subject 1 + subject 2
concat_l <- reshape2::melt(concat, id.vars = c("vol_new", "label1", "label2"))
subj_ab_i <- ggplot(data = concat_l, aes(x = vol_new, y = value, group = variable, colour = variable))  + ylim(c(-1.6, 2.6))
subj_ab_i + geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#FF2600", "#0433FF")) + theme(legend.position = "none") 
ggsave(file.path(fig_dir,"subj_ab_i.jpg"), width = 4, height = 1, units = "cm")

###### figure illustrating ROI2ROI ISFC ######

# read in ROI data #
ROI <- read.csv(file.path(anal_dir,"subject_pair_roi.csv"))
ROI$vol <- c(1:dim(ROI)[1])

# create individual timeseries plots #

# subject 1 roi 1
subj_a_roi1 <- ggplot(data = ROI, aes(x = vol, y = sub_a_roi1))
subj_a_roi1 + geom_line(col = "#FF2600", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a_roi1.jpg"), width = 4, height = 1, units = "cm")

# subject 1 roi 2
subj_a_roi2 <- ggplot(data = ROI, aes(x = vol, y = sub_a_roi2))
subj_a_roi2 + geom_line(col = "#011893", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a_roi2.jpg"), width = 4, height = 1, units = "cm")

# subject 2 roi 1
subj_b_roi1 <- ggplot(data = ROI, aes(x = vol, y = sub_b_roi1))
subj_b_roi1 + geom_line(col = "#C00000", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_b_roi1.jpg"), width = 4, height = 1, units = "cm")

# subject 2 roi 2
subj_b_roi2 <- ggplot(data = ROI, aes(x = vol, y = sub_a_roi2))
subj_b_roi2 + geom_line(col = "#0432FF", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_b_roi2.jpg"), width = 4, height = 1, units = "cm")

# create combined timeseries plots #

# data roi 1
roi_pair <- ROI[,c("vol", "sub_a_roi1", "sub_b_roi2")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

# plot roi 1
subj_a1b2 <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_a1b2 +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#FF2600", "#0432FF")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a1b2.jpg"), width = 4, height = 1, units = "cm")

# data roi 2
roi_pair <- ROI[,c("vol", "sub_a_roi1", "sub_b_roi2")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

# plot roi 2
subj_a2b1 <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_a2b1 +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#011893", "#C00000")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a2b1.jpg"), width = 4, height = 1, units = "cm")

rm(roi_pair, roi_pair_l)

###### figure illustrating seed ISFC ######

# note: same ROI 1 as above

# read in voxel data
sample <- read.csv(file.path(anal_dir,"subject_sample_vis.csv"), header = F)
sample <- data.frame(t(sample))
sample$vol <- c(1:dim(sample)[1])

# create individual timeseries plots #

# subject b voxel x
subj_b_x <- ggplot(data = sample, aes(x = vol, y = X7))
subj_b_x + geom_line(col = "#5B9BD5", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_b_x.jpg"), width = 4, height = 1, units = "cm")

# subject b voxel y
subj_b_y <- ggplot(data = sample, aes(x = vol, y = X8))
subj_b_y + geom_line(col = "#2E75B6", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_b_y.jpg"), width = 4, height = 1, units = "cm")

# subject b voxel z
subj_b_z <- ggplot(data = sample, aes(x = vol, y = X9))
subj_b_z + geom_line(col = "#BDD7EE", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_b_z.jpg"), width = 4, height = 1, units = "cm")

# subject a voxel x
subj_a_x <- ggplot(data = sample, aes(x = vol, y = X10))
subj_a_x + geom_line(col = "#4472C4", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a_x.jpg"), width = 4, height = 1, units = "cm")

# subject a voxel y
subj_a_y <- ggplot(data = sample, aes(x = vol, y = X11))
subj_a_y + geom_line(col = "#2F5597", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a_y.jpg"), width = 4, height = 1, units = "cm")

# subject a voxel z
subj_a_z <- ggplot(data = sample, aes(x = vol, y = X12))
subj_a_z + geom_line(col = "#B4C7E7", size=0.2) + theme_void() + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_a_z.jpg"), width = 4, height = 1, units = "cm")

# create combined timeseries plots #

# combine roi and sample
df <- merge(ROI, sample)

# roi a voxel bx
roi_pair <- df[,c("vol", "sub_a_roi1", "X7")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

subj_abx <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_abx +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#FF2600", "#5B9BD5")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_abx.jpg"), width = 4, height = 1, units = "cm")

# roi a voxel by
roi_pair <- df[,c("vol", "sub_a_roi1", "X8")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

subj_aby <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_aby +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#FF2600", "#2E75B6")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_aby.jpg"), width = 4, height = 1, units = "cm")

# roi a voxel bz
roi_pair <- df[,c("vol", "sub_a_roi1", "X9")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

subj_abz <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_abz +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#FF2600", "#BDD7EE")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_abz.jpg"), width = 4, height = 1, units = "cm")

# roi b voxel ax
roi_pair <- df[,c("vol", "sub_b_roi1", "X10")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

subj_bax <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_bax +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#C00000", "#4472C4")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_bax.jpg"), width = 4, height = 1, units = "cm")

# roi b voxel ay
roi_pair <- df[,c("vol", "sub_b_roi1", "X11")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

subj_bay <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_bay +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#C00000", "#2F5597")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_bay.jpg"), width = 4, height = 1, units = "cm")

# roi b voxel az
roi_pair <- df[,c("vol", "sub_b_roi1", "X12")]
roi_pair_l <- reshape2::melt(roi_pair, id.vars = "vol")

subj_baz <- ggplot(data = roi_pair_l, aes(x = vol, y = value, group = variable, colour = variable))
subj_baz +  geom_line(size=0.2) + theme_void() + scale_color_manual(values=c("#C00000", "#B4C7E7")) + theme(legend.position = "none") + ylim(c(-1.6, 2.6))
ggsave(file.path(fig_dir,"subj_baz.jpg"), width = 4, height = 1, units = "cm")


###### figure illustrating similarity matrices ######

# brain data #

# read in data and rename columns
group <- read.csv(file.path(anal_dir,"subject_group_vis.csv"), header = F)
group <- data.frame(t(group))
names(group) <- c("S1","S2", "S3", "S4", "S5", "Sn")

# multiply values and create correlation matrix
vis_cor <- 3*round(cor(group), 2)
ggcorrplot::ggcorrplot(vis_cor, type = "lower",
                       outline.col = "white", show.diag = T,
                       colors = c("blue", "white", "red"), hc.order = F, show.legend = F)
# save plot
ggsave(file.path(fig_dir,"brain_cor.jpg"), height = 6, width = 6, units = "cm")

# behavioural data #

# transform curiosity data to wide format
cur_long <- dataLong[,c("ID", "stimID", "curiosity_cwc")]
cur_wide <- reshape2::dcast(cur_long, stimID ~ ID)
cur_wide$stimID <- NULL
cur_wide <- cur_wide[,1:6]
names(cur_wide) <- c("S1","S2", "S3", "S4", "S5", "Sn")

# create correlation matrix, multiply and plot
cur_cor <- 2*round(cor(cur_wide), 2)
ggcorrplot::ggcorrplot(cur_cor, type = "lower",
                       outline.col = "white", show.diag = T,
                       colors = c("blue", "white", "red"), show.legend = F)
# save plot
ggsave(file.path(fig_dir,"behav_cor.jpg"), height = 6, width = 6, units = "cm")

###### figure illustrating behavioural similarities ######
 
# extract behavioural data #

# subject a
data_a <- subset(dataLong, ID == 1) # select data from sub 1
data_a <- data_a[with(data_a, order(vidFileName)),] # order by video
data_a <- data_a[,c("vidFileName", "curiosity_cwc", "recognition")] # select relevant columns
names(data_a) <- c("Video", "rating (A)", "encoding (A)")

# subject b
data_b <- subset(dataLong, ID == 2) # select data from sub 2
data_b <- data_b[with(data_b, order(vidFileName)),] # order by video
data_b <- data_b[,c("vidFileName", "curiosity_cwc", "recognition")] # select relevant columns
names(data_b) <- c("Video", "rating (B)", "encoding (B)")

# merge data from both subjects
data <- merge(data_a, data_b, by = "Video")

# create data in long format
data <- reshape2::melt(data, id.vars = "Video") 

# create data for curiosity, memory, and CMlE
data_cur <- data[c(grepl("rating", data$variable)), ]
data_mem <- data[c(grepl("encoding", data$variable)), ]
data_a <- data[c(grepl("(A)", data$variable)), ]
data_b <- data[c(grepl("(B)", data$variable)), ]

# create plots #

# curiosity timeseries
cur <- ggplot(data_cur, aes(x=Video,y=value, group = variable, col = variable)) +
  theme_void() + 
  geom_line(linetype = "dashed") +
  #scale_x_continuous(breaks = seq(1, 36, by=1), labels = seq(1, 36, by=1)) + #scale_y_discrete() +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  ylab("Curiosity (mean centered)") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "none") + coord_cartesian(ylim = c(-4, 4)) + theme(legend.title = element_blank()) +
  scale_color_manual(values=c("#FF2600", "#0433FF"))
# print plot
print(cur)
# save plot 
ggsave(file.path(fig_dir,"timeseries_cur.jpg"), width = 12, height = 3, units = "cm")


# encoding timeseries
mem <- ggplot(data_mem, aes(x=Video,y=value, group = variable, col = variable)) +
  theme_void() + 
  geom_line(linetype = "solid") +
  #scale_x_continuous(breaks = seq(1, 36, by=1), labels = seq(1, 36, by=1)) + scale_y_discrete() +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  ylab("Encoding (dummy coded)") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "none") + coord_cartesian(ylim = c(-4, 4)) + theme(legend.title = element_blank()) +
  scale_color_manual(values=c("#FF2600", "#0433FF"))
# print plot
print(mem)
# save plot 
ggsave(file.path(fig_dir,"timeseries_mem.jpg"), width = 12, height = 3, units = "cm")

# CMLE in subject a
a <- ggplot(data_a, aes(x=Video,y=value, group = variable)) +
  theme_void() + 
  geom_line(aes(linetype=variable), col = "#FF2600") + scale_linetype_manual(values=c("dashed", "solid")) +
  #geom_line(col = "#FF2600") +
  #geom_point(aes(shape=variable), col = "#FF2600") + scale_shape_manual(values = c(18,20))+
  #scale_x_continuous(breaks = seq(1, 36, by=1), labels = seq(1, 36, by=1)) + scale_y_discrete() +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  ylab("") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "none") + coord_cartesian(ylim = c(-4, 4)) + theme(legend.title = element_blank())
# print plot
print(a)
# save plot 
ggsave(file.path(fig_dir,"subj_a_cmle.jpg"), width = 12, height = 3, units = "cm")


# CMLE in subject b
b <- ggplot(data_b, aes(x=Video,y=value, group = variable)) +
  theme_void() + 
  geom_line(aes(linetype=variable), col = "#0433FF") + scale_linetype_manual(values=c("dashed", "solid")) +
  #scale_x_continuous(breaks = seq(1, 36, by=1), labels = seq(1, 36, by=1)) + scale_y_discrete() +
  #ylab("") + xlab("") + labs(title = "Volumes covering the magic tricks in Subject 1") +
  ylab("") + theme(axis.title.x=element_blank()) + 
  theme(legend.position = "none") + coord_cartesian(ylim = c(-4, 4)) + theme(legend.title = element_blank())
# print plot
print(b)
# save plot 
ggsave(file.path(fig_dir,"subj_b_cmle.jpg"), width = 12, height = 3, units = "cm")





