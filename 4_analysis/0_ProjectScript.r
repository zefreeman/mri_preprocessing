##################################################################################
#                                                                                #
#                    Baseline STAR participants fMRI analysis                    #
#                                                                                #
#                         Ze Freeman                                             #
#                                                                                #
#                          Last updated: 17 Jan 2024                             #
#                                                                                #
##################################################################################

library(ARTofR)
library(dplyr)
library(stringr)
library(ggplot2)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          1. load packages & data                         ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
home <- Sys.getenv("HOME") # use "USERPROFILE" if on windows?
datapath <- "Documents/STAR_baseline_fmri/data/"
fullpath <- file.path(home, datapath)

Ldat <- read.delim(paste0(fullpath, "L_group_bold.tsv"))
Ndat <- read.delim(paste0(fullpath, "N_group_bold.tsv"))
Mdat <- read.delim(paste0(fullpath, "M_group_bold.tsv"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                2. clean data                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ldat <- Ldat %>%
  mutate(site = "L", 
        ID = substr(bids_name, 1, 8), 
        task = sub('.*task-(.*?)_.*', '\\1', bids_name)) %>%
  select(ID, site, task, everything())

Ndat <- Ndat %>%
  mutate(site = "N", 
        ID = substr(bids_name, 1, 8), 
        task = sub('.*task-(.*?)_.*', '\\1', bids_name)) %>%
  select(ID, site, task, everything())

Mdat <- Mdat %>%
  mutate(site = "M", 
        ID = substr(bids_name, 1, 8), 
        task = sub('.*task-(.*?)_.*', '\\1', bids_name)) %>%
  select(ID, site, task, everything())

fulldat <- Ldat %>%
  bind_rows(Ndat[-1, ], Mdat[-1, ])
write.csv(fulldat, "~/Documents/STAR_baseline_fMRI/data/processed/allsites_mriqc.csv", row.names=FALSE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  3. plots                                ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dattrm1 <- fulldat %>%
  filter(str_detect(task, "TRMrun1"))
write.csv(dattrm1, "~/Documents/STAR_baseline_fMRI/data/processed/TRM1_mriqc.csv", row.names=FALSE)

dattrm2 <- fulldat %>%
  filter(str_detect(task, "TRMrun2"))
write.csv(dattrm2, "~/Documents/STAR_baseline_fMRI/data/processed/TRM2_mriqc.csv", row.names=FALSE)

dattrm3 <- fulldat %>%
  filter(str_detect(task, "TRMrun3"))
write.csv(dattrm3, "~/Documents/STAR_baseline_fMRI/data/processed/TRM3_mriqc.csv", row.names=FALSE)

summary(filtereddat3)

# Create box plot
TRM1plot <- ggplot(filtereddat, aes(x = site, y = fd_mean, fill = site)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.3), size = 2, alpha = 0.7) +
  labs(title = "Box Plot of FD mean for TRMrun1",
       x = "Site",
       y = "Framewise diaplacement mean per participant",
       fill = "Site") +
  ylim(0.5, max(filtereddat$fd_mean) + 0.1) +  
  theme_minimal()
ggsave("~/Documents/STAR_baseline_fMRI/data/processed/boxplot_trmrun1.jpg", plot = TRM1plot, device = "jpg", width = 8, height = 6, units = "in")


### I want to plot the scores 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  4. tables                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





fulldatsummary <- fulldat %>%
  group_by(site, task = sub("^.*?(TRMrun[0-9]+).*", "\\1", task)) %>%
  summarise(mean_fd = mean(fd_mean), count = n())

#   site  task              mean_fd  count
# 1 L     TRMrun1           0.379    38
# 2 L     TRMrun2           0.403    37
# 3 L     TRMrun3           0.414    36
# 4 L     postresting       0.330    36
# 5 L     preresting        0.281    39
# 6 L     socialthreat      0.425    33
# 7 L     socialthreatold   0.430     1
# 8 M     TRMrun1           0.240    13
# 9 M     TRMrun2           0.257    13
#10 M     TRMrun3           0.245    11
#11 M     postresting       0.266    14
#12 M     preresting        0.282    20
#13 M     socialthreat      0.420    14
#14 N     TRMrun1           0.328    12
#15 N     TRMrun2           0.305    12
#16 N     TRMrun3           0.301    12
#17 N     postresting       0.278    14
#18 N     preresting        0.276    14
#19 N     socialthreat      0.345    14



########### explanation of '.*task-(.*?)_.*': for interest ###########
# .*: Matches any character (except for a newline character) zero or more times.
# task-: Matches the literal characters "task-".
# (.*?): This is a non-greedy capture group ((.*?)), which captures any characters (even an empty string) between "task-" and the next underscore "_". 
#        The non-greedy qualifier ? ensures that it captures the shortest possible match.
# _.*: Matches an underscore "_" followed by any characters zero or more times.
# \\1: In the replacement part of sub, \\1 refers to the content captured by the first (and only) capture group (.*?) in the pattern.
#######################################################################
