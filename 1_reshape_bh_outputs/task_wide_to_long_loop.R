## TRANSFORM TASK FILES FMRI FROM WIDE TO LONG (looped) - STAR TRIAL - EME STUDY

# Install and load packages (need writexls if writexl?)
install.packages("stringr")
install.packages("reshape2")
install.packages("tidyr")
install.packages("WriteXLS")
install.packages("readxl")
install.packages("dplyr")
install.packages("writexl")
library(stringr)
library(reshape2)
library(tidyr)
library(WriteXLS)
library(readxl)
library(dplyr)
library(writexl)

# Set the working directory to where the Excel files are saved
setwd("/Users/Inez/VU/RMCDP/Year 2/Research Project 2/STAR/EME")

# Define the number of runs
runs <- c("run1", "run2", "run3")

# Iterate over each run
for (run in runs) {
  # Load the Excel file into R
  file_pattern <- paste0("*", run, "_P040255_2.csv")
  file_list <- list.files(pattern = file_pattern)
  
  # Loop over the found files
  for (file_name in file_list) {
    #Load csv file
    TR <- read.csv(file_name)
    # Remove bottom couple lines with general task info
    TR <- TR[!is.na(TR$TimeAtStartOfTrial),]
  
    # Convert task ended onset time from s to ms (will change back to ms in long dataformat)
    TR$TimeAtStartOfTrial <- ifelse(TR$TrialNo == "Task ended", TR$TimeAtStartOfTrial * 1000, TR$TimeAtStartOfTrial)
  
    # Rename columns to correct format for BIDS
    names(TR)[names(TR) == "IsOwn"] <- "is_own"
    names(TR)[names(TR) == "RatingRT"] <- "response_time"
    names(TR)[names(TR) == "Rating"] <- "rating"
  
    # Wide to long transformation
    long_TR <- melt(TR, measure.vars = c('TimeAtStartOfTrial', 'ElaborateCueOnset', 'PostElaborateFixOnset', 'RatingOnset', 'GroundCueOnset', 'RestJitterOnset'),
                  variable.name = 'trial_type',
                  value.name = 'onset')
  
    # Rename variables within new column corresponding to BIDS format
    levels(long_TR$trial_type)[levels(long_TR$trial_type) == 'TimeAtStartOfTrial'] <- 'cue'
    levels(long_TR$trial_type)[levels(long_TR$trial_type) == 'ElaborateCueOnset'] <- 'elaborate'
    levels(long_TR$trial_type)[levels(long_TR$trial_type) == 'PostElaborateFixOnset'] <- 'pelaboratefix'
    levels(long_TR$trial_type)[levels(long_TR$trial_type) == 'RatingOnset'] <- 'rating'
    levels(long_TR$trial_type)[levels(long_TR$trial_type) == 'GroundCueOnset'] <- 'groundcue'
    levels(long_TR$trial_type)[levels(long_TR$trial_type) == 'RestJitterOnset'] <- 'rest'
  
    # Display response time and rating only for rating trials
    long_TR$response_time[long_TR$trial_type != 'rating'] <- 'n/a'
    long_TR$rating[long_TR$trial_type != 'rating'] <- 'n/a'
    
    # Create a new folder for the SubjectID if it doesn't exist
    subject_folder <- unique(long_TR$SubjectID)
    if (!dir.exists(subject_folder)) {
      dir.create(subject_folder)
    }
    
    # Set the working directory to the subject folder
    setwd(subject_folder)

    # Save long dataset into a new Excel file
    file_name_big <- paste0("long_", long_TR$SubjectID, "_", long_TR$Visit, "_", run, "_big.xlsx")
    write_xlsx(long_TR, file_name_big)
  
    # Delete rows with unnecessary task information
    long_TR<-long_TR[- grep("Fixation", long_TR$TrialType),]
    long_TR<-long_TR[- grep("OtherRating", long_TR$TrialType),]
    long_TR$Visit<-as.numeric(long_TR$Visit)
  
    # Change onset from ms to s
    long_TR$onset<-as.numeric(long_TR$onset)
    long_TR$onset<-long_TR$onset/1000
    long_TR<-long_TR[!is.na(long_TR$onset),]
    
    # Create and save the dataset with necessary columns only
    long_TR_small <- subset(long_TR, select = c(SubjectID, Visit, onset, trial_type, is_own, response_time, rating))
    file_name_small <- paste0("long_", long_TR_small$SubjectID, "_", long_TR_small$Visit, "_", run, "_small", ".xlsx")
    write_xlsx(long_TR_small, file_name_small)
    
    # Set the working directory back to the original directory
    setwd("..")
  }  
}

## To do in Excel manually for 'small' files:

#1. Sort 'onset' column from smallest to largest and expand for whole dataset 
#2. Create new column to the right of 'onset' column and label 'duration'
#3. Calculate duration between cue onset by subtracting the following onset cell by the preceding onset cell
#4. Drag cells down to calculate duration between cues for all cells
