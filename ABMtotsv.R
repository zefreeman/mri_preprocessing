library(tidyr)
library(tidyverse)

# -------------------------------------------------------------------------

#   root of folders for when zsh wasn't setting the wd
#root <- "/Users/z/Documents/STAR/scans"

#   make list of STAR participant numbers for parent folders
#folders <- list.files(path="/Users/z/Documents/STAR/scans")

                    # load in file (safe version)
                    #setwd("/Users/z/Documents/STAR/scans/STARP040181/")
                    #fullcsv<- read.csv(Sys.glob(file.path("T2/task","*ZTRM*")))
  

# loop that I gave up on: -------------------------------------------------

#for (i in 1:length(folders)) {
#  workingd <- (paste(root,"/",folders[i], "/",sep=""))
#  setwd(workingd)

# -------------------------------------------------------------------------


fullcsv<- read.csv(Sys.glob(file.path("T2/task","*ZTRM*")))
  
# chop off the last rows associated with the mood questions and the scan runtime information
abbreviated <- fullcsv[c(1:27),]


############# get onset times #################################################
# take onset columns to then reshape into a single column
allonsets <- abbreviated[,c("TimeAtStartOfTrial","ElaborateCueOnset","PostElaborateFixOnset","RatingOnset","RestJitterOnset")]
# make all columns numeric so you can add the columns together
lapply(allonsets,as.numeric)

# reshape to a single onset times column for the whole run
d1 <- data.frame(a=unlist(allonsets, use.names = FALSE)) #no idea what this does but if I remove, it breaks

d2 <- d1[order(as.numeric(as.character(d1$a))),] #changing the type so they can be ordered by more than the index
as.numeric(d2) # PIPE THIS?

d3 <- d2[c(4:134)] # subsetting to remove the 0s from the end of run which were reordered to beginning
onset <- as.numeric(d3) # PIPE THIS? 


############# make a durations vector #########################################
# get durations as vector by asking the difference between cell 1 and cell 2, etc.
durations <- diff(onset)
durations <- append(durations, 5000) # add the 5s of fixation at the end but before the mood questions to match rows in onset column
#View(data.frame(durations)) #sense check in same format as example tsv file


############# adding a trial_type vector ######################################
# make trial type column with repeating cue, elaborate, pelaborate fix, rating
trial_type <- rep(c("cue", "elaborate", "pelaboratefix", "rating", "rest"),times=26)
trial_type <- append(trial_type,"fixation")


############# making is_own and ratings vectors ###############################
# take 'is own', repeat each line 3 more times
ratingcolumns <- abbreviated[,c("IsOwn","Rating","RatingRT")]
rcolumns <- ratingcolumns[rep(seq_len(nrow(ratingcolumns)), each=5),] #repeat each line 5 times to get is_own column complete
rcolumns<-rcolumns[c(1:131),]
#View(data.frame(rcolumns)) #sense check in same format as example tsv file

booleanrating <- rep(c(FALSE, FALSE, FALSE, TRUE, FALSE),times=26) #boolean of 'is this row a rating section of the task or not?'
booleanrating<-append(booleanrating,FALSE)
rcolumns$ratingrow <- booleanrating

# only allowing reaction times and ratings on rows which pertain to rating phases of the task
# NB: ifelse in base R does not need the two columns to be the same type, but if_else from dplyr does need that
rcolumns <- rcolumns %>% 
  mutate(RatingRT = ifelse(ratingrow == FALSE, "n/a", RatingRT)) 
rcolumns <- rcolumns %>%  
  mutate(Rating = ifelse(ratingrow == FALSE, "n/a", Rating)) 

tsv <- data.frame(onset, durations, trial_type, rcolumns)
tsv <- as_tibble(tsv)


############# NaN and . replaced as n/a #######################################
tsv[tsv == "NaN"] <- "n/a"
tsv[tsv == "."] <- "n/a"


############# convert from ms to s ############################################
tsv$onset<-as.numeric(tsv$onset) / 1000
tsv$durations<-as.numeric(tsv$durations) / 1000
tsv$RatingRT<-as.numeric(tsv$RatingRT) / 1000


############# change column names and tidy formatting #########################
# change names to is_own, rating, and response_time
tsv<-rename(tsv, duration=durations, is_own=IsOwn, rating=Rating, response_time=RatingRT)

# change boolean column to TRUE only when rating actually happened
tsv$ratingrow[tsv$ratingrow == TRUE & tsv$duration<1] <- FALSE


############# output as tsv file ##############################################
write_tsv(tsv, "abm_events.tsv")

