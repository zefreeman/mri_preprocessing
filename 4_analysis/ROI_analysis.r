# Read the CSV from Desktop 
dat <- read.csv("/Users/z/Documents/2_STAR_MRI/data/20250415_thesis_checkin_tables.csv", header = TRUE)

# Preview the first two columns
head(dat[, 1:2])

# Fit the linear regression w CAPS_B and amygdala output
fit_a <- lm(dat[[2]] ~ dat[[1]], data = dat)
summary(fit_a)

# Fit the linear regression w CAPS_B and hippocampus output
fit_h <- lm(dat[[3]] ~ dat[[1]], data = dat)
summary(fit_h)

# Fit the linear regression w CAPS_B and vmPFC output
fit_v <- lm(dat[[4]] ~ dat[[1]], data = dat)
summary(fit_v)
