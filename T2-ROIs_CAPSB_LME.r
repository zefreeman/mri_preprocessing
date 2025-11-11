library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(broom.mixed)
library(purrr)

data_T2 <- read.csv("/Users/zefreeman/Documents/Prepost_ROIs_20251023.csv", header = TRUE, stringsAsFactors = FALSE)

colnames(data_T2)[1] <- "subject"
colnames(data_T2)[3] <- "site"
colnames(data_T2)[5] <- "treatment_group"
colnames(data_T2)[8] <- "age"
colnames(data_T2)[9] <- "CAPSB"
colnames(data_T2)[10] <- "run1vv"
colnames(data_T2)[11] <- "run2vv"
colnames(data_T2)[12] <- "run3vv"
colnames(data_T2)[13] <- "vmpfc_v_total"
colnames(data_T2)[22] <- "run1vd" # not using
colnames(data_T2)[23] <- "run2vd" # not using
colnames(data_T2)[24] <- "run3vd" # not using
colnames(data_T2)[25] <- "vmpfc_d_total" # not using
colnames(data_T2)[34] <- "run1al"
colnames(data_T2)[35] <- "run2al"
colnames(data_T2)[36] <- "run3al"
colnames(data_T2)[37] <- "amyg_l_total"
colnames(data_T2)[46] <- "run1ar"
colnames(data_T2)[47] <- "run2ar"
colnames(data_T2)[48] <- "run3ar"
colnames(data_T2)[49] <- "amyg_r_total"
colnames(data_T2)[58] <- "run1hl"
colnames(data_T2)[59] <- "run2hl"
colnames(data_T2)[60] <- "run3hl"
colnames(data_T2)[61] <- "hipp_l_total"
colnames(data_T2)[70] <- "run1hr"
colnames(data_T2)[71] <- "run2hr"
colnames(data_T2)[72] <- "run3hr"
colnames(data_T2)[73] <- "hipp_r_total"
colnames(data_T2)[90] <- "avo_1"
colnames(data_T2)[91] <- "diss_1"
colnames(data_T2)[92] <- "anx_1"
colnames(data_T2)[93] <- "avo_2"
colnames(data_T2)[94] <- "diss_2"
colnames(data_T2)[95] <- "anx_2"
colnames(data_T2)[96] <- "avo_3"
colnames(data_T2)[97] <- "diss_3"
colnames(data_T2)[98] <- "anx_3"
colnames(data_T2)[101] <- "now_own_1"
colnames(data_T2)[102] <- "now_own_2"
colnames(data_T2)[103] <- "now_own_3"
colnames(data_T2)[104] <- "now_own_mean" # mean own memory nownness score across runs
colnames(data_T2)[108] <- "run1vv_T2"
colnames(data_T2)[109] <- "run2vv_T2"
colnames(data_T2)[110] <- "run3vv_T2"
colnames(data_T2)[112] <- "vmpfc_v_total_T2"
# colnames(data_T2)[113] <- "run1vd_T2" # not using
# colnames(data_T2)[114] <- "run2vd_T2" # not using
# colnames(data_T2)[115] <- "run3vd_T2" # not using
# colnames(data_T2)[116] <- "vmpfc_d_total_T2" # not using
colnames(data_T2)[132] <- "run1al_T2"
colnames(data_T2)[133] <- "run2al_T2"
colnames(data_T2)[134] <- "run3al_T2"
colnames(data_T2)[135] <- "amyg_l_total_T2"
colnames(data_T2)[144] <- "run1ar_T2"
colnames(data_T2)[145] <- "run2ar_T2"
colnames(data_T2)[146] <- "run3ar_T2"
colnames(data_T2)[147] <- "amyg_r_total_T2"
colnames(data_T2)[156] <- "run1hl_T2"
colnames(data_T2)[157] <- "run2hl_T2"
colnames(data_T2)[158] <- "run3hl_T2"
colnames(data_T2)[159] <- "hipp_l_total_T2"
colnames(data_T2)[168] <- "run1hr_T2"
colnames(data_T2)[169] <- "run2hr_T2"
colnames(data_T2)[170] <- "run3hr_T2"
colnames(data_T2)[171] <- "hipp_r_total_T2"
colnames(data_T2)[184] <- "run1hl_beta_own" # hippocampus left beta value to own cue in run 1
colnames(data_T2)[185] <- "run2hl_beta_own"
colnames(data_T2)[186] <- "run3hl_beta_own"
colnames(data_T2)[187] <- "run1hl_beta_other" # hippocampus left beta value to other cue in run 1
colnames(data_T2)[188] <- "run2hl_beta_other"
colnames(data_T2)[189] <- "run3hl_beta_other"
colnames(data_T2)[190] <- "run1hr_beta_own" # hippocampus right beta value to own cue in run 1
colnames(data_T2)[191] <- "run2hr_beta_own"
colnames(data_T2)[192] <- "run3hr_beta_own"
colnames(data_T2)[193] <- "run1hr_beta_other" # hippocampus right beta value to other cue in run 1
colnames(data_T2)[194] <- "run2hr_beta_other"
colnames(data_T2)[195] <- "run3hr_beta_other"
colnames(data_T2)[196] <- "run1whbr_beta_own" # whole brain beta value to own cue in run 1 at baseline
colnames(data_T2)[197] <- "run2whbr_beta_own"
colnames(data_T2)[198] <- "run3whbr_beta_own"
colnames(data_T2)[199] <- "whbr_own_beta_total"
colnames(data_T2)[200] <- "run1whbr_beta_own_T2" # whole brain beta value to own cue in run 1 at follow up
colnames(data_T2)[201] <- "run2whbr_beta_own_T2"
colnames(data_T2)[202] <- "run3whbr_beta_own_T2"
colnames(data_T2)[203] <- "whbr_own_beta_total_T2"
colnames(data_T2)[204] <- "run1hl_beta_own_T2"
colnames(data_T2)[205] <- "run2hl_beta_own_T2"
colnames(data_T2)[206] <- "run3hl_beta_own_T2"
colnames(data_T2)[207] <- "run1hl_beta_other_T2"
colnames(data_T2)[208] <- "run2hl_beta_other_T2"
colnames(data_T2)[209] <- "run3hl_beta_other_T2"
colnames(data_T2)[210] <- "run1hr_beta_own_T2"
colnames(data_T2)[211] <- "run2hr_beta_own_T2"
colnames(data_T2)[212] <- "run3hr_beta_own_T2"
colnames(data_T2)[213] <- "run1hr_beta_other_T2"
colnames(data_T2)[214] <- "run2hr_beta_other_T2"
colnames(data_T2)[215] <- "run3hr_beta_other_T2"

data_T2 <- data_T2 %>%
  mutate(
    # vmPFC v
    run1vv_ch = run1vv_T2 - run1vv,
    run2vv_ch = run2vv_T2 - run2vv,
    run3vv_ch = run3vv_T2 - run3vv,
    # vmPFC d
    # run1vd_ch = run1vd_T2 - run1vd,
    # run2vd_ch = run2vd_T2 - run2vd,
    # run3vd_ch = run3vd_T2 - run3vd,
    # Amygdala left
    run1al_ch = run1al_T2 - run1al,
    run2al_ch = run2al_T2 - run2al,
    run3al_ch = run3al_T2 - run3al,
    # Amygdala right
    run1ar_ch = run1ar_T2 - run1ar,
    run2ar_ch = run2ar_T2 - run2ar,
    run3ar_ch = run3ar_T2 - run3ar,
    # Hippocampus left
    run1hl_ch = run1hl_T2 - run1hl,
    run2hl_ch = run2hl_T2 - run2hl,
    run3hl_ch = run3hl_T2 - run3hl,
    # Hippocampus right
    run1hr_ch = run1hr_T2 - run1hr,
    run2hr_ch = run2hr_T2 - run2hr,
    run3hr_ch = run3hr_T2 - run3hr
  )
# does avoidance score explain variation in ROI signal?
# wide to long including T1 and T2 change scores!
data_long_runs_T2 <- data_T2 %>%
  pivot_longer(
    cols = c(run1vv_ch, run2vv_ch, run3vv_ch, 
             run1al_ch, run2al_ch, run3al_ch, 
             run1ar_ch, run2ar_ch, run3ar_ch, 
             run1hl_ch, run2hl_ch, run3hl_ch, 
             run1hr_ch, run2hr_ch, run3hr_ch),
            #  run1vv_T2, run2vv_T2, run3vv_T2,
            #  run1al_T2, run2al_T2, run3al_T2,
            #  run1ar_T2, run2ar_T2, run3ar_T2,
            #  run1hl_T2, run2hl_T2, run3hl_T2,
            #  run1hr_T2, run2hr_T2, run3hr_T2),
    names_to = c("run", "ROI"),
    names_pattern = "run([0-9]+)([a-zA-Z]+)",
    values_to = "value"
  ) %>%
  mutate(
    run = as.factor(run)
    # timepoint = ifelse(is.na(timepoint), 1, 2)
    # ROI = case_when(
    #   ROI %in% c("vv") ~ "vmpfc_v",
    #   ROI %in% c("al","ar") ~ "amygdala",
    #   ROI %in% c("hl","hr") ~ "hippocampus"
    # )
  ) %>%
  select(subject, site, CAPS.POST.PRE, caps_b_t3, age, run, ROI, value, now_own_mean, STAR_ID)


# ------------------------------------------------------------
# # wide to long for avoidance scores
# data_long_avoid <- data %>%
#   pivot_longer(
#     cols = c(avo_1, avo_2, avo_3),
#     names_to = "run",
#     names_pattern = "avo_([0-9]+)",
#     values_to = "avoid_score"
#   ) %>%
#     mutate(run = as.character(run)) %>%
#   select(subject, run, avoid_score)

# # Make sure 'run' columns are character
# data_long_runs <- data_long_runs %>%
#   mutate(run = as.character(run))

# data_long_avoid <- data_long_avoid %>%
#   mutate(run = as.character(run)) %>%
#   distinct(subject, run, .keep_all = TRUE)  # remove any duplicates

# # Left join: replicate avoid_score across all ROI rows per subject/run
# data_long_runs <- data_long_runs %>%
#   left_join(data_long_avoid, by = c("subject", "run"))


# # wide to long for anxiety
# data_long_anxiety <- data %>%
#   pivot_longer(
#     cols = c(anx_1, anx_2, anx_3),
#     names_to = "run",
#     names_pattern = "anx_([0-9]+)",
#     values_to = "anxiety_score"
#   ) %>%
#     mutate(run = as.character(run)) %>%
#   select(subject, run, anxiety_score)

# # wide to long for dissociation
# data_long_diss <- data %>%
#   pivot_longer(
#     cols = c(diss_1, diss_2, diss_3),
#     names_to = "run",
#     names_pattern = "diss_([0-9]+)",
#     values_to = "diss_score"
#   ) %>%
#   mutate(run = as.character(run)) %>%
#   select(subject, run, diss_score)

# # wide to long for nowness
# data_long_now <- data %>%
#   pivot_longer(
#     cols = c(now_own_1, now_own_2, now_own_3),
#     names_to = "run",
#     names_pattern = "now_own_([0-9]+)",
#     values_to = "now_score"
#   ) %>%
#   mutate(run = as.character(run)) %>%
#   select(subject, run, now_score)

# data_long_runs <- data_long_runs %>%
#   mutate(run = as.character(run)) %>%
#   left_join(data_long_now, by = c("subject", "run"))



#######################################################################################
--------------------------------------------------------------------------
# 1. mixed effects with random slope for run, stepwise across multiple USING THESE <<<<<<<<<<<<<<<<<<<<<<<
---------------------------------------------------------------------------
#######################################################################################

# data_long_runs_T2 <- data_long_runs_T2 %>%
#   mutate(run = as.numeric(run))
# library(performance)

# Combine left/right amygdala into one dataset with hemi variable
amygdala_data_T2 <- data_long_runs_T2 %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run)
  )

amygdala_data_T2 <- amygdala_data_T2 %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

# Run the mixed modelREPORTING THIS FOR AMYGDALA
results_amygdala_T2 <- amygdala_data_T2 %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygdala\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * CAPS.POST.PRE + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: amygdala\n")
      return(NULL)
    }
    cat("succeeded for: amygdala\n")
    return(m1)
  }) %>%
  set_names("amygdala")

map(results_amygdala_T2, summary)

# Run the mixed model                                                                   REPORTING THIS FOR AMYGDALA
results_amygdala_T2 <- list(
  hippocampus = try(
    lmer(value ~ hemi + site + age + run * CAPS.POST.PRE + (1 | subject),
         data = amygdala_data_T2,
         REML = FALSE),
    silent = TRUE
  )
)
map(results_amygdala_T2, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amygdala_T2, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amygdala_ch-scores.csv", row.names = FALSE)

# Run the mixed model                                                               NO CAPS CHECK FOR AMYGDALA
results_amygdala_T2_nocaps <- list(
  hippocampus = try(
    lmer(value ~ hemi + site + age + run + (1 | subject),
         data = amygdala_data_T2,
         REML = FALSE),
    silent = TRUE
  )
)
map(results_amygdala_T2_nocaps, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amygdala_T2_nocaps, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amygdala_ch-scores_nocaps.csv", row.names = FALSE)

# Combine left/right hippocampus into one dataset with hemi variable
hipp_data_T2 <- data_long_runs_T2 %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)
  )

hipp_data_T2 <- hipp_data_T2 %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

# Run the mixed model                                                                   REPORTING THIS FOR HIPPOCAMPUS
results_hipp_T2 <- list(
  hippocampus = try(
    lmer(value ~ hemi + site + age + run * CAPS.POST.PRE + (1 | subject),
         data = hipp_data_T2,
         REML = FALSE),
    silent = TRUE
  )
)
map(results_hipp_T2, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp_T2, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp_ch-scores.csv", row.names = FALSE)

# Run the mixed model                                                               NO CAPS CHECK FOR hippocampus
results_hipp_T2_nocaps <- list(
  hippocampus = try(
    lmer(value ~ hemi + site + age + run + (1 | subject),
         data = hipp_data_T2,
         REML = FALSE),
    silent = TRUE
  )
)
map(results_hipp_T2_nocaps, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp_T2_nocaps, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp_ch-scores_nocaps.csv", row.names = FALSE)


data_long_runs_T2 <- data_long_runs_T2 %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

results_m2_T2 <- data_long_runs_T2 %>%                               #####             REPORTING FOR VMPFC - ALSO NULL
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs_T2$ROI)) %>%
  map(function(data_long_runs_T2) {
    cat("Running ROI:", unique(data_long_runs_T2$ROI), "\n")
    m2 <- try(
      lmer(value ~ CAPS.POST.PRE * run + site + age + (1 | subject),
           data = data_long_runs_T2, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m2, "try-error")) {
      cat("failed for:", unique(data_long_runs_T2$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs_T2$ROI), "\n")
    return(m2)
  })
#results_m2
map(results_m2_T2, summary)


# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m2_T2, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_vmpfc_ch-scores.csv", row.names = FALSE)






#-----------------------------------------------------------------------------
# JUST POST TIMEPOINT TO REPLICATE THE BASELINE ANALYSES 
#-----------------------------------------------------------------------------
# does avoidance score explain variation in ROI signal?                       #####          CHANGING TO T2 ONLY
# wide to long including T1 and T2
data_long_runs_T2only <- data_T2 %>%
  pivot_longer(
    cols = c(
            #  run1vv_ch, run2vv_ch, run3vv_ch, 
            #  run1al_ch, run2al_ch, run3al_ch, 
            #  run1ar_ch, run2ar_ch, run3ar_ch, 
            #  run1hl_ch, run2hl_ch, run3hl_ch, 
            #  run1hr_ch, run2hr_ch, run3hr_ch),
             run1vv_T2, run2vv_T2, run3vv_T2,
             run1al_T2, run2al_T2, run3al_T2,
             run1ar_T2, run2ar_T2, run3ar_T2,
             run1hl_T2, run2hl_T2, run3hl_T2,
             run1hr_T2, run2hr_T2, run3hr_T2),
    names_to = c("run", "ROI"),
    names_pattern = "run([0-9]+)([a-zA-Z]+)",
    values_to = "value"
  ) %>%
  mutate(
    run = as.factor(run)
    # timepoint = ifelse(is.na(timepoint), 1, 2)
    # ROI = case_when(
    #   ROI %in% c("vv") ~ "vmpfc_v",
    #   ROI %in% c("al","ar") ~ "amygdala",
    #   ROI %in% c("hl","hr") ~ "hippocampus"
    # )
  ) %>%
  select(subject, site, caps_b_t3, age, run, ROI, value, treatment_group, CAPS.POST.PRE, STAR_ID)

data_long_runs_T2only <- data_long_runs_T2only %>%
  mutate(caps_b_t3 = as.numeric(caps_b_t3))

data_long_runs_T2only_CBT <- subset(data_long_runs_T2only, treatment_group == "2")

results_m2_T2only_group <- data_long_runs_T2only %>%                               #####           T2 ONLY REPORTING VMPFC
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs_T2only$ROI)) %>%
  map(function(data_long_runs_T2only) {
    cat("Running ROI:", unique(data_long_runs_T2only$ROI), "\n")
    m2 <- try(
      lmer(value ~ treatment_group * run + site + age + (1 | subject),
           data = data_long_runs_T2only, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m2, "try-error")) {
      cat("failed for:", unique(data_long_runs_T2only$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs_T2only$ROI), "\n")
    return(m2)
  })
#results_m2
map(results_m2_T2only_group, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m2_T2only_group, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_vmpfc_T2only_group.csv", row.names = FALSE)




# Combine left/right hippocampus into one dataset with hemi variable
hipp_data_T2only <- data_long_runs_T2only %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                                  T2 ONLY REPORTING THIS FOR HIPPOCAMPUS
results_hipp_T2only_group <- hipp_data_T2only %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * treatment_group + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp_T2only_group, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp_T2only_group, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp_T2only_group.csv", row.names = FALSE)


hipp_data_T2only <- hipp_data_T2only %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

hipp_data_T2only_CBT <- subset(hipp_data_T2only, treatment_group == "2")

# Run the mixed model                                                                  T2 ONLY REPORTING THIS FOR HIPPOCAMPUS
results_hipp_T2only_CBT <- hipp_data_T2only_CBT %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * caps_b_t3 + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp_T2only_CBT, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp_T2only_CBT, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp_T2only_CBT.csv", row.names = FALSE)




# Run the mixed model                                                               NO CAPS CHECK FOR hippocampus
results_hipp_T2only_nocaps <- list(
  hippocampus = try(
    lmer(value ~ hemi + site + age + run + (1 | subject),
         data = hipp_data_T2only,
         REML = FALSE),
    silent = TRUE
  )
)
map(results_hipp_T2only_nocaps, summary)

# # Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp_T2_nocaps, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp_T2only_nocaps.csv", row.names = FALSE)


# Combine left/right hippocampus into one dataset with hemi variable
amygdala_data_T2only <- data_long_runs_T2only %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                                  T2 ONLY REPORTING THIS FOR AMYGDALA
results_amyg_T2only_group <- amygdala_data_T2only %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygs\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * treatment_group + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: amygs\n")
      return(NULL)
    }
    cat("succeeded for: amygss\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amyg_T2only_group, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amyg_T2only_group, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amyg_T2only_group.csv", row.names = FALSE)


amygdala_data_T2only <- amygdala_data_T2only %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

amygdala_data_T2only_CBT <- subset(amygdala_data_T2only, treatment_group == "2")

# Run the mixed model                                                                  T2 ONLY REPORTING FOR AMYGDALA CBT
results_amyg_T2only_CBT <- amygdala_data_T2only_CBT %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygs\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * caps_b_t3 + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: amygs\n")
      return(NULL)
    }
    cat("succeeded for: amygss\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amyg_T2only_CBT, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amyg_T2only_CBT, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amyg_T2only_CBT.csv", row.names = FALSE)


# Run the mixed model                                                               NO CAPS CHECK FOR AMYGDALA
results_amygdala_T2only_nocaps <- list(
  hippocampus = try(
    lmer(value ~ hemi + site + age + run + (1 | subject),
         data = amygdala_data_T2only,
         REML = FALSE),
    silent = TRUE
  )
)
map(results_amygdala_T2only_nocaps, summary)


# # Extract fixed effects from all ROIs
results_table <- map_dfr(results_amygdala_T2only_nocaps, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amygdala_T2only_nocaps.csv", row.names = FALSE)





# does avoidance score explain variation in ROI signal?                       #####          CHANGING TO T1 AND T2
# wide to long including T1 and T2
data_long_runs_T1and2 <- data_T2 %>%
  pivot_longer(
    cols = c(
      run1vv, run2vv, run3vv,
      run1al, run2al, run3al,
      run1ar, run2ar, run3ar,
      run1hl, run2hl, run3hl,
      run1hr, run2hr, run3hr,
      run1vv_T2, run2vv_T2, run3vv_T2,
      run1al_T2, run2al_T2, run3al_T2,
      run1ar_T2, run2ar_T2, run3ar_T2,
      run1hl_T2, run2hl_T2, run3hl_T2,
      run1hr_T2, run2hr_T2, run3hr_T2
    ),
    names_to = c("run", "ROI", "tp"),
    names_pattern = "run([0-9]+)([a-zA-Z]+)(_T2)?",
    values_to = "value"
  ) %>%
  mutate(
    run = as.factor(run),
    timepoint = ifelse(tp == "_T2", 2L, 1L)
  ) %>%
  select(subject, site, CAPSB, caps_b_t3, CAPS.POST.PRE, treatment_group, age,
         run, ROI, timepoint, value, STAR_ID)

View(data_long_runs_T1and2)
write.csv(data_long_runs_T1and2,
          file = "data_long_runs_T1and2.csv",
          row.names = FALSE)

data_long_runs_T1and2 <- data_long_runs_T1and2 %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

results_m2_T1and2 <- data_long_runs_T1and2 %>%                               #####           T1 and 2 REPORTING VMPFC
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs_T1and2$ROI)) %>%
  map(function(data_long_runs_T1and2) {
    cat("Running ROI:", unique(data_long_runs_T1and2$ROI), "\n")
    m2 <- try(
      lmer(value ~ timepoint * treatment_group * run + site + age + (1 | subject),
           data = data_long_runs_T1and2, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m2, "try-error")) {
      cat("failed for:", unique(data_long_runs_T1and2$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs_T1and2$ROI), "\n")
    return(m2)
  })
#results_m2
map(results_m2_T1and2, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m2_T1and2, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_vmpfc-and-unilat_T1and2int.csv", row.names = FALSE)


CBTvmpfcT1and2 <- subset(data_long_runs_T1and2, treatment_group == "2")
CBTvmpfcT1and2 <- subset(CBTvmpfcT1and2, ROI == "vv")
results_m2_T1and2_CBT <- CBTvmpfcT1and2 %>%                               #####           T1 and 2 REPORTING VMPFC CBT
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(CBTvmpfcT1and2$ROI)) %>%
  map(function(CBTvmpfcT1and2) {
    cat("Running ROI:", unique(CBTvmpfcT1and2$ROI), "\n")
    results_m2_T1and2_CBT <- try(
      lmer(value ~ timepoint * run + site + age + (1 | subject),
           data = CBTvmpfcT1and2, REML = FALSE),
      silent = TRUE
    )
    if (inherits(results_m2_T1and2_CBT, "try-error")) {
      cat("failed for:", unique(CBTvmpfcT1and2$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(CBTvmpfcT1and2$ROI), "\n")
    return(results_m2_T1and2_CBT)
  })
#results_m2
map(results_m2_T1and2_CBT, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m2_T1and2_CBT, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_vmpfc-and-unilat_T1and2_CBT.csv", row.names = FALSE)


m1_vv_CBT <- results_m2_T1and2_CBT$vv
# Then run emmeans
vv_emm <- emmeans(m1_vv_CBT, ~ timepoint * run * site)
# Convert to dataframe for plotting
vv_plot_data <- as.data.frame(vv_emm)

# 4️⃣ Plot (same visual style as hippocampus/vv)
vv_plot <- ggplot(vv_plot_data, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ timepoint) +  # facet by timepoint (e.g., pre vs post)
  coord_cartesian(ylim = c(-1.25, 1)) +
  theme_minimal(base_size = 8) +
  labs(
    title = "vmPFC (CBT Group): Predicted ROI Own Cue > Other Cue Signal by Run, Site, and Timepoint",
    x = "Run",
    y = "Estimated marginal means (± SEM)",
    fill = "Site"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("vmPFC_CBT_run_site_timepoint.png", vv_plot, width = 7, height = 4, dpi = 300)



TAUvmpfcT1and2 <- subset(data_long_runs_T1and2, timepoint == "1")
TAUvmpfcT1and2 <- subset(TAUvmpfcT1and2, ROI == "vv")
#TAUvmpfcT1and2 <- subset(data_long_runs_T1and2, treatment_group == "1")

results_m2_T1and2_TAU <- TAUvmpfcT1and2 %>%                               #####           T2 ONLY REPORTING VMPFC CBT
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(TAUvmpfcT1and2$ROI)) %>%
  map(function(TAUvmpfcT1and2) {
    cat("Running ROI:", unique(TAUvmpfcT1and2$ROI), "\n")
    results_m2_T1and2_TAU <- try(
      lmer(value ~ timepoint * run + site + age + (1 | subject),
           data = TAUvmpfcT1and2, REML = FALSE),
      silent = TRUE
    )
    if (inherits(results_m2_T1and2_TAU, "try-error")) {
      cat("failed for:", unique(TAUvmpfcT1and2$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(TAUvmpfcT1and2$ROI), "\n")
    return(results_m2_T1and2_TAU)
  })
#results_m2
map(results_m2_T1and2_TAU, summary)

m1_vv_TAU <- results_m2_T1and2_TAU$vv
# Then run emmeans
vv_emm_TAU <- emmeans(m1_vv_TAU, ~ timepoint * run * site)
# Convert to dataframe for plotting
vv_plot_data_TAU <- as.data.frame(vv_emm_TAU)

# 4️⃣ Plot (same visual style as hippocampus/vv)
vv_plot <- ggplot(vv_plot_data_TAU, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ timepoint) +  # facet by timepoint (e.g., pre vs post)
  coord_cartesian(ylim = c(-0.75, 0.5)) +
  theme_minimal(base_size = 8) +
  labs(
    title = "vmPFC (TAU Group): Predicted ROI Own Cue > Other Cue Signal by Run, Site, and Timepoint",
    x = "Run",
    y = "Estimated marginal means (± SEM)",
    fill = "Site"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("vmPFC_TAU_run_site_timepoint.png", vv_plot, width = 7, height = 4, dpi = 300)



vmpfcT1_T2script <- subset(data_long_runs_T1and2, timepoint == "1")
vmpfcT1_T2script <- subset(vmpfcT1_T2script, ROI == "vv")

results_m2_T1 <- vmpfcT1_T2script %>%                               #####           T1 ONLY REplicating
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(vmpfcT1_T2script$ROI)) %>%
  map(function(vmpfcT1_T2script) {
    cat("Running ROI:", unique(vmpfcT1_T2script$ROI), "\n")
    results_m2_T1 <- try(
      lmer(value ~ run + site + age + (1 | subject),
           data = vmpfcT1_T2script, REML = FALSE),
      silent = TRUE
    )
    if (inherits(results_m2_T1, "try-error")) {
      cat("failed for:", unique(vmpfcT1_T2script$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(vmpfcT1_T2script$ROI), "\n")
    return(results_m2_T1)
  })
#results_m2
map(results_m2_T1, summary)

m1_vv_T1 <- results_m2_T1$vv                                #                  REPLICATING m1_vv IN SCRIPT 1
# Then run emmeans
vv_emm_T1 <- emmeans(m1_vv_T1, ~ run * site)
# Convert to dataframe for plotting
vv_plot_data_T1 <- as.data.frame(vv_emm_T1)

# 4️⃣ Plot (same visual style as hippocampus/vv)
vv_plot <- ggplot(vv_plot_data_T1, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  #facet_wrap(~ timepoint) +  # facet by timepoint (e.g., pre vs post)
  coord_cartesian(ylim = c(-1.0, 1.5)) +
  theme_minimal(base_size = 8) +
  labs(
    title = "vmPFC (T1): Predicted ROI Own Cue > Other Cue Signal by Run, Site",
    x = "Run",
    y = "Estimated marginal means (± SEM)",
    fill = "Site"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("vmPFC_T1_run_site.png", vv_plot, width = 7, height = 4, dpi = 300)


CBTvmpfcT2 <- subset(CBTvmpfcT1and2, timepoint == "2")

results_m2_T2_CBT <- CBTvmpfcT2 %>%                               #####           T2 ONLY REPORTING VMPFC CBT
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(CBTvmpfcT2$ROI)) %>%
  map(function(CBTvmpfcT2) {
    cat("Running ROI:", unique(CBTvmpfcT2$ROI), "\n")
    m2 <- try(
      lmer(value ~ run * CAPS.POST.PRE + site + age + (1 | subject),
           data = CBTvmpfcT2, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m2, "try-error")) {
      cat("failed for:", unique(CBTvmpfcT2$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(CBTvmpfcT2$ROI), "\n")
    return(m2)
  })
#results_m2
map(results_m2_T2_CBT, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m2_T2_CBT, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_vmpfc-and-unilat_T2_CBT.csv", row.names = FALSE)



# Combine left/right into one dataset with hemi variable
amygdala_data_T1and2 <- data_long_runs_T1and2 %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run)
  )
# Run the mixed model                                                                  T1 and T2 REPORTING THIS FOR AMYGDALA
results_amyg_T1and2 <- amygdala_data_T1and2 %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygs\n")
    m1 <- try(
      lmer(value ~ hemi + timepoint * treatment_group + run + site + age + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: amygs\n")
      return(NULL)
    }
    cat("succeeded for: amygss\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amyg_T1and2, summary)

# Run the mixed model                                                                  T1 and T2 REPORTING THIS FOR AMYGDALA
results_amyg_T1and2 <- amygdala_data_T1and2 %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygs\n")
    m1 <- try(
      lmer(value ~ hemi + timepoint * treatment_group * run + site + age + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: amygs\n")
      return(NULL)
    }
    cat("succeeded for: amygss\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amyg_T1and2, summary)

CBTamygT1and2 <- subset(amygdala_data_T1and2, treatment_group == "2")

results_amyg_T1and2_CBT <- CBTamygT1and2 %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygs\n")
    m1 <- try(
      lmer(value ~ hemi + timepoint * run + site + age + (1 | subject),
      #lmer(value ~ hemi + timepoint + CAPS.POST.PRE * run + site + age + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: amygs\n")
      return(NULL)
    }
    cat("succeeded for: amygss\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amyg_T1and2_CBT, summary)


m1_amyg_CBT <- results_amyg_T1and2_CBT$amygdala
# Then run emmeans
amyg_emm <- emmeans(m1_amyg_CBT, ~ timepoint * run * site)
# Convert to dataframe for plotting
amyg_plot_data <- as.data.frame(amyg_emm)

# 4️⃣ Plot (same visual style as hippocampus/vv)
amyg_plot <- ggplot(amyg_plot_data, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ timepoint) +  # facet by timepoint (e.g., pre vs post)
  coord_cartesian(ylim = c(-0.75, 1.5)) +
  theme_minimal(base_size = 8) +
  labs(
    title = "Amygdala (CBT Group): Predicted ROI Own Cue > Other Cue Signal by Run, Site, and Timepoint",
    x = "Run",
    y = "Estimated marginal means (± SEM)",
    fill = "Site"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# 5️⃣ Save the plot
ggsave("amygdala_CBT_run_site_timepoint.png", amyg_plot, width = 7, height = 4, dpi = 300)



# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amyg_T1and2_CBT, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amyg_T1and2_CBT.csv", row.names = FALSE)

# Combine left/right hippocampus into one dataset with hemi variable
hipp_data_T1and2 <- data_long_runs_T1and2 %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                                  T1 and T2 REPORTING THIS FOR HIPPOCAMPUS
results_hipp_T1and2 <- hipp_data_T1and2 %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(value ~ hemi + timepoint * treatment_group + run + site + age + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp_T1and2, summary)


CBThippT1and2 <- subset(hipp_data_T1and2, treatment_group == "2")

results_hipp_T1and2_CBT <- CBThippT1and2 %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(value ~ hemi + timepoint * run + site + age + (1 | subject),
      #lmer(value ~ hemi + timepoint * treatment_group * run + site + age + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp_T1and2_CBT, summary)


m1_hipp_CBT <- results_hipp_T1and2_CBT$hippocampus
# Then run emmeans
hipp_emm <- emmeans(m1_hipp_CBT, ~ timepoint * run * site)
# Convert to dataframe for plotting
hipp_plot_data <- as.data.frame(hipp_emm)

hipp_plot <- ggplot(hipp_plot_data, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ timepoint) +  # facet by timepoint (e.g., pre vs post)
  coord_cartesian(ylim = c(-0.25, 1.25)) +
  theme_minimal(base_size = 8) +
  labs(
    title = "Hippocampus (CBT Group): Predicted ROI Own Cue > Other Cue Signal by Run, Site, and Timepoint",
    x = "Run",
    y = "Estimated marginal means (± SEM)",
    fill = "Site"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("hippocampus_CBT_run_site_timepoint.png", hipp_plot, width = 7, height = 4, dpi = 300)



# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp_T1and2_CBT, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp_T1and2_CBT.csv", row.names = FALSE)



#-----------------------------------------------------------------------------
# USING BETAS TO DO STATISTICS ON THE WHOLE BRAIN ACTIVITY AS WELL AS ROI
#-----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lme4)


cols_to_gather <- c(
  # Hippocampus left/right T1
  "run1hl_beta_own","run2hl_beta_own","run3hl_beta_own",
  "run1hl_beta_other","run2hl_beta_other","run3hl_beta_other",
  "run1hr_beta_own","run2hr_beta_own","run3hr_beta_own",
  "run1hr_beta_other","run2hr_beta_other","run3hr_beta_other",
  # Wholebrain T1
  "run1whbr_beta_own","run2whbr_beta_own","run3whbr_beta_own",
  # Hippocampus and wholebrain T2
  "run1hl_beta_own_T2","run2hl_beta_own_T2","run3hl_beta_own_T2",
  "run1hl_beta_other_T2","run2hl_beta_other_T2","run3hl_beta_other_T2",
  "run1hr_beta_own_T2","run2hr_beta_own_T2","run3hr_beta_own_T2",
  "run1hr_beta_other_T2","run2hr_beta_other_T2","run3hr_beta_other_T2",
  "run1whbr_beta_own_T2","run2whbr_beta_own_T2","run3whbr_beta_own_T2"
)


data_long_T1T2_betas <- data_T2 %>%
  pivot_longer(
    cols = all_of(cols_to_gather),
    names_to = c("run","ROI","condition","tp"),
    names_pattern = "run([123])([a-z]+)_beta_(own|other)(_T2)?",
    values_to = "value"
  ) %>%
  mutate(
    run = factor(run, levels = c("1","2","3")),
    condition = factor(condition, levels = c("own","other")),
    timepoint = ifelse(tp=="_T2", 2L, 1L),
    hemi = case_when(
      ROI %in% c("hl","hr") ~ ifelse(ROI=="hl","L","R"),
      TRUE ~ NA_character_
    ),
    ROI = case_when(
      ROI %in% c("hl","hr") ~ "hippocampus",
      ROI %in% c("whbr") ~ "wholebrain",
      TRUE ~ ROI
    )
  ) %>%
  select(subject, site, CAPS.POST.PRE, caps_b_t3, age,
         run, ROI, hemi, condition, timepoint, treatment_group, value, STAR_ID)


hippo_long <- data_long_T1T2_betas %>%
  filter(ROI == "hippocampus", hemi %in% c("L","R"))
hippo_long <- hippo_long %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))


library(tidyverse)

# filter hippocampus data set factors
hipp_data <- hippo_long %>%
  filter(ROI == "hippocampus") %>%
  mutate(
    run = as.factor(run),
    condition = as.factor(condition),
    hemi = as.factor(hemi),
    subject = as.factor(subject)
  )

hipp_data <- subset(hipp_data, timepoint == "2")

# compute mean ± SEM per run × condition × hemi
hipp_summary <- hipp_data %>%
  group_by(run, condition, hemi) %>%
  summarise(
    mean_signal = mean(value, na.rm = TRUE),
    sem = sd(value, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# plot grp means
hipp_plot <- ggplot(hipp_summary, aes(x = run, y = mean_signal, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_signal - sem, ymax = mean_signal + sem),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ hemi) +
  coord_cartesian(ylim = c(-0.5, 0.75)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Hippocampus: Mean Beta Values by Run and Condition",
    x = "Run",
    y = "Mean estimated activation (± SEM)",
    fill = "Condition"
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# 4️⃣ Save the plot
ggsave("T2_hippocampus_run_condition_raw.png", hipp_plot, width = 6, height = 4, dpi = 300)


hippo_long_CBT <- subset(hippo_long, treatment_group == "2")

hippo_long_TAU <- subset(hippo_long, treatment_group == "1")

hipp_betas<- lmer(
  value ~ hemi + condition + CAPS.POST.PRE + run * timepoint + treatment_group +age + site + (1 | subject),
  data = hippo_long
)
summary(hipp_betas)

results_table <- broom.mixed::tidy(hipp_betas) %>%
  select(term, estimate, std.error, df, statistic, p.value)

write.csv(results_table, "results_hipp_beta_T1and2.csv", row.names = FALSE)


library(dplyr)

# Create a separate dataset for wholebrain
wb_long <- data_long_T1T2_betas %>%
  filter(ROI == "wholebrain") %>%
  group_by(subject, site, CAPS.POST.PRE, caps_b_t3, age, condition, timepoint, treatment_group, STAR_ID) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ROI = "wholebrain",
    hemi = NA  # or "WB" if you prefer
  ) 
  

wb_long <- wb_long %>%
  mutate(CAPS.POST.PRE = as.numeric(CAPS.POST.PRE))

model_betas_whbr <- lmer(
  value ~ timepoint * CAPS.POST.PRE + site + age + (1 | subject),
  data = wb_long
)
summary(model_betas_whbr)


results_table <- broom.mixed::tidy(model_betas_whbr) %>%
  select(term, estimate, std.error, df, statistic, p.value)

write.csv(results_table, "2_results_whbr_beta_T1and2.csv", row.names = FALSE)


wb_long_TAU <- subset(wb_long, treatment_group == "1")

model_betas_whbr_CBT <- lmer(
  value ~ timepoint * CAPS.POST.PRE + site + age + (1 | subject),
  data = wb_long_CBT
)
summary(model_betas_whbr_CBT)


results_table <- broom.mixed::tidy(model_betas_whbr_TAU) %>%
  select(term, estimate, std.error, df, statistic, p.value)

write.csv(results_table, "results_whbr_beta_T1and2_TAU.csv", row.names = FALSE)


# Extract fixed effects from all ROIs
results_table <- map_dfr(models_betas_whbr, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_whbr_T1and2_CBT.csv", row.names = FALSE)






































#------------------------------------------------------------------------------           T1 REPORTING TRIAL BY TRIAL NOWNESS
# nowness by trial for vmpfc

trial_data <- read.csv("/Users/zefreeman/Documents/T1_nowness_ratings.csv", header = TRUE, stringsAsFactors = FALSE)

trial_data <- trial_data %>%
  mutate(Run = as.integer(Run)) %>%
  rename(STAR_ID = SubjectID)

vv_data <- data_long_runs %>%
  filter(ROI == "vv") %>%                    # keep only vv ROI
  select(subject, site, CAPSB, age, run, value, STAR_ID) %>%  # select desired columns
  rename(
    vv_value = value,                        # rename value to vv_value
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer

trial_data_vv <- trial_data %>%
  left_join(vv_data, by = c("STAR_ID" = "STAR_ID", "Run" = "Run"))

trial_data_vv <- trial_data_vv %>%
  mutate(
    CAPSB = as.numeric(CAPSB),
    run = as.factor(Run),
    site = as.factor(site)
  )
model_trial_ratingbyrun <-lmer(vv_value ~ Rating * run + CAPSB + Trial + site + age + (1 | subject),
                             data = trial_data_vv)
summary(model_trial_ratingbyrun)
model_trial <- lmer(vv_value ~ CAPSB + Rating + run + Trial + site + age + (1 | subject),
                    data = trial_data_vv)
summary(model_trial)

results_table <- tidy(model_trial_ratingbyrun) %>%
  select(term, estimate, std.error, df, statistic, p.value)

# Write to CSV
write.csv(results_table, "results_vmpfc_runbynowness.csv", row.names = FALSE)

model_data <- model.frame(model_trial)
model_data$Run <- trial_data_vv$Run[as.numeric(rownames(model_data))]

model_data$predicted <- predict(model_trial, re.form = NA)

# plot
library(ggplot2)
plot_vv <- ggplot(model_data, aes(x = Trial, y = vv_value, color = factor(Run))) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = predicted, group = Run), size = 1) +
  labs(
    title = "Trial-by-trial vv_value with Run-specific predictions",
    x = "Trial",
    y = "vv_value",
    color = "Run"
  ) +
  theme_minimal()

print(plot_vv)
ggsave("vv_trial_by_trial.png", plot_vv, width = 8, height = 6, dpi = 300)


# Ensure uniqueness per subject x Run for AR and AL
ar_data <- ar_data %>%
  distinct(STAR_ID, Run, .keep_all = TRUE)

al_data <- al_data %>%
  distinct(STAR_ID, Run, .keep_all = TRUE)

# Merge AR and AL together (many-to-many allowed)
arl_data <- ar_data %>%
  left_join(al_data, by = c("STAR_ID", "Run", "subject", "site", "CAPSB", "age"), relationship = "many-to-many")

arl_data <- arl_data %>%
  mutate(Run = as.integer(Run))

trial_data_arl <- trial_data %>%
  left_join(arl_data, by = c("STAR_ID" = "STAR_ID", "Run" = "Run"))


trial_data_arl_long <- trial_data_arl %>%
  pivot_longer(
    cols = c(al_value, ar_value),
    names_to = "Hemisphere",
    values_to = "amyg_value"
  ) %>%
  mutate(
    Hemisphere = factor(Hemisphere, levels = c("al_value", "ar_value")),  # left as reference
    CAPSB = as.numeric(CAPSB),
    run = as.factor(Run),
    site = as.factor(site)
  )
model_arl_ratingbyrun <- lmer(
amyg_value ~ Rating * run + CAPSB + Trial + site + age + Hemisphere + (1 | subject),
  data = trial_data_arl_long
)
summary(model_arl_ratingbyrun)

results_table <- tidy(model_arl_ratingbyrun) %>%
  select(term, estimate, std.error, df, statistic, p.value)

# Write to CSV
write.csv(results_table, "results_amygdala_runbynowness.csv", row.names = FALSE)

model_arl <- lmer(
  amyg_value ~ CAPSB + Rating + run + Trial + site + age + Hemisphere + (1 | subject),
  data = trial_data_arl_long
)
summary(model_arl)

# 1️⃣ Convert to long format for plotting
plot_data_long <- trial_data_arl %>%
  pivot_longer(
    cols = c(al_value, ar_value),
    names_to = "Hemisphere",
    values_to = "value"
  ) %>%
  mutate(Hemisphere = factor(Hemisphere, levels = c("al_value", "ar_value")))

# 2️⃣ Compute model predictions for fixed effects only
# Use only the rows actually used in the model
model_rows <- model.frame(model)  # or model_trial if that's your model
preds <- predict(model, newdata = model_rows, re.form = NA)

# 3️⃣ Add predictions to the long-format dataframe
# Make sure number of rows matches
plot_data_long <- plot_data_long[1:length(preds), ]
plot_data_long$predicted <- preds

# 4️⃣ Plot
plot_arl<-ggplot(plot_data_long, aes(x = Trial, y = value, color = Hemisphere)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(y = predicted), size = 1) +
  facet_wrap(~Run) +
  labs(
    title = "Amygdala Values (AR vs AL) Across Trials and Runs",
    x = "Trial",
    y = "Amygdala Signal",
    color = "Hemisphere"
  ) +
  theme_minimal()
  ggsave("arl_trial_by_trial.png", plot_arl, width = 8, height = 6, dpi = 300)

#------------------------------------------------------------------------------
# AVOID SCORE INCLUDED, INTERACTION W/CAPS BY RUN
#      lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + (1 | subject),

results_m3 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m3 <- try(
      lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m3, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m3)
  })
#results_m2
map(results_m3, summary)

vif_results3 <- results_m3 %>%
  compact() %>%  # remove NULLs if any
  map(get_vif_lmer)
vif_results3

#------------------------------------------------------------------------------
# AVOID SCORE INCLUDED, INTERACTION W/CAPS BY RUN                                     NOT REPORTING EXPLORATORY 1
#      lmer(roi_signal ~  avoid_score * run + CAPSB + site + age + (1 | subject),

results_m3.5 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m3.5 <- try(
      lmer(roi_signal ~  avoid_score * run + CAPSB + site + age + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m3.5, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m3.5)
  })
#results_m2
map(results_m3.5, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m3.5, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_m3.5.csv", row.names = FALSE)

results_hipp <- hipp_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(roi_signal ~ avoid_score * run + CAPSB + site + age + hemi + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp, summary)
# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp.csv", row.names = FALSE)


results_amygdala <- amygdala_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: AMYG\n")
    m1 <- try(
      lmer(roi_signal ~ avoid_score * run + CAPSB + site + age + hemi + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: amyg\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amygdala, summary)
# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amygdala, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amyg.csv", row.names = FALSE)




# --------------------------------------------------------------------------------
# Save summary table of results to .csv
# --------------------------------------------------------------------------------
          library(readr)

          # Combine all fixed effects into one table
          results_table <- map_dfr(names(results_m3.5), function(roi_name) {
            model <- results_m3.5[[roi_name]]
            if (!is.null(model)) {
              broom.mixed::tidy(model, effects = "fixed") %>%
                select(term, estimate, std.error, df, statistic, p.value) %>%
                mutate(
                  across(where(is.numeric), ~round(., 4)),
                  ROI = roi_name
                )
            } else {
              tibble(
                ROI = roi_name,
                term = NA,
                estimate = NA,
                std.error = NA,
                df = NA,
                statistic = NA,
                p.value = NA
              )
            }
          })

          # Add model fit indices (AIC, BIC, logLik, nobs) per ROI
          fit_table <- map_dfr(names(results_m3.5), function(roi_name) {
            model <- results_m3.5[[roi_name]]
            if (!is.null(model)) {
              fit <- broom.mixed::glance(model)
              tibble(
                ROI = roi_name,
                AIC = round(fit$AIC, 3),
                BIC = round(fit$BIC, 3),
                logLik = round(fit$logLik, 3),
                nobs = fit$nobs
              )
            } else {
              tibble(ROI = roi_name, AIC = NA, BIC = NA, logLik = NA, nobs = NA)
            }
          })

          # Merge both tables (so fit info repeats for each term)
          final_table <- left_join(results_table, fit_table, by = "ROI") %>%
            relocate(ROI, .before = term)

          # Save to .tsv
          write_csv(final_table, "results_m3.5_summary.csv")


vif_results3.5 <- results_m3.5 %>%
  compact() %>%  # remove NULLs if any
  map(get_vif_lmer)
vif_results3.5
      

#------------------------------------------------------------------------------
# PENULTIMATE, same as above INCLUDING RANDOM SLOPE FOR RUN
# lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + (1 + run | subject),              REPORT PROPERLY AS SECOND MODEL

results_m4 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m4 <- try(
    #   lmer(roi_signal ~ avoid_score + CAPSB * run + site + age + (1 | subject/run),
    #  data = data_long_runs, REML = FALSE),
      lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m4, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m4)
  })
#results_m2
map(results_m4, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m4, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_m4.csv", row.names = FALSE)

results_hipp <- hipp_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + hemi +(1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp, summary)
# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp.csv", row.names = FALSE)


results_amygdala <- amygdala_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: AMYG\n")
    m1 <- try(
      lmer(roi_signal ~ avoid_score + CAPSB * run + site + age + hemi +(1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: amyg\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amygdala, summary)
# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amygdala, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amyg.csv", row.names = FALSE)


      
#-------------------------------------------------------------------------------
# LAST MODEL WITH THE MOST INTERACTIONS AND RUN EFFECTS
# lmer(roi_signal ~ CAPSB * avoid_score * run + site + age + (1 + run | subject),


results_m5 <- data_long_runs %>%
  mutate(run = as.numeric(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m5 <- try(
      lmer(roi_signal ~ CAPSB * avoid_score * run + site + age + (1 + run | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m5, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m5)
  })
results_m5
map(results_m5, summary)


      #lmer(roi_signal ~ CAPSB * avoid_score * run + site + age + (1 + subject),

      #lmer(roi_signal ~ CAPSB * avoid_score * run + site + age + (1 + run | subject),

    #  lmer(roi_signal ~  NOWNESS_SCORE + CAPSB * run + site + age + (1 + run | subject), #


#------------------------------------------------------------------------------
# NOWNESS SCORE INCLUDED, INTERACTION W/NOWNESS BY RUN PLUS CAPS
#      lmer(roi_signal ~  now_own_mean * run + CAPSB + site + age + (1 | subject),
# Convert to numeric (before running the model)
data_long_runs$now_own_mean <- as.numeric(as.character(data_long_runs$now_own_mean))

results_m6 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m6 <- try(
      lmer(roi_signal ~  now_own_mean * run + CAPSB + site + age + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m6, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m6)
  })
#results_m2
map(results_m6, summary)
# results nicer
walk2(results_m6, names(results_m3.5), function(model, roi_name) {
  if (!is.null(model)) {
    cat("\n\n============================\n")
    cat("ROI:", roi_name, "\n")
    cat("============================\n")
    broom.mixed::tidy(model, effects = "fixed") %>%
      select(term, estimate, std.error, df, statistic, p.value) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      print(n = Inf)
  } else {
    cat("\nROI:", roi_name, "— model failed.\n")
  }
})


#------------------------------------------------------------------------------
# NOWNESS SCORE INCLUDED, INTERACTION W/CAPS BY RUN
#      lmer(roi_signal ~  now_score * run + CAPSB + site + age + (1 | subject),

results_m7 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m7 <- try(
      lmer(roi_signal ~  now_score * run + CAPSB + site + age + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m7, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m7)
  })
#results_m2
map(results_m7, summary)
# results nicer
walk2(results_m7, names(results_m7), function(model, roi_name) {
  if (!is.null(model)) {
    cat("\n\n============================\n")
    cat("ROI:", roi_name, "\n")
    cat("============================\n")
    broom.mixed::tidy(model, effects = "fixed") %>%
      select(term, estimate, std.error, df, statistic, p.value) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      print(n = Inf)
  } else {
    cat("\nROI:", roi_name, "— model failed.\n")
  }
})



###########################
# plots 20251008
########################### 

# Extract subjects used in the right hippocampus model
subjects_in_model <- unique(m_hr@frame$subject)
# filter to those only
data_hr_fit <- data_hr %>% filter(subject %in% subjects_in_model)

# Add fitted values from the model
data_hr_fit <- data_hr_fit %>%
  mutate(
    fitted = predict(m_hr, newdata = data_hr_fit),  # create fitted column first
    run = as.factor(run)                             # convert run to factor
  )



# Plot fitted vs. avoidance
hr_fittedvsavoid <- ggplot(data_hr_fit, aes(x = avoid_score, y = fitted, color = run)) +
  geom_point(alpha = 0.6) +                         # scatter points
  geom_smooth(method = "lm", se = TRUE) +          # linear trend lines per run
  labs(
    x = "Avoidance Score",
    y = "Fitted ROI Signal (hr)",
    color = "Run",
    title = "Fitted ROI Signal vs. Avoidance Score by Run"
  ) +
  theme_minimal()
ggsave(
  filename = "hr_fitted_vs_avoidance.png",  # file name
  plot = hr_fittedvsavoid,
  width = 7,      # width in inches
  height = 5,     # height in inches
  dpi = 300       # resolution
)


-----------------------------------------------------------------------------------
### 2. models per ROI: CAPSB ~ ROI + avg_avoid + ROI:avg_avoid + site (covariate)
# ---------------------------------------------------------------------------------

# # linear models per ROI total score with avoidance averaged !
# fit_roi_lm <- function(roi_col, df) {
#   formula_text <- paste0("CAPSB ~ ", roi_col, " * avg_avoid + site")
#   mod <- lm(as.formula(formula_text), data = df)
#   tidy_mod <- broom::tidy(mod)
#   glance_mod <- broom::glance(mod)
#   list(model = mod, tidy = tidy_mod, glance = glance_mod)
# }

fit_roi_lmer_df <- function(data_long_runs) {
  # Mixed model per ROI
mod <- lmer(
  roi_signal ~ CAPSB * avoid_score * factor(run) + site + age + (1 | subject),
  data = data_long_runs
)
  # tidy fixed effects
  broom.mixed::tidy(mod, effects = "fixed") %>%
    mutate(
      r.squared = NA,      # optional: compute later if desired
      adj.r.squared = NA
    )
}
results <- data_long_runs %>%
  group_by(ROI) %>%
  group_modify(~ fit_roi_lmer_df(.x))
results


##------------------------------------------------------------------------------------
## not needed
##-----------------------------------------------------------------------------------

# compute average avoidance across runs
data <- data %>%
  mutate(avg_avoid = rowMeans(select(., avo_1, avo_2, avo_3), na.rm = TRUE))
# list of ROI column names (replace with your actual ROI variable names)
roi_cols <- c("vmpfc_v_total", "amyg_l_total", "amyg_r_total", 
              "hipp_l_total", "hipp_r_total")
# apply model-fitting function to each ROI
roi_results <- map(roi_cols, ~ fit_roi_lm(.x, data)) %>%
  set_names(roi_cols)
# extract  tidy summaries and add ROI name
roi_tidy <- map_dfr(roi_results, "tidy", .id = "ROI")
# model-level summaries (R², AIC, etc.)
roi_glance <- map_dfr(roi_results, "glance", .id = "ROI")
print(roi_tidy)
print(roi_glance)



--------------------------------------------------------------------------------------
### 3. models per ROI and run: CAPSB ~ ROI + avoid_score + ROI:avoid_score + site (covariate)
# ------------------------------------------------------------------------------------

# Ensure run is treated as a factor (categorical)
data_long_runs <- data_long_runs %>%
  mutate(run = factor(run, levels = c("1", "2", "3")))

# Function to fit an LME model per ROI
fit_roi_lme_by_run <- function(roi_name, df) {
  model <- lmer(avoid_score ~ roi_signal * run + site + CAPSB+ (1 | subject), 
                data = df %>% filter(ROI == roi_name))
  tidy_model <- tidy(model)
  tidy_model$ROI <- roi_name
  return(tidy_model)
}

# Apply the model to each ROI
roi_models_by_run <- data_long_runs %>%
  split(.$ROI) %>%
  map_dfr(~ fit_roi_lme_by_run(unique(.x$ROI), .x))

# Keep only fixed effects
roi_results_by_run <- roi_models_by_run %>%
  filter(effect == "fixed")

# View the key terms
roi_results_by_run %>%
  select(ROI, term, estimate, std.error, statistic, p.value)




# ---------------------------------------------------------------------------------
### 4. plots
# ---------------------------------------------------------------------------------

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)

# Make folder for plots
if (!dir.exists("ROI_scatterplots")) dir.create("ROI_scatterplots")

# Get list of unique ROIs
unique_rois <- unique(data_long_runs$ROI)

for (roi_name in unique_rois) {
  df_roi <- data_long_runs %>% filter(ROI == roi_name)
  
  # Fit LME model: roi_signal ~ CAPSB + site + (1|subject)
  lme_model <- lmer(roi_signal ~ CAPSB + site + (1 | subject), data = df_roi)
  
  # Create predicted values for plotting
  df_roi <- df_roi %>%
    mutate(predicted = predict(lme_model, newdata = df_roi, re.form = NA)) # fixed effects only
  
  # Plot with points and LME-predicted lines per run
  p <- ggplot(df_roi, aes(x = CAPSB, y = roi_signal, color = factor(run))) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(aes(y = predicted, group = interaction(run)), size = 1) +
    labs(
      title = paste0("ROI: ", roi_name),
      x = "CAPSB Score",
      y = "ROI Activation",
      color = "Run"
    ) +
    theme_minimal(base_size = 14)
  
  # Print interactively
  print(p)
  
  # Save to file
  ggsave(
    filename = paste0("ROI_scatterplots/ROI_scatter_", roi_name, ".png"),
    plot = p,
    width = 6,
    height = 4
  )
}




# Extract the interaction terms and p-values for ROI:avg_avoid
interaction_results <- map_dfr(names(results_lm), function(rc) {
  tid <- results_lm[[rc]]$tidy
  inter_row <- tid %>% filter(term == paste0(rc, ":avg_avoid") | term == paste0("avg_avoid:", rc))
  if(nrow(inter_row) == 0) {
    # sometimes lm term order differs; try pattern match
    inter_row <- tid %>% filter(grepl("avg_avoid.*", term) & grepl(rc, term))
  }
  if(nrow(inter_row) == 0) {
    return(tibble(ROI = rc, estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_))
  } else {
    inter_row %>% transmute(ROI = rc, estimate, std.error, statistic, p.value)
  }
})


--------------------------------------------------------------------------------
### model per ROI: CAPSB ~ ROI + avg_avoid + ROI:avg_avoid + site_no (covariate)
# ------------------------------------------------------------------------------

### total ROIs t-testss
data_totals_filtered <- data_long %>% filter(run == "total")
ROIs_totals <- unique(data_totals_filtered$ROI)
# t-tests
for(r in ROIs_totals) {
  cat("T-test for ROI:", r, "\n")
  roi_values <- data_totals_filtered %>% 
                filter(ROI == r) %>% 
                pull(value)
  print(t.test(roi_values, mu = 0))
}

### total ROIs vs CAPSB (now including average avoidance scores)
model_totals <- lmer(value ~ ROI * CAPSB + avg_avoid + site_no + (1 | subject), 
                     data = data_long %>% filter(run=="total"))
summary(model_totals)

# run (1,2,3) specific ROIs vs CAPSB (now including avoidance scores)<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
model <- lmer(value ~ run + CAPSB + avoid_score + run:CAPSB + run:avoid_score + (1 | subject), 
              data = data_long %>% filter(run != "total"))
summary(model)
# run3 shows the strongest avoidance-related effects 
# might indicate that avoidance responses become more pronounced...

## which ROI is driving this effect???
# run (1,2,3) specific ROIs vs CAPSB, separate models per ROI (now including avoidance scores)
ROIs <- unique(data_long_runs$ROI)
models_list <- list()
for(r in ROIs) {
  cat("Processing ROI:", r, "\n")
  mod <- lmer(value ~ run * CAPSB + avoid_score + run:avoid_score + (1 | subject), #add int. of CAPSB:avoid
              data = data_long_runs %>% filter(ROI == r))
  print(summary(mod))
  models_list[[r]] <- mod
}


####################### results:
#### vmpfc:
# main effect of avoidance: β = -0.22, p = 0.039* - higher avoidance scores predict lower VMPFCv activation in run1
# run2 × Avoidance interaction: β = 0.37, p = 0.008** - avoidance effect becomes more positive in run2
# run3 × Avoidance interaction: β = 0.49, p = 0.0006*** - strongest positive avoidance effect in run3
# run1: Higher avoidance = Lower activation
# run2: Avoidance effect becomes neutral/slightly positive
# run3: Higher avoidance = sig higher activation
# noise?

# amygdala left (al): sig run effects (runs 2&3 < run1) but no avoid effects
# amygdala right (ar): no CAPS no avoid
# left hipp: no avoid
# r hipp: sig run3 × CAPSB interaction (p = 0.030*), but no avoid effects
#########################


# additional analysis: models with all three bh scores
# run (1,2,3) specific ROIs vs CAPSB with all bh measures
model_full <- lmer(value ~ run + CAPSB + avoid_score + anxiety_score + diss_score + 
                   run:CAPSB + run:avoid_score + (1 | subject), 
                   data = data_long %>% filter(run != "total"))
summary(model_full)

# separate models per ROI with all bh measures
models_list_full <- list()
for(r in ROIs) {
  cat("Processing ROI with full model:", r, "\n")
  mod_full <- lmer(value ~ run * CAPSB + avoid_score + anxiety_score + diss_score + 
                   run:avoid_score + (1 | subject),
                   data = data_long_runs %>% filter(ROI == r))
  print(summary(mod_full))
  models_list_full[[r]] <- mod_full
}


# # wide to long for run-specific ROIs
# data_long_runs <- data %>%
#   pivot_longer(
#     cols = c(run1vv, run2vv, run3vv, 
#              run1al, run2al, run3al, 
#              run1ar, run2ar, run3ar, 
#              run1hl, run2hl, run3hl, 
#              run1hr, run2hr, run3hr),
#     names_to = c("run", "ROI"),
#     names_pattern = "run([0-9]+)([a-zA-Z]+)",
#     values_to = "value"
#   ) %>%
#   select(subject, site_no, CAPSB, run, ROI, value)   # keep covariates

# # wide to long for total ROIs
# data_totals <- data %>%
#   pivot_longer(
#     cols = c(vmpfc_v_total,
#              amyg_l_total, amyg_r_total, 
#              hipp_l_total, hipp_r_total),
#     names_to = "ROI",
#     names_pattern = "(.+)_total",
#     values_to = "value"
#   ) %>%
#   mutate(run = "total") %>%
#   select(subject, site_no, CAPSB, run, ROI, value)


# # combine both long
# data_long <- bind_rows(data_long_runs, data_totals)
# print(data_long)

# ### total ROIs t-tests
# data_totals_filtered <- data_long %>% filter(run == "total")
# ROIs_totals <- unique(data_totals_filtered$ROI)
# # t-tests
# for(r in ROIs_totals) {
#   cat("T-test for ROI:", r, "\n")
#   roi_values <- data_totals_filtered %>% 
#                 filter(ROI == r) %>% 
#                 pull(value)
#   print(t.test(roi_values, mu = 0))
# }

# ### total ROIs vs CAPSB
# model_totals <- lmer(value ~ ROI * CAPSB + site_no (1 | subject), data = data_long %>% filter(run=="total"))
# summary(model_totals)

# # run (1,2,3) specific ROIs vs CAPSB
# model <- lmer(value ~ run + CAPSB + run:CAPSB + (1 | subject), data = data_long %>% filter(run != "total"))
# summary(model)

# # run (1,2,3) specific ROIs vs CAPSB, separate models per ROI
# ROIs <- unique(data_long_runs$ROI)
# models_list <- list()
# for(r in ROIs) {
#   cat("Processing ROI:", r, "\n")
#   mod <- lmer(value ~ run * CAPSB + (1 | subject),
#               data = data_long_runs %>% filter(ROI == r))
#   print(summary(mod))
#   models_list[[r]] <- mod
# } # add the three avoidance scores per run --

# ####### PLOTS #######

# df_vv <- data_long_runs %>% filter(ROI == "vv")
# # ROI values against CAPSB
# ggplot(df_vv, aes(x = CAPSB, y = value, color = run)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = TRUE, size = 1.2) +
#   labs(
#     title = "Region of interest values for contrast of own > control cue in vmPFC vs CAPS B subscale",
#     x = "CAPS B subscale score",
#     y = "ROI contrast value",
#     color = "Run"
#   ) +
#   theme_minimal(base_size = 14)


# df_al <- data_long_runs %>% filter(ROI == "al")
# # ROI values against CAPSB
# ggplot(df_al, aes(x = CAPSB, y = value, color = run)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = TRUE, size = 1.2) +
#   labs(
#     title = "Region of interest values for contrast of own > control cue in left amygdala vs CAPS B subscale",
#     x = "CAPS B subscale score",
#     y = "ROI contrast value",
#     color = "Run"
#   ) +
#   theme_minimal(base_size = 14)

# df_ar <- data_long_runs %>% filter(ROI == "ar")
# # ROI values against CAPSB
# ggplot(df_ar, aes(x = CAPSB, y = value, color = run)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = TRUE, size = 1.2) +
#   labs(
#     title = "Region of interest values for contrast of own > control cue in right amygdala vs CAPS B subscale",
#     x = "CAPS B subscale score",
#     y = "ROI contrast value",
#     color = "Run"
#   ) +
#   theme_minimal(base_size = 14)

# df_hl <- data_long_runs %>% filter(ROI == "hl")
# # ROI values against CAPSB
# ggplot(df_hl, aes(x = CAPSB, y = value, color = run)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = TRUE, size = 1.2) +
#   labs(
#     title = "Region of interest values for contrast of own > control cue in left hippocampus vs CAPS B subscale",
#     x = "CAPS B subscale score",
#     y = "ROI contrast value",
#     color = "Run"
#   ) +
#   theme_minimal(base_size = 14)

# df_hr <- data_long_runs %>% filter(ROI == "hr")
# # ROI values against CAPSB
# ggplot(df_hr, aes(x = CAPSB, y = value, color = run)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = TRUE, size = 1.2) +
#   labs(
#     title = "Region of interest values for contrast of own > control cue in right hippocampus vs CAPS B subscale",
#     x = "CAPS B subscale score",
#     y = "ROI contrast value",
#     color = "Run"
#   ) +
#   theme_minimal(base_size = 14)


# #betas afterwards for plotting only, use contrasts for stats
# #bar for own and control, having two bars rather than contrast single bars - mean of the betas not same as mean of diff but wtv
