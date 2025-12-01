library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(broom.mixed)
library(purrr)

data <- read.csv("/Users/zefreeman/Documents/All_ROIs_20251017.csv", header = TRUE, stringsAsFactors = FALSE)

newnames <- c(
  "subject" = 1, "site" = 3, "age" = 8, "CAPSB" = 9,
  "run1vv" = 10, "run2vv" = 11, "run3vv" = 12, "vmpfc_v_total" = 13,
  "run1vd" = 22, "run2vd" = 23, "run3vd" = 24, "vmpfc_d_total" = 25,
  "run1al" = 34, "run2al" = 35, "run3al" = 36, "amyg_l_total" = 37,
  "run1ar" = 46, "run2ar" = 47, "run3ar" = 48, "amyg_r_total" = 49,
  "run1hl" = 58, "run2hl" = 59, "run3hl" = 60, "hipp_l_total" = 61,
  "run1hr" = 70, "run2hr" = 71, "run3hr" = 72, "hipp_r_total" = 73,
  "avo_1" = 90, "diss_1" = 91, "anx_1" = 92,
  "avo_2" = 93, "diss_2" = 94, "anx_2" = 95,
  "avo_3" = 96, "diss_3" = 97, "anx_3" = 98,
  "now_own_1" = 101, "now_own_2" = 102, "now_own_3" = 103,
  "now_own_mean" = 104, "STAR_ID" = 105
)
names(data)[newnames] <- names(newnames)

# colnames(data)[1] <- "subject"
# colnames(data)[3] <- "site"
# colnames(data)[8] <- "age"
# colnames(data)[9] <- "CAPSB"
# colnames(data)[10] <- "run1vv"
# colnames(data)[11] <- "run2vv"
# colnames(data)[12] <- "run3vv"
# colnames(data)[13] <- "vmpfc_v_total"
# colnames(data)[22] <- "run1vd" # not using
# colnames(data)[23] <- "run2vd" # not using
# colnames(data)[24] <- "run3vd" # not using
# colnames(data)[25] <- "vmpfc_d_total" # not using
# colnames(data)[34] <- "run1al"
# colnames(data)[35] <- "run2al"
# colnames(data)[36] <- "run3al"
# colnames(data)[37] <- "amyg_l_total"
# colnames(data)[46] <- "run1ar"
# colnames(data)[47] <- "run2ar"
# colnames(data)[48] <- "run3ar"
# colnames(data)[49] <- "amyg_r_total"
# colnames(data)[58] <- "run1hl"
# colnames(data)[59] <- "run2hl"
# colnames(data)[60] <- "run3hl"
# colnames(data)[61] <- "hipp_l_total"
# colnames(data)[70] <- "run1hr"
# colnames(data)[71] <- "run2hr"
# colnames(data)[72] <- "run3hr"
# colnames(data)[73] <- "hipp_r_total"
# colnames(data)[90] <- "avo_1"
# colnames(data)[91] <- "diss_1"
# colnames(data)[92] <- "anx_1"
# colnames(data)[93] <- "avo_2"
# colnames(data)[94] <- "diss_2"
# colnames(data)[95] <- "anx_2"
# colnames(data)[96] <- "avo_3"
# colnames(data)[97] <- "diss_3"
# colnames(data)[98] <- "anx_3"
# colnames(data)[101] <- "now_own_1"
# colnames(data)[102] <- "now_own_2"
# colnames(data)[103] <- "now_own_3"
# colnames(data)[104] <- "now_own_mean" # mean own memory nownness score across runs
# colnames(data)[105] <- "STAR_ID"

# does avoidance score explain variation in ROI signal?

# wide to long for run-specific ROIs
data_long_runs <- data %>%
  pivot_longer(
    cols = c(run1vv, run2vv, run3vv, 
             run1al, run2al, run3al, 
             run1ar, run2ar, run3ar, 
             run1hl, run2hl, run3hl, 
             run1hr, run2hr, run3hr),
    names_to = c("run", "ROI"),
    names_pattern = "run([0-9]+)([a-zA-Z]+)",
    values_to = "value"
  ) %>%
  mutate(run = as.character(run)) %>%
  filter(!is.na(subject), !is.na(value)) %>%   # remove rows with missing subject or value
  select(subject, site, CAPSB, age, run, ROI, value, now_own_mean, STAR_ID)
  
# wide to long for avoidance scores
data_long_avoid <- data %>%
  pivot_longer(
    cols = c(avo_1, avo_2, avo_3),
    names_to = "run",
    names_pattern = "avo_([0-9]+)",
    values_to = "avoid_score"
  ) %>%
    mutate(run = as.character(run)) %>%
  select(subject, run, avoid_score)

# Make sure 'run' columns are character
data_long_runs <- data_long_runs %>%
  mutate(run = as.character(run))

data_long_avoid <- data_long_avoid %>%
  mutate(run = as.character(run)) %>%
  distinct(subject, run, .keep_all = TRUE)  # remove any duplicates

# Left join: replicate avoid_score across all ROI rows per subject/run
data_long_runs <- data_long_runs %>%
  left_join(data_long_avoid, by = c("subject", "run"))


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

# wide to long for nowness
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


-----------------------------------------------------------------------------------
### 1. LME per ROI without CAPS and without random slopes, overall not by run
----------------------------------------------------------------------------------- # nolint

results <- data_long_runs %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map_dfr(function(data_long_runs) {
    m <- try(lmer(roi_signal ~ avoid_score + site + (1 | subject), data = data_long_runs), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    tidy(m) %>%
      filter(term == "avoid_score")
  }, .id = "ROI")
results



#######################################################################################
--------------------------------------------------------------------------
# 2. mixed effects with random slope for run, stepwise across multiple USING THESE <<<<<<<<<<<<<<<<<<<<<<<
# ---------------------------------------------------------------------------

data_long_runs <- data_long_runs %>% # nolint
  mutate(run = as.numeric(run))

library(performance)

get_vif_lmer <- function(model) {
  performance::check_collinearity(model)
}
# 20251021 amygdala combined model

# Combine left/right amygdala into one dataset with hemi variable
amygdala_data <- data_long_runs %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                             REPORTING THIS FOR AMYGDALA
results_amygdala <- amygdala_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: amygdala\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * CAPSB + (1 | subject),
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
map(results_amygdala, summary)



library(broom.mixed)
library(emmeans)
m1_amygdala <- results_amygdala$amygdala
# Compute estimated marginal means (predicted means per run × site)
emm_amygdala <- emmeans(m1_amygdala, ~ site * run)
aplot_data <- as.data.frame(emm_amygdala)

a_plot <- ggplot(aplot_data, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  coord_cartesian(ylim = c(-0.75, 1.5)) +  # <-- y-axis limits
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Amygdala: Predicted ROI Signal by Run and Site",
    x = "Run",
    y = "Fitted Value (± SEM)",
    fill = "Site"
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )
ggsave("amyg_run*site.png", a_plot, width = 6, height = 4, dpi = 300)



# Extract fixed effects from all ROIs
results_table <- map_dfr(results_amygdala, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_amygdala.csv", row.names = FALSE)


# Combine left/right hippocampus into one dataset with hemi variable
hipp_data <- data_long_runs %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                                   REPORTING THIS FOR HIPPOCAMPUS
results_hipp <- hipp_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * CAPSB + (1 | subject),
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

m1_hipp<- results_hipp$hippocampus
# Compute estimated marginal means (predicted means per run × site)
emm_hipp <- emmeans(m1_hipp, ~ site * run)
hplot_data <- as.data.frame(emm_hipp)

h_plot <- ggplot(hplot_data, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  coord_cartesian(ylim = c(-0.5, 1)) +  # <-- y-axis limits
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Hippocampus: Predicted ROI Signal by Run and Site",
    x = "Run",
    y = "Fitted Value (± SEM)",
    fill = "Site"
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("hipp_run*site.png", h_plot, width = 6, height = 4, dpi = 300)





# Extract fixed effects from all ROIs
results_table <- map_dfr(results_hipp, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_hipp.csv", row.names = FALSE)






#-----------------------------------------------------------------------------
# SIMPLEST, no interaction, no avoidance
# lmer(roi_signal ~ CAPSB + run + site + age + (1 + subject),

results_m1 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m1 <- try(
      lmer(roi_signal ~ site + age + run * CAPSB + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m1)
  })
#results_m1
map(results_m1, summary)

vif_results <- results_m1 %>%
  compact() %>%  # remove NULLs if any
  map(get_vif_lmer)

#------------------------------------------------------------------------------
# SIMPLEST INCLUDING INTERACTION W/CAPS BY RUN, no avoidance                          REPORTING THIS FOR VMPFC
#      lmer(roi_signal ~ CAPSB * run + site + age + (1 | subject),

results_m2 <- data_long_runs %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs$ROI)) %>%
  map(function(data_long_runs) {
    cat("Running ROI:", unique(data_long_runs$ROI), "\n")
    m2 <- try(
      lmer(roi_signal ~ CAPSB * run + site + age + (1 | subject),
           data = data_long_runs, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m2, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m2)
  })
#results_m2
map(results_m2, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_m2, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_m2.csv", row.names = FALSE)


#------------------------------------------------------------------------------
# plots of these
#------------------------------------------------------------------------------
#################################

vv_data <- subset(data_long_runs, ROI == "vv") %>%
  mutate(run = as.factor(run))

v#v_data<-subset(vv_data, treatment)
# the mixed model for vv ROI
m1_vv <- lmer(value ~ run + site + age + (1 | subject), # * CAPSB
              data = vv_data, REML = FALSE)
summary(m1_vv)

# 3️⃣ Compute predicted marginal means for run × site
vv_emm <- emmeans(m1_vv, ~ site * run)
vv_plot_data <- as.data.frame(vv_emm)

# 4️⃣ Plot (same style as hippocampus/amygdala)
vv_plot <- ggplot(vv_plot_data, aes(x = run, y = emmean, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  coord_cartesian(ylim = c(-0.75, 0.5)) +  # y-axis range
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "vmPFC: Predicted ROI Signal by Run and Site",
    x = "Run",
    y = "Fitted Value (± SEM)",
    fill = "Site"
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# 5️⃣ Save the plot
ggsave("vmpfc_run*site.png", vv_plot, width = 6, height = 4, dpi = 300)



# ---------------------------
# 3. Grouped bar plot: Run × Site
# ---------------------------
run_site_summary <- data_fitted %>%
  group_by(run, site) %>%
  summarise(meanF = mean(fitted, na.rm = TRUE),
            semF  = sem(fitted)) %>%
  ungroup()

plot_run_site <- ggplot(run_site_summary, aes(x = run, y = meanF, fill = site)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = meanF - semF, ymax = meanF + semF),
                position = position_dodge(width = 0.8), width = 0.2) +
  ylab("Fitted ROI signal (mean ± SEM)") + xlab("Run") +
  theme_minimal()
ggsave("bar_Run_by_Site.png", plot_run_site, width = 6, height = 4, dpi = 300)

###########


#------------------------------------------------------------------------------           REPORTING NUMBER OF VOXELS PRESENT AS A COVARIATE


# 1. Read your CSV with the overlap data
roi_overlap <- read.csv("/Users/zefreeman/Documents/roi_all_overlap_summary_20251120.csv")
# Make sure it has columns like: subject, ROI_long, overlap

# 2. Map long ROI names to short codes
roi_map <- tribble(
  ~ROI, ~ROI_short,
  "lamygdala", "al",
  "ramygdala", "ar",
  "lhippocampus", "hl",
  "rhippocampus", "hr",
  "vmpfc", "vv"
)
roi_overlap <- roi_overlap %>%
  left_join(roi_map, by = "ROI")

data_long_runs_voxels <- data_long_runs %>%
  left_join(roi_overlap %>% select(SubjectID, ROI_short, ProportionOverlap), 
            by = c("subject" = "SubjectID", "ROI" = "ROI_short"))

results_voxels <- data_long_runs_voxels %>%
  mutate(run = as.factor(run)) %>%
  group_split(ROI) %>%
  set_names(unique(data_long_runs_voxels$ROI)) %>%
  map(function(data_long_runs_voxels) {
    cat("Running ROI:", unique(data_long_runs_voxels$ROI), "\n")
    m <- try(
      lmer(value ~ CAPSB * run + site + age + ProportionOverlap + (1 | subject),
           data = data_long_runs_voxels, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m, "try-error")) {
      cat("failed for:", unique(data_long_runs_voxels$ROI), "\n")
      return(NULL)
    }
    cat("succeeded for:", unique(data_long_runs_voxels$ROI), "\n")
    return(m)
  })
#results_m2
map(results_voxels, summary)

# Extract fixed effects from all ROIs
results_table <- map_dfr(results_voxels, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column

# Write to CSV
write.csv(results_table, "results_vv_voxels.csv", row.names = FALSE)

hipp_data_voxels <- data_long_runs_voxels %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                                   REPORTING THIS FOR HIPPOCAMPUS
results_hipp_voxels <- hipp_data_voxels %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * CAPSB + ProportionOverlap + (1 | subject),
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
map(results_hipp_voxels, summary)

results_table <- map_dfr(results_hipp_voxels, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column
write.csv(results_table, "results_h_voxels.csv", row.names = FALSE)

amyg_data_voxels <- data_long_runs_voxels %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run)
  )

# Run the mixed model                                                                   REPORTING THIS FOR HIPPOCAMPUS
results_amyg_voxels <- amyg_data_voxels %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    cat("Running ROI: a\n")
    m1 <- try(
      lmer(value ~ hemi + site + age + run * CAPSB + ProportionOverlap + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    if (inherits(m1, "try-error")) {
      cat("failed for: h\n")
      return(NULL)
    }
    cat("succeeded for: a\n")
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_amyg_voxels, summary)

results_table <- map_dfr(results_amyg_voxels, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column
write.csv(results_table, "results_a_voxels.csv", row.names = FALSE)




#------------------------------------------------------------------------------           REPORTING TRIAL BY TRIAL NOWNESS
# nowness by trial for vmpfc

trial_data <- read.csv("/Users/zefreeman/Documents/T1_nowness_ratings.csv", header = TRUE, stringsAsFactors = FALSE)

# Filter to only "Yes" trials
df_yes <- trial_data %>%
  filter(IsOwn == "Yes")

  # Filter to other trials
df_other <- trial_data %>%
  filter(IsOwn == "No")

# Basic histogram of raw ratings across all trials, runs, participants                    HISTOGRAM OF RAW OWN NOWNESS RATINGS REPORTING
nowrate <- ggplot(df_yes, aes(x = Rating)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Raw Own Nowness Ratings",
    x = "Raw Rating",
    y = "Count"
  )
ggsave(filename = "nowness_ratings_own.png", plot = nowrate, width = 8, height = 6, dpi = 300)

summary_table_wide <- trial_data %>%
  group_by(ScanSite, Run, IsOwn) %>%
  summarise(
    mean = mean(Rating, na.rm = TRUE),
    sd   = sd(Rating, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean, sd)) %>%
  select(ScanSite, Run, IsOwn, mean_sd) %>%
  tidyr::pivot_wider(
    names_from = Run,
    values_from = mean_sd,
    names_prefix = "Run "
  )
summary_table_wide
write.csv(summary_table_wide, "T1_trialbynowness.csv", row.names = FALSE)


summary_table_long <- df_yes %>%
  filter(!is.na(Run), Run != "NA") %>%      # Remove NA or "NA" run values
  group_by(ScanSite, Run) %>%
  summarise(
    mean = mean(Rating, na.rm = TRUE),
    sd   = sd(Rating, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean), !is.na(sd)) 

meannow <- ggplot(summary_table_long,                 #   REPORTING THIS PLOT FOR OWN NOWNESS RATING MEAN SCORES PER RUN
       aes(x = factor(Run), y = mean, fill = ScanSite)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2,
                position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean Own Nowness Rating by Run and Site (± SD)",
    x = "Run",
    y = "Mean Rating",
    fill = "Site"
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

ggsave(filename = "m_nowness_ratings_own.png", plot = meannow, width = 8, height = 6, dpi = 300)


########################################
## CORRELATION OF NOWNESS WITH CAPS AND AVOIDANCE SCORES
#########################################

nowness_run_avgs <- trial_data %>%              #
  filter(!is.na(Run), Run != "NA") %>%
  group_by(subject = STAR_ID, run = Run) %>%
  summarise(
    mean_nowness = mean(Rating, na.rm = TRUE),
    .groups = "drop"
  )
avoid_caps_runs <- data_long_runs %>%
  filter(run %in% c(1,2,3)) %>%
  group_by(STAR_ID, run) %>%
  summarise(
    mean_avoid = mean(avoid_score, na.rm = TRUE),
    mean_CAPS  = mean(CAPSB, na.rm = TRUE),
    .groups = "drop"
  )

merged_runs <- nowness_run_avgs %>%
  left_join(avoid_caps_runs, by = c("subject" = "STAR_ID", "run" = "run"))
vars_for_cor <- merged_runs %>%
  select(subject, run, mean_nowness, mean_avoid, mean_CAPS)
cor_by_run <- function(run_num) {
  df <- vars_for_cor %>%
    filter(run == run_num) %>%
    select(mean_nowness, mean_avoid, mean_CAPS)

  cor(df, use = "complete.obs")
}
corr_run1 <- cor_by_run(1)
corr_run2 <- cor_by_run(2)
corr_run3 <- cor_by_run(3)

# # Create a function to extract correlations in long format
# cor_to_long <- function(cor_matrix, run_number) {
#   melt(cor_matrix) %>%
#     mutate(run = paste0("Run ", run_number))
# }
# library(reshape2)
# # Convert all three runs
# cor_long_all <- bind_rows(
#   cor_to_long(corr_run1, 1),
#   cor_to_long(corr_run2, 2),
#   cor_to_long(corr_run3, 3)
# )
# heatmap<-ggplot(cor_long_all, aes(Var1, Var2, fill = value)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = sprintf("%.2f", value)), size = 4) +
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = c(-1, 1),
#     name = "Correlation"
#   ) +
#   facet_wrap(~ run, nrow = 1) +    # side-by-side like a time series
#   theme_minimal(base_size = 11) +
#   labs(
#     title = "Correlation Matrices Across Runs (Nowness, Avoidance, CAPS-B)",
#     x = "",
#     y = ""
#   ) +
#   theme(
#     strip.text = element_text(size = 14, face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   )
# ggsave(filename = "heatmap_own.png", plot = heatmap, width = 8, height = 6, dpi = 300)




# # Assign Run based on trial ranges
# df_change <- df_yes %>%
#   group_by(STAR_ID, ScanSite) %>%
#   arrange(Trial, .by_group = TRUE) %>%
#   mutate(
#     baseline = first(Rating),
#     change_from_baseline = Rating - baseline,
#     Run = case_when(
#       Trial <= 14 ~ "Run 1",
#       Trial >= 15 & Trial <= 32 ~ "Run 2",
#       Trial >= 33 ~ "Run 3"
#     )
#   ) %>%
#   ungroup() %>%
#   mutate(Run = factor(Run, levels = c("Run 1", "Run 2", "Run 3")))
# # Define site-specific base colors and shading per run
# site_colours <- list(
#   "Newcastle" = c("Run 1" = "#1f77b4", "Run 2" = "#6baed6", "Run 3" = "#c6dbef"),
#   "Manchester" = c("Run 1" = "#2ca02c", "Run 2" = "#1ca02c", "Run 3" = "#98df8a"),
#   "London" = c("Run 1" = "#d62728", "Run 2" = "#ff9896", "Run 3" = "#e7969c")
# )
# # Function to generate plot per site
# make_site_plot <- function(site_name) {

#   df_site <- df_change %>% filter(ScanSite == site_name)

#   ggplot(df_site,
#          aes(x = Trial,
#              y = change_from_baseline,
#              group = STAR_ID,
#              color = factor(Run))) +
#     geom_line(alpha = 0.6) +
#     geom_point(alpha = 0.7) +
#     scale_color_manual(values = site_colours[[site_name]]) +
#     labs(
#       title = paste("Nowness Rating Change From Baseline –", site_name),
#       x = "Trial",
#       y = "Change From Baseline (Rating)",
#       color = "Run"
#     ) +
#     theme_minimal(base_size = 14)
# }
# # Create the three plots
# plot_london     <- make_site_plot("London")
# plot_manchester <- make_site_plot("Manchester")
# plot_newcastle  <- make_site_plot("Newcastle")
# ggsave(filename = "T1_nowness_london.png", plot = plot_london, width = 8, height = 6, dpi = 300)
# ggsave(filename = "T1_nowness_manchester.png", plot = plot_manchester, width = 8, height = 6, dpi = 300)
# ggsave(filename = "T1_nowness_newcastle.png", plot = plot_newcastle, width = 8, height = 6, dpi = 300)

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

model_vtrial <- lmer(vv_value ~ CAPSB  * run + Rating + Trial + site + age + (1 | subject),
                    data = trial_data_vv)
summary(model_vtrial)

model_vtrialbynow <- lmer(vv_value ~ Trial  * run + Rating + CAPSB + site + age + (1 | subject),
                    data = trial_data_vv)
summary(model_vtrialbynow)

results_table <- tidy(model_vtrialbynow) %>%
  select(term, estimate, std.error, df, statistic, p.value)

# Write to CSV
write.csv(results_table, "results_vmpfc_trialbynowness.csv", row.names = FALSE)

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

##############################################################################################################
##############################################################################################################

# AMYGDALA VERSION WITH HEMISPHERES                                               REPORTING TRIAL AMYG.
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

model_arl <- lmer(
  amyg_value ~ CAPSB * run + Rating + Trial + site + age + Hemisphere + (1 | subject),
  data = trial_data_arl_long
)
summary(model_arl)

model_arlbynow <- lmer(
  amyg_value ~ Trial * run + Rating + CAPSB + site + age + Hemisphere + (1 | subject),
  data = trial_data_arl_long
)
summary(model_arlbynow)

results_table <- tidy(model_arlbynow) %>%
  select(term, estimate, std.error, df, statistic, p.value)

# Write to CSV
write.csv(results_table, "results_amygdala_trialbynowness.csv", row.names = FALSE)



# SAME TRIAL BY TRIAL INCLUSION FOR HIPPOCAMPUS WITH HEMISPHERES                             REPORTING THIS FOR HIPPOCAMPUS TRIAL

hr_data <- data_long_runs %>%
  filter(ROI == "hr") %>%                    # keep only hr
  select(subject, site, CAPSB, age, run, value, STAR_ID) %>%  # select desired columns
  rename(
    hr_value = value,                        
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer

hr_data <- hr_data %>%
  distinct(STAR_ID, Run, .keep_all = TRUE)

hl_data <- data_long_runs %>%
  filter(ROI == "hl") %>%                    # keep only hr
  select(subject, site, CAPSB, age, run, value, STAR_ID) %>%  # select desired columns
  rename(
    hl_value = value,                        
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer

hl_data <- hl_data %>%
  distinct(STAR_ID, Run, .keep_all = TRUE)

# Merge HR AND HL 
hrl_data <- hr_data %>%
  left_join(hl_data, by = c("STAR_ID", "Run", "subject", "site", "CAPSB", "age"), relationship = "many-to-many")

hrl_data <- hrl_data %>%
  mutate(Run = as.integer(Run))

trial_data_hrl <- trial_data %>%
  left_join(hrl_data, by = c("STAR_ID" = "STAR_ID", "Run" = "Run"))


trial_data_hrl_long <- trial_data_hrl %>%
  pivot_longer(
    cols = c(hl_value, hr_value),
    names_to = "Hemisphere",
    values_to = "hipp_value"
  ) %>%
  mutate(
    Hemisphere = factor(Hemisphere, levels = c("hl_value", "hr_value")),  # left as reference
    CAPSB = as.numeric(CAPSB),
    run = as.factor(Run),
    site = as.factor(site)
  )

model_hrl <- lmer(
  hipp_value ~ CAPSB * run + Rating + Trial + site + age + Hemisphere + (1 | subject),
  data = trial_data_hrl_long
)
summary(model_hrl)

model_hrlbynow <- lmer(
  hipp_value ~ Trial * run + Rating + CAPSB + site + age + Hemisphere + (1 | subject),
  data = trial_data_hrl_long
)
summary(model_hrlbynow)

model_hrl_ratingbyrun <- lmer(
hipp_value ~ Rating * run + CAPSB + Trial + site + age + Hemisphere + (1 | subject),
  data = trial_data_hrl_long
)
summary(model_hrl_ratingbyrun)

results_table <- tidy(model_hrl_ratingbyrun) %>%
  select(term, estimate, std.error, df, statistic, p.value)

# Write to CSV
write.csv(results_table, "results_hippocampus_runbynowness.csv", row.names = FALSE)


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
### 4. scatter plots
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




