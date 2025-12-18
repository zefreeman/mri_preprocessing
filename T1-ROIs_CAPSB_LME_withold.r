# ---------------------------------------------------------------
# Ze Freeman, 20251017 ----

install.packages("pacman")
library(pacman)
pacman::p_load(tidyr, tidyverse, dplyr, lme4, lmerTest, performance,
               ggplot2, broom.mixed, purrr, emmeans) #ggpredict

# to do:
# - why did some people not run in the baseline group level model?

# 0. Data loading -------------------------------------------------------------

data <- read.csv("/Users/zefreeman/Documents/All_ROIs_20251017.csv",
                            header = TRUE, stringsAsFactors = FALSE)
newnames <- c(
  "subject" = 1, "site" = 3, "age" = 8, "CAPSB" = 9,
  # "run1vv" = 10, "run2vv" = 11, "run3vv" = 12, "vmpfc_v_total" = 13,
  # "run1vd" = 22, "run2vd" = 23, "run3vd" = 24, "vmpfc_d_total" = 25,
  # "run1al" = 34, "run2al" = 35, "run3al" = 36, "amyg_l_total" = 37,
  # "run1ar" = 46, "run2ar" = 47, "run3ar" = 48, "amyg_r_total" = 49,
  # "run1hl" = 58, "run2hl" = 59, "run3hl" = 60, "hipp_l_total" = 61,
  # "run1hr" = 70, "run2hr" = 71, "run3hr" = 72, "hipp_r_total" = 73,
  "avo_1" = 90, "diss_1" = 91, "anx_1" = 92,
  "avo_2" = 93, "diss_2" = 94, "anx_2" = 95,
  "avo_3" = 96, "diss_3" = 97, "anx_3" = 98,
  "now_own_1" = 101, "now_own_2" = 102, "now_own_3" = 103,
  "now_own_mean" = 104, "STAR_ID" = 105,
  "run1vv" = 106, "run2vv" = 107, "run3vv" = 108,
  "run1al" = 109, "run2al" = 110, "run3al" = 111,
  "run1ar" = 112, "run2ar" = 113, "run3ar" = 114,
  "run1hl" = 115, "run2hl" = 116, "run3hl" = 117,
  "run1hr" = 118, "run2hr" = 119, "run3hr" = 120)
names(data)[newnames] <- names(newnames)

trial_data <- read.csv("/Users/zefreeman/Documents/T1_nowness_ratings.csv",
                                    header = TRUE, stringsAsFactors = FALSE)

roi_overlap <- read.csv("/Users/zefreeman/Documents/roi_all_overlap_summary_20251120.csv")
roi_map <- tribble(
  ~ROI, ~ROI_short,
  "lamygdala", "al",
  "ramygdala", "ar",
  "lhippocampus", "hl",
  "rhippocampus", "hr",
  "vmpfc", "vv")
roi_overlap <- roi_overlap %>%
  left_join(roi_map, by = "ROI")


# 0. Reshaping main data ------------------------------------------------------

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
    values_to = "roi_signal") %>%
  mutate(run = as.character(run)) %>%
  filter(!is.na(subject), !is.na(roi_signal)) %>%  # remove rows w missing p or val
  select(subject, site, CAPSB, age, run, ROI, roi_signal, now_own_mean, STAR_ID)

# wide to long for avoidance scores
data_long_avoid <- data %>%
  pivot_longer(
    cols = c(avo_1, avo_2, avo_3),
    names_to = "run",
    names_pattern = "avo_([0-9]+)",
    values_to = "avoid_score") %>%
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

# data_long_runs <- data_long_runs %>%
#   mutate(run = as.character(run)) %>%
#   left_join(data_long_now, by = c("subject", "run"))


########################################
## CORRELATION OF NOWNESS WITH CAPS AND AVOIDANCE SCORES
#########################################

nowness_run_avgs <- trial_data %>%
  filter(!is.na(Run), Run != "NA") %>%
  group_by(subject = SubjectID, run = Run) %>%
  summarise(
    mean_nowness = mean(Rating, na.rm = TRUE),
    .groups = "drop")
avoid_caps_runs <- data_long_runs %>%
  filter(run %in% c(1,2,3)) %>%
  group_by(STAR_ID, run) %>%
  summarise(
    mean_avoid = mean(avoid_score, na.rm = TRUE),
    mean_CAPS  = mean(CAPSB, na.rm = TRUE),
    .groups = "drop")
nowness_run_avgs <- nowness_run_avgs %>%
  mutate(run = as.numeric(run))

avoid_caps_runs <- avoid_caps_runs %>%
  mutate(run = as.numeric(run))
merged_runs <- nowness_run_avgs %>%
  left_join(avoid_caps_runs, by = c("subject" = "STAR_ID", "run" = "run"))

vars_for_cor <- merged_runs %>%
  select(subject, run, mean_nowness, mean_avoid, mean_CAPS)
cor_by_run <- function(run_num) {
  df <- vars_for_cor %>%
    filter(run == run_num) %>%
    select(mean_nowness, mean_avoid, mean_CAPS)

  cor(df, use = "complete.obs")}
corr_run1 <- cor_by_run(1)
corr_run2 <- cor_by_run(2)
corr_run3 <- cor_by_run(3)

# correlation heatmap by run with correlation values and all runs
corr_long <- bind_rows(
  melt(corr_run1) %>% mutate(run = "Run 1"),
  melt(corr_run2) %>% mutate(run = "Run 2"),
  melt(corr_run3) %>% mutate(run = "Run 3"))
ggplot(corr_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.2f", value)),
    size = 3) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation") +
  facet_wrap(~ run, nrow = 1) +
  labs(
    x = "",
    y = "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold"))
ggsave("correlation_nowness_avoid_CAPS_by_run.png", width = 10, height = 4, dpi = 300)


# predicting avoidance with caps by run interaction ------------------------------
data_long_runs <- data_long_runs %>%
  mutate(
    run = factor(run, levels = c("1", "2", "3"))
  )

m_avoid <- lmer(
  avoid_score ~ run * CAPSB + site + age + (1 | subject),
  data = data_long_runs,
  REML = FALSE)

summary(m_avoid)
write.csv(
  broom.mixed::tidy(m_avoid) %>%
    select(term, estimate, std.error, df, statistic, p.value),
  "results_avoid_predicted_by_CAPSbyrun.csv",
  row.names = FALSE
)
# mixed effects w random slope for run --------------------------- FIRST MODEL A

data_long_runs <- data_long_runs %>% 
  mutate(run = as.numeric(run))

get_vif_lmer <- function(model) {
  performance::check_collinearity(model)}
# 20251021 amygdala combined model

# Combine left/right amygdala into one dataset with hemi variable
amygdala_data <- data_long_runs %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run))

results_m1_amyg <- list(
  amygdala = lmer(
    roi_signal ~ hemi + site + age + run * CAPSB + (1 | subject),
    data = amygdala_data,
    REML = FALSE))
map(results_m1_amyg, summary)

results_m1_amyg_table <- map_dfr(results_m1_amyg, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column
write.csv(results_m1_amyg_table, "results_m1_amyg.csv", row.names = FALSE)

plot_data <- augment(results_m1_amyg$amygdala) %>%
  mutate(
    run  = factor(run),
    site = factor(site))

a_plot <- ggplot(plot_data, aes(x = run, y = .fitted, fill = site)) +
  geom_bar(
    stat = "summary",
    fun = mean,
    position = position_dodge(width = 0.8),
    width = 0.7) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.2) +
  coord_cartesian(ylim = c(-0.25, 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Amygdala: Predicted ROI Signal by Run and Site",
    x = "Run",
    y = "Fitted Value (± SEM)",
    fill = "Site") +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank())
ggsave("amyg_run*site.png", a_plot, width = 6, height = 4, dpi = 300)

# plot the run effects collapsing by site ---------------------------------
a2_plot <- ggplot(plot_data, aes(x = run, y = .fitted)) +
  geom_bar(
    stat = "summary",
    fun = mean,
    width = 0.7,
    fill = "skyblue") +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    width = 0.2) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Amygdala: Predicted ROI Signal by Run (Collapsed Across Sites)",
    x = "Run",
    y = "Fitted Value (± SEM)") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank())
ggsave("amyg_run.png", a2_plot, width = 6, height = 4, dpi = 300)



# scatter plot of fitted roi by CAPS Z scored -------------------------------
mf <- model.frame(m1_amygdala)
mf <- mf %>%
  mutate(
    CAPSB_z = as.numeric(scale(CAPSB)),
    fitted_roi = fitted(m1_amygdala),
    run = factor(run))
ggplot(mf, aes(x = CAPSB_z, y = fitted_roi, colour = run)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Amygdala: Fitted ROI Signal by CAPS-B and Run",
    x = "CAPS-B (z-scored)",
    y = "Fitted ROI Signal",
    colour = "Run")
ggsave("amyg_roifitted_byCAPS_run.png", width = 8, height = 6, dpi = 300)


# same for hippocampus version with hemispheres ------------------ FIRST MODEL H

hipp_data <- data_long_runs %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run))

results_m1_hipp <- list(
  hippocampus = lmer(
    roi_signal ~ hemi + site + age + run * CAPSB + (1 | subject),
    data = hipp_data,
    REML = FALSE))
map(results_m1_hipp, summary)

results_m1_hipp_table <- map_dfr(results_m1_hipp, function(model) {
  if (is.null(model)) return(NULL)  
  tidy(model) %>%
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")
write.csv(results_m1_hipp_table, "results_m1_hipp.csv", row.names = FALSE)

plot_data <- augment(results_m1_hipp$hippocampus) %>%
  mutate(
    run  = factor(run),
    site = factor(site))

h_plot <- ggplot(plot_data, aes(x = run, y = .fitted, fill = site)) +
  geom_bar(
    stat = "summary",
    fun = mean,
    position = position_dodge(width = 0.8),
    width = 0.7) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.2) +
  coord_cartesian(ylim = c(-0.25, 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Hippocampus: Predicted ROI Signal by Run and Site",
    x = "Run",
    y = "Fitted Value (± SEM)",
    fill = "Site") +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank())
ggsave("hipp_run*site.png", h_plot, width = 6, height = 4, dpi = 300)

# plot the run effects collapsing by site ---------------------------------
h2_plot <- ggplot(plot_data, aes(x = run, y = .fitted)) +
  geom_bar(
    stat = "summary",
    fun = mean,
    width = 0.7,
    fill = "skyblue") +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    width = 0.2) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Hippocampus: Predicted ROI Signal by Run (Collapsed Across Sites)",
    x = "Run",
    y = "Fitted Value (± SEM)") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank())
ggsave("hipp_run.png", h2_plot, width = 6, height = 4, dpi = 300)

# scatter plot of fitted roi by CAPS Z scored -------------------------------
mf <- model.frame(m1_hipp)
mf <- mf %>%
  mutate(
    CAPSB_z = as.numeric(scale(CAPSB)),
    fitted_roi = fitted(m1_hipp),
    run = factor(run))

ggplot(mf, aes(x = CAPSB_z, y = fitted_roi, colour = run)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Hippocampus: Fitted ROI Signal by CAPS-B and Run",
    x = "CAPS-B (z-scored)",
    y = "Fitted ROI Signal",
    colour = "Run")
ggsave("hipp_roifitted_byCAPS_run.png", width = 8, height = 6, dpi = 300)



# vmpfc model --------------------------------------------------- FIRST MODEL V
vmpfc_data <- data_long_runs %>%
  filter(ROI == "vv") %>%
  mutate(run = as.factor(run))
results_m1_vmpfc <- lmer(
  roi_signal ~ site + age + run * CAPSB + (1 | subject),
  data = vmpfc_data,
  REML = FALSE)
summary(results_m1_vmpfc)

results_table_vmpfc <- broom.mixed::tidy(results_m1_vmpfc) %>%
  select(term, estimate, std.error, df, statistic, p.value) %>%
  mutate(ROI = "vv")
write.csv(results_table_vmpfc, "results_m1_vv.csv", row.names = FALSE)

vif_results <- results_m1 %>%
  compact() %>%  # remove NULLs if any
  map(get_vif_lmer)


# plot the fitted values ------------------------------------------------

plot_data <- augment(results_m1_vmpfc) %>%
  mutate(
    run  = factor(run),
    site = factor(site))

vv_plot <- ggplot(plot_data, aes(x = run, y = .fitted, fill = site)) +
  geom_bar(
    stat = "summary",
    fun = mean,
    position = position_dodge(width = 0.8),
    width = 0.7) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(width = 0.8),
    width = 0.2) +
  coord_cartesian(ylim = c(-0.25, 0.75)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "vmPFC: Predicted ROI Signal by Run and Site",
    x = "Run",
    y = "Fitted Value (± SEM)",
    fill = "Site") +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank())
ggsave("vmpfc_run*site.png", vv_plot, width = 6, height = 4, dpi = 300)

# plot the run effects collapsing by site ---------------------------------
vv2_plot <- ggplot(plot_data, aes(x = run, y = .fitted)) +
  geom_bar(
    stat = "summary",
    fun = mean,
    width = 0.7,
    fill = "skyblue") +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    width = 0.2) +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "vmPFC: Predicted ROI Signal by Run (Collapsed Across Sites)",
    x = "Run",
    y = "Fitted Value (± SEM)") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank())
ggsave("vmpfc_run.png", vv2_plot, width = 6, height = 4, dpi = 300)


# scatter plot of fitted roi by CAPS Z scored -------------------------------
mf <- model.frame(results_m1_vmpfc)
mf <- mf %>%
  mutate(
    CAPSB_z = as.numeric(scale(CAPSB)),
    fitted_roi = fitted(results_m1),
    run = factor(run))

ggplot(mf, aes(x = CAPSB_z, y = fitted_roi, colour = run)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "vmPFC: Fitted ROI Signal by CAPS-B and Run",
    x = "CAPS-B (z-scored)",
    y = "Fitted ROI Signal",
    colour = "Run")
ggsave("vmpfc_roifitted_byCAPS_run.png", width = 8, height = 6, dpi = 300)

# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_m1_vv.csv")
amyg_csv <- read.csv("results_m1_amyg.csv")
hipp_csv <- read.csv("results_m1_hipp.csv")

combined_results <- bind_rows(
  amyg_csv %>% mutate(model = "amygdala"),
  hipp_csv %>% mutate(model = "hippocampus"),
  vv_csv %>% mutate(model = "vmpfc"))

write.csv(combined_results, "results_m1_combined.csv", row.names = FALSE)




# avoidance mixed effects ------------------------------------- SECOND (AVOID) MODEL

vif_results3 <- results_m3 %>%
  compact() %>%  # remove NULLs if any
  map(get_vif_lmer)
vif_results3

#----------------------------------------------------------------------------
results_m2_vmpfc <- lmer(
  roi_signal ~ avoid_score * run + CAPSB + site + age + (1 | subject),
  data = vmpfc_data,
  REML = FALSE)
summary(results_m2_vmpfc)
m2_vmpfc <- results_m2_vmpfc

results_m2_vmpfc <- broom.mixed::tidy(results_m2_vmpfc) %>%
  select(term, estimate, std.error, df, statistic, p.value) %>%
  mutate(ROI = "vv")
write.csv(results_m2_vmpfc, "results_vmpfc_avoidbyrun.csv", row.names = FALSE)

# results_m3.5 <- data_long_runs %>%
#   mutate(run = as.factor(run)) %>%
#   group_split(ROI) %>%
#   set_names(unique(data_long_runs$ROI)) %>%
#   map(function(data_long_runs) {
#   cat("Running ROI:", unique(data_long_runs$ROI), "\n")
#   m3.5 <- lmer(roi_signal ~  avoid_score * run + CAPSB + site + age + (1 | subject),
#            data = data_long_runs, REML = FALSE),
#       silent = TRUE})
# map(results_m3.5, summary)

# vif_results3.5 <- results_m3.5 %>%
#   compact() %>%  # remove NULLs if any
#   map(get_vif_lmer)
# results_table <- map_dfr(results_m3.5, function(model) {
#   if (is.null(model)) return(NULL)  # skip failed models
#   tidy(model) %>% 
#     select(term, estimate, std.error, df, statistic, p.value)
# }, .id = "ROI")  
# # just the vmpfc vv column
# resultsm3.5vv_table <- tidy(results_m3.5$vv) %>%
#   select(term, estimate, std.error, df, statistic, p.value)
# write.csv(resultsm3.5vv_table, "results_vmpfc_avoidbyrun.csv", row.names = FALSE)


# scatter plot of fitted roi by CAPS Z scored -------------------------------
mf <- model.frame(m2_vmpfc)
mf <- mf %>%
  mutate(
    avoid_score = as.numeric(scale(avoid_score)),
    fitted_roi = fitted(m2_vmpfc),
    run = factor(run))

ggplot(mf, aes(x = avoid_score, y = fitted_roi, colour = run)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "vmPFC: Fitted ROI Signal by Avoidance Score and Run",
    x = "Avoidance Score (z-scored)",
    y = "Fitted ROI Signal",
    colour = "Run")
ggsave("vmpfc_roifitted_byAvoid_run.png", width = 8, height = 6, dpi = 300)




results_m2_hipp <- hipp_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    m2 <- try(
      lmer(roi_signal ~ avoid_score * run + CAPSB + site + age + hemi + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE
    )
    return(m2)
  }) %>%
  set_names("hippocampus")
map(results_m2_hipp, summary)
m2_hipp <- results_m2_hipp$hippocampus

results_m2_hipp <- map_dfr(results_m2_hipp, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")
write.csv(results_m2_hipp, "results_hipp_avoidbyrun.csv", row.names = FALSE)

# scatter plot of fitted roi by CAPS Z scored -------------------------------
mf <- model.frame(m2_hipp)
mf <- mf %>%
  mutate(
    avoid_score = as.numeric(scale(avoid_score)),
    fitted_roi = fitted(m2_hipp),
    run = factor(run))

ggplot(mf, aes(x = avoid_score, y = fitted_roi, colour = run)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Hippocampus: Fitted ROI Signal by Avoidance Score and Run",
    x = "Avoidance Score (z-scored)",
    y = "Fitted ROI Signal",
    colour = "Run")
ggsave("hipp_roifitted_byAvoid_run.png", width = 8, height = 6, dpi = 300)



results_m2_amygdala <- amygdala_data %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    m1 <- try(
      lmer(roi_signal ~ avoid_score * run + CAPSB + site + age + hemi + (1 | subject), 
      # originally run * caps
           data = df, REML = FALSE),
      silent = TRUE
    )
    return(m1)
  }) %>%
  set_names("amygdala")
map(results_m2_amygdala, summary)

m2_amyg <- results_m2_amygdala$amygdala

results_m2_amygdala <- map_dfr(results_m2_amygdala, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # add ROI as a column
write.csv(results_m2_amygdala, "results_amyg_avoidbyrun.csv", row.names = FALSE)

# scatter plot of fitted roi by CAPS Z scored -------------------------------
mf <- model.frame(m2_amyg)
mf <- mf %>%
  mutate(
    avoid_score = as.numeric(scale(avoid_score)),
    fitted_roi = fitted(m2_amyg),
    run = factor(run))

ggplot(mf, aes(x = avoid_score, y = fitted_roi, colour = run)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Amygdala: Fitted ROI Signal by Avoidance Score and Run",
    x = "Avoidance Score (z-scored)",
    y = "Fitted ROI Signal",
    colour = "Run")
ggsave("amyg_roifitted_byAvoid_run.png", width = 8, height = 6, dpi = 300)


# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_vmpfc_avoidbyrun.csv")
amyg_csv <- read.csv("results_amyg_avoidbyrun.csv")
hipp_csv <- read.csv("results_hipp_avoidbyrun.csv")

combined_results <- bind_rows(
  hipp_csv %>% mutate(model = "hippocampus"),
  amyg_csv %>% mutate(model = "amygdala"),
  vv_csv %>% mutate(model = "vmpfc")) %>%
  mutate(model = factor(model,
      levels = c("hippocampus", "amygdala", "vmpfc"))) %>%
  arrange(model)

write.csv(combined_results, "results_avoidbyrun_combined.csv", row.names = FALSE)

results <- read.csv("results_avoidbyrun_combined.csv")
term_labels <- c(
  "(Intercept)"        = "Intercept",
  "hemi"               = "Hemisphere",
  "Hemispherear_value" = "Hemisphere",
  "Hemispherehr_value" = "Hemisphere",
  "CAPSB"              = "CAPS-B main effect",
  "avoid_score"        = "Avoidance rating",
  "Rating"             = "Nowness rating",
  "trial"              = "Trial effect",
  "run2"               = "Run 2",
  "run3"               = "Run 3",
  "avoid_score:run2"   = "Run 2 x avoidance score",
  "avoid_score:run3"   = "Run 3 x avoidance score",
  "Rating:run2"        = "Run 2 x nowness rating",
  "Rating:run3"        = "Run 3 x nowness rating",
  "siteM"              = "Manchester",
  "siteN"              = "Newcastle",
  "age"                = "Age")
formatted_table <- combined_results %>%
  filter(!str_detect(term, "^sd__")) %>%
  mutate(
    ROI = case_when(
      model == "hippocampus" ~ "Hippocampus",
      model == "amygdala"    ~ "Amygdala",
      model == "vmpfc"       ~ "vmPFC"
    ),
    ROI = factor(
      ROI,
      levels = c("Hippocampus", "Amygdala", "vmPFC")
    ),
    term = recode(term, !!!term_labels),
    β  = round(estimate, 3),
    SE = round(std.error, 3),
    df = round(df, 3),
    t  = round(statistic, 3),
    p  = case_when(
      p.value < .001 ~ "<0.001***",
      p.value < .01  ~ "<0.01**",
      p.value < .05  ~ "<0.05*",
      TRUE           ~ as.character(round(p.value, 3))
    )
  ) %>%
  arrange(ROI) %>%   # <-- enforces final row order
  select(ROI, term, β, SE, df, t, p)
term_order <- c(
  "Intercept", "Hemisphere", "CAPS-B main effect", "Avoidance rating",
  "Trial effect", "Run 2", "Run 3", "Run 2 x avoidance score",
  "Run 3 x avoidance score", "Manchester", "Newcastle", "Age")

formatted_table <- formatted_table %>%
  mutate(term = factor(term, levels = term_order)) %>%
  arrange(ROI, term)
write.csv(formatted_table, "results_avoidbyrun_table_formatted.csv",
  row.names = FALSE)


#### GOT TO HERE 10:30PM TUESDAY 16TH DECEMBER ### ------------------------------------------------------------------------------

# trial by trial nowness -----------------------------------------
# descriptives

# Filter to only own ("Yes") trials
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

meannow <- ggplot(summary_table_long %>% 
  filter(Run %in% c("1", "2", "3")),  # Exclude "practice" runs
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

ggsave(filename = "m_nowness_ratings_own.png", 
plot = meannow, width = 8, height = 6, dpi = 300)


# TRIAL BY TRIAL NOWNESS ----------------------------------------------- third (nowness) model
trial_data <- trial_data %>%
  mutate(Run = as.integer(Run)) %>%
  rename(STAR_ID = SubjectID)

vv_data <- data_long_runs %>%
  filter(ROI == "vv") %>%                    # keep only vv ROI
  select(subject, site, CAPSB, age, run, roi_signal, STAR_ID) %>%  
  rename(
    vv_value = roi_signal,                        
    STAR_ID = STAR_ID,                     
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

model_vtrial <- lmer(vv_value ~ CAPSB * run + Rating + site + age + (1 | subject),
                    data = trial_data_vv)
summary(model_vtrial)
results_table <- tidy(model_vtrial) %>%
  select(term, estimate, std.error, df, statistic, p.value)
write.csv(results_table, "results_vmpfc_capsbyrunwnowness.csv", row.names = FALSE)


model_vratingbyrun <- lmer(vv_value ~ Rating * run + CAPSB + site + age + (1 | subject),
                    data = trial_data_vv)
summary(model_vratingbyrun)
results_table <- tidy(model_vratingbyrun) %>%
  select(term, estimate, std.error, df, statistic, p.value)
write.csv(results_table, "results_vmpfc_runbynowness.csv", row.names = FALSE)

# extract predicted values for plotting ---------------------
model_data <- model.frame(model_vtrial)
model_data$Run <- trial_data_vv$Run[as.numeric(rownames(model_data))]
model_data$predicted <- predict(model_vtrial, re.form = NA)

plot_vv <- ggplot(model_data, aes(x = Run, y = vv_value, color = factor(Run))) +  # this does not work - think it through
  geom_point(alpha = 0.3) +
  geom_line(aes(y = predicted, group = Run), size = 1) +
  labs(
    title = "Trial-by-trial vmPFC with run-specific predictions",
    x = "Run",
    y = "vmPFC ROI Signal",
    color = "Run"
  ) +
  theme_minimal()
print(plot_vv)
ggsave("vv_trial_by_trial.png", plot_vv, width = 8, height = 6, dpi = 300)



# AMYGDALA VERSION WITH HEMISPHERES                                               REPORTING TRIAL AMYG.

ar_data <- data_long_runs %>%
  filter(ROI == "ar") %>%                    # keep only ar
  select(subject, site, CAPSB, age, run, roi_signal, STAR_ID) %>%  # select desired columns
  rename(
    ar_value = roi_signal,                        
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer
 ar_data <- ar_data %>%
   distinct(STAR_ID, Run, .keep_all = TRUE)

al_data <- data_long_runs %>%
  filter(ROI == "al") %>%                    # keep only al
  select(subject, site, CAPSB, age, run, roi_signal, STAR_ID) %>%  # select desired columns
  rename(
    al_value = roi_signal,                        
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer
 al_data <- al_data %>%
   distinct(STAR_ID, Run, .keep_all = TRUE)

# merge AR and AL together (many-to-many allowed)
arl_data <- ar_data %>%
  left_join(al_data, by = c("STAR_ID", "Run", "subject", "site", "CAPSB", "age"), 
  relationship = "many-to-many") %>%
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
  amyg_value ~ CAPSB * run + Rating + site + age + Hemisphere + (1 | subject),
  data = trial_data_arl_long)
summary(model_arl)
results_table <- tidy(model_arl) %>%
  select(term, estimate, std.error, df, statistic, p.value)
write.csv(results_table, "results_amygdala_capsbyrunwnowness.csv", row.names = FALSE)


model_arl_ratingbyrun <- lmer(
amyg_value ~ Rating * run + CAPSB + site + age + Hemisphere + (1 | subject),
  data = trial_data_arl_long)
summary(model_arl_ratingbyrun)
results_table <- tidy(model_arl_ratingbyrun) %>%
  select(term, estimate, std.error, df, statistic, p.value)
write.csv(results_table, "results_amygdala_runbynowness.csv", row.names = FALSE)




hr_data <- data_long_runs %>%
  filter(ROI == "hr") %>%                    # keep only hr
  select(subject, site, CAPSB, age, run, roi_signal, STAR_ID) %>%  # select desired columns
  rename(
    hr_value = roi_signal,
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer

hr_data <- hr_data %>%
  distinct(STAR_ID, Run, .keep_all = TRUE)

hl_data <- data_long_runs %>%
  filter(ROI == "hl") %>%                    # keep only hr
  select(subject, site, CAPSB, age, run, roi_signal, STAR_ID) %>%  # select desired columns
  rename(
    hl_value = roi_signal,
    STAR_ID = STAR_ID,                     # match trial data
    Run = run
  ) %>%
  mutate(Run = as.integer(Run))              # run as integer

hl_data <- hl_data %>%
  distinct(STAR_ID, Run, .keep_all = TRUE)

# Merge HR AND HL 
hrl_data <- hr_data %>%
  left_join(hl_data, by = c("STAR_ID", "Run", "subject", "site", "CAPSB", "age"), 
  relationship = "many-to-many") %>%
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
  hipp_value ~ CAPSB * run + Rating + site + age + Hemisphere + (1 | subject),
  data = trial_data_hrl_long)
summary(model_hrl_ratingbyrun)
results_table <- tidy(model_hrl) %>%
  select(term, estimate, std.error, df, statistic, p.value)
write.csv(results_table, "results_hippocampus_capsbyrunwnowness.csv", row.names = FALSE)


model_hrl_ratingbyrun <- lmer(
hipp_value ~ Rating * run + CAPSB + site + age + Hemisphere + (1 | subject),
  data = trial_data_hrl_long)
summary(model_hrl)
results_table <- tidy(model_hrl_ratingbyrun) %>%
  select(term, estimate, std.error, df, statistic, p.value)
write.csv(results_table, "results_hippocampus_runbynowness.csv", row.names = FALSE)



# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_vmpfc_runbynowness.csv")
amyg_csv <- read.csv("results_amygdala_runbynowness.csv")
hipp_csv <- read.csv("results_hippocampus_runbynowness.csv")

combined_results <- bind_rows(
  amyg_csv %>% mutate(model = "amygdala"),
  hipp_csv %>% mutate(model = "hippocampus"),
  vv_csv %>% mutate(model = "vmpfc"))

write.csv(combined_results, "results_runbynowness_combined.csv", row.names = FALSE)


# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_vmpfc_capsbyrunwnowness.csv")
amyg_csv <- read.csv("results_amygdala_capsbyrunwnowness.csv")
hipp_csv <- read.csv("results_hippocampus_capsbyrunwnowness.csv")

combined_results <- bind_rows(
  amyg_csv %>% mutate(model = "amygdala"),
  hipp_csv %>% mutate(model = "hippocampus"),
  vv_csv %>% mutate(model = "vmpfc"))
write.csv(combined_results, "results_capsbyrunwnowness_combined.csv", row.names = FALSE)



# proportion of voxels present as covariate -------------------------      FOURTH (voxels) MODEL

data_long_runs_voxels <- data_long_runs %>%
  left_join(roi_overlap %>% select(SubjectID, ROI_short, ProportionOverlap), 
            by = c("subject" = "SubjectID", "ROI" = "ROI_short"))


vmpfc_data <- data_long_runs_voxels %>%
  filter(ROI == "vv") %>%
  mutate(run = as.factor(run))
results_voxels_vmpfc <- lmer(
  roi_signal ~ CAPSB + ProportionOverlap * run + site + age + (1 | subject), # CAPSB * run before
  data = data_long_runs_voxels,
  REML = FALSE)
summary(results_voxels_vmpfc)

results_table_vmpfc <- broom.mixed::tidy(results_voxels_vmpfc) %>%
  select(term, estimate, std.error, df, statistic, p.value) %>%
  mutate(ROI = "vv")  # add ROI column
write.csv(results_table_vmpfc, "rresults_vv_voxels_poverlapbyrun.csv", row.names = FALSE)


# hippocampus version with hemispheres and voxel covariate -------------------------------      REPORTING THIS FOR HIPPOCAMPUS WITH VOXEL COVARIATE
hipp_data_voxels <- data_long_runs_voxels %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run))
results_hipp_voxels <- hipp_data_voxels %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
  map(function(df) {
    m1 <- try(
      lmer(roi_signal ~ hemi + site + age + ProportionOverlap * run + CAPSB + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE)
    return(m1)
  }) %>%
  set_names("hippocampus")
map(results_hipp_voxels, summary)

results_table <- map_dfr(results_hipp_voxels, function(model) {
  if (is.null(model)) return(NULL)  # skip failed models
  tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  
write.csv(results_table, "results_h_voxels_poverlapbyrun.csv", row.names = FALSE)

# amygdala version with hemispheres and voxel covariate -------------------------------      REPORTING THIS FOR AMYGDALA WITH VOXEL COVARIATE
amyg_data_voxels <- data_long_runs_voxels %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run))
results_amyg_voxels <- amyg_data_voxels %>%
  #mutate(run = as.numeric(run)) %>%
  group_split(ROI_group = "amygdala") %>%  # just for structure, single group
  map(function(df) {
    m1 <- try(
      lmer(roi_signal ~ hemi + site + age + ProportionOverlap * run + CAPSB + (1 | subject),
           data = df, REML = FALSE),
      silent = TRUE)
    if (inherits(m1, "try-error")) {
      cat("failed for: a\n")
      return(NULL)}
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
write.csv(results_table, "results_a_voxels_poverlapbyrun.csv", row.names = FALSE)

# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_vv_voxels_poverlapbyrun.csv")
amyg_csv <- read.csv("results_a_voxels_poverlapbyrun.csv")
hipp_csv <- read.csv("results_h_voxels_poverlapbyrun.csv")

combined_results <- bind_rows(
  amyg_csv %>% mutate(model = "amygdala"),
  hipp_csv %>% mutate(model = "hippocampus"),
  vv_csv %>% mutate(model = "vmpfc"))

write.csv(combined_results, "results_voxels_poverlapbyrun_combined.csv", row.names = FALSE)




# #------------------------------------------------------------------------------
# # PENULTIMATE, same as above INCLUDING RANDOM SLOPE FOR RUN MAKES SINGULAR
# results_m4 <- data_long_runs %>%
#   mutate(run = as.factor(run)) %>%
#   group_split(ROI) %>%
#   set_names(unique(data_long_runs$ROI)) %>%
#   map(function(data_long_runs) {
#     cat("Running ROI:", unique(data_long_runs$ROI), "\n")
#     m4 <- try(
#     #   lmer(roi_signal ~ avoid_score + CAPSB * run + site + age + (1 | subject/run),
#     #  data = data_long_runs, REML = FALSE),
#       lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + (1 | subject),
#            data = data_long_runs, REML = FALSE),
#       silent = TRUE
#     )
#     if (inherits(m4, "try-error")) {
#       cat("failed for:", unique(data_long_runs$ROI), "\n")
#       return(NULL)
#     }
#     cat("succeeded for:", unique(data_long_runs$ROI), "\n")
#     return(m4)
#   })

# map(results_m4, summary)

# results_table <- map_dfr(results_m4, function(model) {
#   if (is.null(model)) return(NULL)  # skip failed models
#   tidy(model) %>% 
#     select(term, estimate, std.error, df, statistic, p.value)
# }, .id = "ROI")  # add ROI as a column

# resultsavoidvv_table <- results_table %>%
#   dplyr::filter(ROI == "vv") %>%
#   dplyr::select(term, estimate, std.error, df, statistic, p.value)
# write.csv(resultsavoidvv_table, "results_vv_avoid.csv", row.names = FALSE)

# results_hipp <- hipp_data %>%
#   #mutate(run = as.numeric(run)) %>%
#   group_split(ROI_group = "hippocampus") %>%  # just for structure, single group
#   map(function(df) {
#     cat("Running ROI: hippocampus\n")
#     m1 <- try(
#       lmer(roi_signal ~  avoid_score + CAPSB * run + site + age + hemi +(1 | subject),
#            data = df, REML = FALSE),
#       silent = TRUE
#     )
#     if (inherits(m1, "try-error")) {
#       cat("failed for: h\n")
#       return(NULL)
#     }
#     cat("succeeded for: hippocampus\n")
#     return(m1)
#   }) %>%
#   set_names("hippocampus")
# map(results_hipp, summary)

# results_table <- map_dfr(results_hipp, function(model) {
#   if (is.null(model)) return(NULL)  # skip failed models
#   tidy(model) %>% 
#     select(term, estimate, std.error, df, statistic, p.value)
# }, .id = "ROI")  # add ROI as a column

# write.csv(results_table, "results_hipp_avoid.csv", row.names = FALSE)


# results_amygdala <- amygdala_data %>%
#   #mutate(run = as.numeric(run)) %>%
#   group_split(ROI_group = "amygdala") %>%  # just for structure, single group
#   map(function(df) {
#     cat("Running ROI: AMYG\n")
#     m1 <- try(
#       lmer(roi_signal ~ avoid_score + CAPSB * run + site + age + hemi +(1 | subject),
#            data = df, REML = FALSE),
#       silent = TRUE
#     )
#     if (inherits(m1, "try-error")) {
#       cat("failed for: h\n")
#       return(NULL)
#     }
#     cat("succeeded for: amyg\n")
#     return(m1)
#   }) %>%
#   set_names("amygdala")
# map(results_amygdala, summary)

# results_table <- map_dfr(results_amygdala, function(model) {
#   if (is.null(model)) return(NULL)  # skip failed models
#   tidy(model) %>% 
#     select(term, estimate, std.error, df, statistic, p.value)
# }, .id = "ROI")  # add ROI as a column

# write.csv(results_table, "results_amyg.csv", row.names = FALSE)









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
      silent = TRUE)
    if (inherits(m6, "try-error")) {
      cat("failed for:", unique(data_long_runs$ROI), "\n")
      return(NULL)}
    cat("succeeded for:", unique(data_long_runs$ROI), "\n")
    return(m6)})

map(results_m6, summary)

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
    cat("\nROI:", roi_name, "— model failed.\n")}})


# NOWNESS SCORE INCLUDED, INTERACTION W/CAPS BY RUN ---------------------
#     lmer(roi_signal ~  now_score * run + CAPSB + site + age + (1 | subject),

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
    cat("\nROI:", roi_name, "— model failed.\n")}})



###########################
# plots 20251008 ------------------------------------------------------------------------

# Extract subjects used in the right hippocampus model
subjects_in_model <- unique(m_hr@frame$subject)
# filter to those only
data_hr_fit <- data_hr %>% filter(subject %in% subjects_in_model)

# Add fitted values from the model
data_hr_fit <- data_hr_fit %>%
  mutate(
    fitted = predict(m_hr, newdata = data_hr_fit),  # create fitted column first
    run = as.factor(run))                          

# Plot fitted vs. avoidance
hr_fittedvsavoid <- ggplot(data_hr_fit, aes(x = avoid_score, y = fitted, color = run)) +
  geom_point(alpha = 0.6) +                         # scatter points
  geom_smooth(method = "lm", se = TRUE) +          # linear trend lines per run
  labs(
    x = "Avoidance Score",
    y = "Fitted ROI Signal (hr)",
    color = "Run",
    title = "Fitted ROI Signal vs. Avoidance Score by Run") +
  theme_minimal()
ggsave(
  filename = "hr_fitted_vs_avoidance.png",  # file name
  plot = hr_fittedvsavoid,
  width = 7,      
  height = 5,     
  dpi = 300)


# models per ROI: CAPSB ~ ROI + avg_avoid + ROI:avg_avoid + site (covariate) ---------------

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




# scatter plots --------------------------------------------------------------------

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
      color = "Run") +
    theme_minimal(base_size = 14)

  print(p)
  ggsave(
    filename = paste0("ROI_scatterplots/ROI_scatter_", roi_name, ".png"),
    plot = p,
    width = 6,
    height = 4)}

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





# number of spikes per run as covariate ----------------------------------- FIFTH (spikes) MODEL

load_spike_file <- function(filepath) {
  read_tsv(filepath, show_col_types = FALSE) %>%
    pivot_longer(cols = starts_with("Run"),
                 names_to = "run",
                 values_to = "spike_count") %>%
    mutate(
      run = str_remove(run, "Run") |> as.integer(),
      spike_count = ifelse(spike_count == "NA", NA, as.numeric(spike_count)))}

L1  <- load_spike_file("L1_spike-colrun-count.tsv") %>% mutate(site = "L")
N1  <- load_spike_file("N1_spike-colrun-count.tsv") %>% mutate(site = "N")
M1  <- load_spike_file("M1_spike-colrun-count.tsv") %>% mutate(site = "M")

spike_all <- bind_rows(L1, N1, M1)
spike_all <- spike_all %>% rename(subject = Participant)
data_long_runs <- data_long_runs %>%
  mutate(run = as.integer(run))

merged_data <- data_long_runs %>%
  left_join(spike_all, by = c("subject", "run", "site"))
print(merged_data %>% head(20))


# analysis with spike_count as covariate for amygdala ------------------------

amygdala_spikes_data <- merged_data %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run))

# Run model including spike_count per run
results_spike_amygdala <- amygdala_spikes_data %>%
  group_split(ROI_group = "amygdala") %>%  #
  map(function(df) {
    cat("Running ROI: amygdala\n")
    m4 <- try(
      lmer(
        roi_signal ~ hemi + site + age + spike_count + run * CAPSB + (1 | subject),
        data = df,
        REML = FALSE),
      silent = TRUE)
    if (inherits(m4, "try-error")) {
      cat("failed for: amygdala\n")
      return(NULL)}
    cat("succeeded for: amygdala\n")
    return(m4)
  }) %>%
  set_names("amygdala")
map(results_spike_amygdala, summary)

results_table_amygdala <- map_dfr(results_spike_amygdala, function(model) {
  if (is.null(model)) return(NULL)  
  broom.mixed::tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # Add ROI as a column
write.csv(results_table_amygdala, "results_amyg_spike_capsbyrun.csv", row.names = FALSE)

# analysis with spike_count as covariate for hippocampus -----------------------

hippocampus_spikes_data <- merged_data %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)   # keep run categorical
  )

results_spike_hippocampus <- hippocampus_spikes_data %>%
  group_split(ROI_group = "hippocampus") %>%  #
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m5 <- try(
      lmer(roi_signal ~ hemi + site + age + spike_count + run * CAPSB + (1 | subject),
        data = df,
        REML = FALSE
      ),
      silent = TRUE
    )
    if (inherits(m5, "try-error")) {
      cat("failed for: hippocampus\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m5)
  }) %>%
  set_names("hippocampus")
map(results_spike_hippocampus, summary)

results_table_hippocampus <- map_dfr(results_spike_hippocampus, function(model) {
  if (is.null(model)) return(NULL)  
  broom.mixed::tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  
write.csv(results_table_hippocampus, "results_hipp_spike_capsbyrun.csv", row.names = FALSE)

# vmpfc model with spike covariate ----------------------------------------
vmpfc_data <- merged_data %>%
  filter(ROI == "vv") %>%
  mutate(run = as.factor(run))

results_spike_vmpfc <- lmer(
  roi_signal ~ site + age + spike_count + run * CAPSB + (1 | subject),
  data = vmpfc_data,
  REML = FALSE)

summary(results_spike_vmpfc)
results_table_vmpfc <- broom.mixed::tidy(results_spike_vmpfc) %>%
  select(term, estimate, std.error, df, statistic, p.value) %>%
  mutate(ROI = "vv")  # add ROI column
write.csv(results_table_vmpfc, "results_vmpfc_spike_capsbyrun.csv", row.names = FALSE)


# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_vmpfc_spike_capsbyrun.csv")
amyg_csv <- read.csv("results_amyg_spike_capsbyrun.csv")
hipp_csv <- read.csv("results_hipp_spike_capsbyrun.csv")

combined_results <- bind_rows(
  amyg_csv %>% mutate(model = "amygdala"),
  hipp_csv %>% mutate(model = "hippocampus"),
  vv_csv %>% mutate(model = "vmpfc"))

write.csv(combined_results, "results_spike_capsbyrun_combined.csv", row.names = FALSE)




# avoidance plus spikes as covariates ----------------------------------- SIXTH (avoidance + spikes) MODEL


# avoidance mixed effects ------------------------------------- SECOND (AVOID) MODEL

vif_results3 <- results_m3 %>%
  compact() %>%  # remove NULLs if any
  map(get_vif_lmer)
vif_results3

#----------------------------------------------------------------------------
# vmpfc model with spike covariate ----------------------------------------
vmpfc_data <- merged_data %>%
  filter(ROI == "vv") %>%
  mutate(run = as.factor(run))

results_spike_vmpfc <- lmer(
  roi_signal ~ site + age + spike_count + run * avoid_score + CAPSB + (1 | subject),
  data = vmpfc_data,
  REML = FALSE)
summary(results_spike_vmpfc)
results_table_vmpfc <- broom.mixed::tidy(results_spike_vmpfc) %>%
  select(term, estimate, std.error, df, statistic, p.value) %>%
  mutate(ROI = "vv")  # add ROI column
write.csv(results_table_vmpfc, "results_vmpfc_spike_avoidbyrun.csv", row.names = FALSE)

# analysis with spike_count as covariate for amygdala ------------------------

amygdala_spikes_data <- merged_data %>%
  filter(ROI %in% c("al", "ar")) %>%
  mutate(
    hemi = ifelse(ROI == "al", "L", "R"),
    run = as.factor(run))

# Run model including spike_count per run
results_spike_amygdala <- amygdala_spikes_data %>%
  group_split(ROI_group = "amygdala") %>%  #
  map(function(df) {
    cat("Running ROI: amygdala\n")
    m4 <- try(
      lmer(
        roi_signal ~ hemi + site + age + spike_count + avoid_score * run + CAPSB + (1 | subject),
        data = df,
        REML = FALSE),
      silent = TRUE)
    if (inherits(m4, "try-error")) {
      cat("failed for: amygdala\n")
      return(NULL)}
    cat("succeeded for: amygdala\n")
    return(m4)
  }) %>%
  set_names("amygdala")
map(results_spike_amygdala, summary)

results_table_amygdala <- map_dfr(results_spike_amygdala, function(model) {
  if (is.null(model)) return(NULL)  
  broom.mixed::tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  # Add ROI as a column
write.csv(results_table_amygdala, "results_amyg_spike_avoidbyrun.csv", row.names = FALSE)

# analysis with spike_count as covariate for hippocampus -----------------------

hippocampus_spikes_data <- merged_data %>%
  filter(ROI %in% c("hl", "hr")) %>%
  mutate(
    hemi = ifelse(ROI == "hl", "L", "R"),
    run = as.factor(run)   # keep run categorical
  )

results_spike_hippocampus <- hippocampus_spikes_data %>%
  group_split(ROI_group = "hippocampus") %>%  #
  map(function(df) {
    cat("Running ROI: hippocampus\n")
    m5 <- try(
      lmer(roi_signal ~ hemi + site + age + spike_count + avoid_score * run + CAPSB + (1 | subject),
        data = df,
        REML = FALSE
      ),
      silent = TRUE
    )
    if (inherits(m5, "try-error")) {
      cat("failed for: hippocampus\n")
      return(NULL)
    }
    cat("succeeded for: hippocampus\n")
    return(m5)
  }) %>%
  set_names("hippocampus")
map(results_spike_hippocampus, summary)

results_table_hippocampus <- map_dfr(results_spike_hippocampus, function(model) {
  if (is.null(model)) return(NULL)  
  broom.mixed::tidy(model) %>% 
    select(term, estimate, std.error, df, statistic, p.value)
}, .id = "ROI")  
write.csv(results_table_hippocampus, "results_hipp_spike_avoidbyrun.csv", row.names = FALSE)

# combine all three results into one csv --------------------------------------
vv_csv   <- read.csv("results_vmpfc_spike_avoidbyrun.csv")
amyg_csv <- read.csv("results_amyg_spike_avoidbyrun.csv")
hipp_csv <- read.csv("results_hipp_spike_avoidbyrun.csv")

combined_results <- bind_rows(
  amyg_csv %>% mutate(model = "amygdala"),
  hipp_csv %>% mutate(model = "hippocampus"),
  vv_csv %>% mutate(model = "vmpfc"))

write.csv(combined_results, "results_spike_avoidbyrun_combined.csv", row.names = FALSE)
