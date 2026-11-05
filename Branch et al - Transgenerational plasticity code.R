### =========================================================
### PCA + Linear Mixed Effects Model analysis
### Transgenerational plasticity (TGP) study
###
### Response variables:
###   - PC1: composite growth axis (flower_total, Diameter, Height,
###          nodes, branches, biomass, days_to_germ)
###   - FD, sla, days_to_germ retained/analyzed as needed based on
###     PCA adequacy diagnostics below
###
### Research questions:
###   Q1: Does TGP exist (does S0_trt / grandparent treatment matter)?
###   Q2: Does it differ across regions?
###   Q3: Is it adaptive?
###   Q4: Does it evolve?
###   Q5: What is the relative contribution?
### =========================================================

library(tidyverse)
library(Hmisc)
library(lme4)
library(lmerTest)   # needed for p-values on lmer models
library(car)
library(ggfortify)
library(ggplot2)
library(visreg)
library(emmeans)
library(effectsize)
library(ggrepel)
library(stringr)
library(psych)      # KMO()
library(corrplot)   # correlation visualization
library(gt)         # nice table formatting, incl. direct .docx export
library(patchwork)  # multi-panel figure composition
library(DHARMa)     # GLMM residual diagnostics (pollen)
library(lattice)    # qqmath()

### =========================================================
### GLOBAL SETTING: Kenward-Roger throughout (not Satterthwaite)
###

emm_options(lmer.df = "kenward-roger")

### =========================================================
### PART 1: DATA CLEANING 
### =========================================================

### -----------------------------------------------------------
### 1. Read data, basic cleaning
### -----------------------------------------------------------

dat3 <- read_csv("Data/dat3.csv")
dat3$ID <- as.character(dat3$ID)

dat3$Pollen_viability <- 100 - dat3$pollen  # pollen abortion -> viability

# sample sizes per treatment combo, before any NA removal
id_counts <- dat3 %>%
  distinct(ID, S0_trt, S1_trt, Treatment, Region) %>%
  count(S0_trt, S1_trt, Treatment, Region, name = "n_ID")
id_counts

### -----------------------------------------------------------
### 2. Trait adequacy diagnostics for PCA
###    Bartlett's test: is there shared correlation structure at all?
###    KMO: item-level sampling adequacy -- flags traits that don't
###    cohere with the rest of the set (MSA < ~0.5 is problematic).
### -----------------------------------------------------------

dat3_pca <- dat3 %>%
  filter(biomass > 0, nodes == round(nodes), nodes <= 30, Diameter <= 10) %>%
  mutate(
    days_to_germ_log = log(days_to_germ),   # right-skewed (resid skew 1.4-2.1) -- log
    biomass_sqrt = sqrt(biomass)             # right-skewed (resid skew 1.5 North) -- sqrt
  )

cor_mat_full <- cor(na.omit(dat3_pca[, pca_traits_full]))
round(cor_mat_full, 2)
corrplot(cor_mat_full, method = "circle", type = "upper")

cortest.bartlett(cor_mat_full, n = nrow(na.omit(dat3_pca[, pca_traits_full])))
KMO(cor_mat_full)

# Result: overall KMO = 0.79
#   FD (0.32) and sla (0.43) -- well below the 0.5 threshold
#   days_to_germ (0.59) -- borderline, retained

pca_traits <- c("flower_total", "days_to_germ_log", "Diameter",
                "Height", "nodes", "branches", "biomass_sqrt")
cor_mat <- cor(na.omit(dat3_pca[, pca_traits]))
round(cor_mat, 2)
KMO(cor_mat)

### -----------------------------------------------------------
### 3. PCA on the 7-trait set
### -----------------------------------------------------------

pc1 <- prcomp(na.omit(dat3_pca[, pca_traits]), scale = TRUE)

summary(pc1)
pc1$rotation

# extract scores and loadings for plotting
scores <- as.data.frame(pc1$x)
loadings <- as.data.frame(pc1$rotation)
loadings$Trait <- rownames(loadings)

loadings$Trait[loadings$Trait == "flower_total"] <- "Flower Number"
loadings$Trait[loadings$Trait == "days_to_germ"] <- "Germination Date"

loadings_scaled <- loadings
loadings_scaled[, 1:2] <- loadings_scaled[, 1:2] * 5

PCA_plot <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.1, color = "skyblue3") +
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = loadings_scaled, aes(x = PC1, y = PC2, label = Trait),
                  size = 4, color = "black") +
  xlab("PC1 (61.1%) - Growth") +
  ylab("PC2 (14.0%) - Germination Date") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(face = "bold"))

PCA_plot   # shown for transparency; PC2 not used as a formal response

var_contrib_pct <- sweep(as.matrix(loadings[,1:7]^2), 2, colSums(loadings[,1:7]^2), "/") * 100
round(var_contrib_pct[, 1:2], 1)

### -----------------------------------------------------------
### 4. Rebuild full dataset with PC1 
### -----------------------------------------------------------

pc_data <- as.data.frame(pc1$x) %>% dplyr::select(PC1)  # keep PC1 only

dat4 <- na.omit(dat3_pca[, c("ID", "Region", "Site", "S0_ID", "S1_ID", "TRTHIS",
                             "S1_trt", "S0_trt", "Treatment", "PrePeak", "Rep",
                             pca_traits)])
dat <- cbind(dat4, pc_data)

dat_subset <- dat %>%
  dplyr::select(ID, Treatment, PC1)

dat <- dat3 %>%
  dplyr::left_join(dat_subset, by = c("ID", "Treatment"))

dat$ID <- as.character(dat$ID)
dat$S0_ID <- as.character(dat$S0_ID)

dat$Rep       <- as.factor(dat$Rep)
dat$S0_trt    <- as.factor(dat$S0_trt)
dat$S1_trt    <- as.factor(dat$S1_trt)
dat$PrePeak   <- as.factor(dat$PrePeak)
dat$Region    <- as.factor(dat$Region)
dat$Treatment <- as.factor(dat$Treatment)

### -----------------------------------------------------------
### 5. Subset to pre-drought (baseline TGP question) 
### -----------------------------------------------------------

dat5 <- filter(dat, PrePeak == "Pre")

## Set sum-to-zero contrasts

contrasts(dat5$Rep)       <- contr.sum
contrasts(dat5$S0_trt)    <- contr.sum
contrasts(dat5$S1_trt)    <- contr.sum
contrasts(dat5$Treatment) <- contr.sum
contrasts(dat5$Region)    <- contr.sum

### =========================================================
### PART 2: 
###
### Q1: Does TGP exist?
### Q2: Des it differ between regions?
### =========================================================

### -----------------------------------------------------------
### 1. Model-fitting pipeline (reusable for PC1, FD, sla,
###    days_to_germ): 
###    Justify additive-S1_trt structure via LRT
### -----------------------------------------------------------
                     
fit_trait_models <- function(response, data) {
  f_add  <- as.formula(paste(response, "~ Region * Treatment * S0_trt + S1_trt + Rep + (1 | S0_ID)"))
  f_full <- as.formula(paste(response, "~ Region * Treatment * S0_trt * S1_trt + Rep + (1 | S0_ID)"))
  
  m_add_ml  <- lmer(f_add,  data = data, REML = FALSE)
  m_full_ml <- lmer(f_full, data = data, REML = FALSE)
  
  lrt <- anova(m_add_ml, m_full_ml)
  
  m_final <- lmer(f_add, data = data)  # REML = TRUE default
  
  list(
    lrt = lrt,
    singular = isSingular(m_final),
    model = m_final,
    summary = summary(m_final),
    anova_typeIII = anova(m_final, ddf = "Kenward-Roger")
  )
}

## --- PC1 (growth composite) ---
res_PC1 <- fit_trait_models("PC1", dat5)
res_PC1$lrt
res_PC1$singular
res_PC1$summary
res_PC1$anova_typeIII

## --- days_to_germ ---
res_germ <- fit_trait_models("log(days_to_germ)", dat5)
res_germ$lrt
res_germ$singular
res_germ$summary
res_germ$anova_typeIII

## --- FD (flowering date) ---
res_FD <- fit_trait_models("FD", dat5)
res_FD$lrt
res_FD$singular
res_FD$summary
res_FD$anova_typeIII

## --- sla ---
# trait-specific filter: exclude the 4 rows with sla >= 700 
dat5_sla <- filter(dat5, sla < 700)
res_sla <- fit_trait_models("sla", dat5_sla)
res_sla$lrt
res_sla$singular
res_sla$summary
res_sla$anova_typeIII


### -----------------------------------------------------------
### 2. Result visualization
### -----------------------------------------------------------

## PC1 (growth) ---
emm_PC1 <- as.data.frame(emmeans(res_PC1$model, ~ S0_trt))

pA <- ggplot() +
  geom_point(data = emm_PC1, aes(x = S0_trt, y = emmean),
             size = 5, color = "black") +
  geom_errorbar(data = emm_PC1, aes(x = S0_trt, ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "PC1 (growth)", title = "") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold")
  )

## Days to germination ---

emm_germ <- as.data.frame(emmeans(res_germ$model, ~ S0_trt, type = "response"))

pB <- ggplot() +
  geom_point(data = emm_germ, aes(x = S0_trt, y = response),
             size = 3, color = "black") +
  geom_errorbar(data = emm_germ, aes(x = S0_trt, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, color = "black") +
  labs(x = "Grandparent treatment", y = "Days to germination", title = "") +
  theme_classic(base_size = 13)

## Flowering date ---
emm_FD <- as.data.frame(emmeans(res_FD$model, ~ S0_trt))

pC <- ggplot() +
  geom_point(data = emm_FD, aes(x = S0_trt, y = emmean),
             size = 3, color = "black") +
  geom_errorbar(data = emm_FD, aes(x = S0_trt, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, color = "black") +
  labs(x = "Grandparent treatment", y = "Flowering date", title = "") +
  theme_classic(base_size = 13)

## SLA, by Region (interaction) ---
# Untransformed (raw sla) -- see res_sla definition for rationale.
emm_sla <- as.data.frame(emmeans(res_sla$model, ~ S0_trt | Region))

pD <- ggplot(emm_sla, aes(x = S0_trt, y = emmean, color = Region, group = Region)) +
  geom_point(size = 5) +
  geom_line(aes(linetype = Region), linewidth = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  scale_linetype_manual(breaks = c("North", "South"), values = c("solid", "dashed")) +
  scale_color_manual(breaks = c("North", "South"), values = c("navyblue", "magenta2")) +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Specific Leaf Area", title = "",
       color = "Region", linetype = "Region") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold")
  )

#Save plots
#ggsave("PC1_S0.png", pA, width = 5, height = 4.5, dpi = 300)
#ggsave("SLA_S0_Region.png", pD, width = 5, height = 4.5, dpi = 300)


## --- Combine into one figure ---
combined_fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

combined_fig

# save for manuscript
#ggsave("results_panel_figure.png", combined_fig, width = 9, height = 8, dpi = 300)

## --- SLA: simple effects (D-W within each Region) and interaction
simple_effects_sla <- emmeans(res_sla$model, ~ S0_trt | Region) %>%
  pairs()
simple_effects_sla

emm_cross_sla <- emmeans(res_sla$model, ~ S0_trt * Region)
interaction_contrast_sla <- contrast(emm_cross_sla, interaction = "pairwise")
interaction_contrast_sla

### -----------------------------------------------------------
### 3. Supplemental figures: all significant terms NOT already
### -----------------------------------------------------------

#main effect
plot_main_effect <- function(model, data, response_var, factor_var,
                             y_lab, x_lab, title_lab, response_scale = FALSE) {
  emm_obj <- if (response_scale) {
    emmeans(model, as.formula(paste("~", factor_var)), type = "response")
  } else {
    emmeans(model, as.formula(paste("~", factor_var)))
  }
  emm <- as.data.frame(emm_obj)
  
  y_col   <- if ("response" %in% names(emm)) "response" else "emmean"
  lcl_col <- if ("asymp.LCL" %in% names(emm)) "asymp.LCL" else "lower.CL"
  ucl_col <- if ("asymp.UCL" %in% names(emm)) "asymp.UCL" else "upper.CL"
  
  ggplot() +
    geom_point(data = emm, aes(x = .data[[factor_var]], y = .data[[y_col]]),
               size = 3, color = "black") +
    geom_errorbar(data = emm, aes(x = .data[[factor_var]], ymin = .data[[lcl_col]], ymax = .data[[ucl_col]]),
                  width = 0.1, color = "black") +
    labs(x = x_lab, y = y_lab, title = title_lab) +
    theme_classic(base_size = 13)
}

## two-way interaction plot (e.g. Region:Treatment) ----
plot_interaction <- function(model, data, response_var, x_var, group_var,
                             y_lab, x_lab, group_lab, title_lab,
                             group_colors = c("navyblue", "magenta2"),
                             group_linetypes = c("solid", "dotted")) {
  emm <- as.data.frame(emmeans(model, as.formula(paste("~", x_var, "|", group_var))))
  
  ggplot() +
    geom_point(data = emm, aes(x = .data[[x_var]], y = emmean, color = .data[[group_var]], group = .data[[group_var]]),
               position = position_dodge(0.5), size = 3) +
    geom_line(data = emm, aes(x = .data[[x_var]], y = emmean, color = .data[[group_var]],
                              group = .data[[group_var]], linetype = .data[[group_var]]),
              position = position_dodge(0.5), linewidth = 0.8) +
    geom_errorbar(data = emm, aes(x = .data[[x_var]], ymin = lower.CL, ymax = upper.CL,
                                  color = .data[[group_var]], group = .data[[group_var]]),
                  width = 0.1, position = position_dodge(0.5)) +
    scale_color_manual(values = group_colors) +
    scale_linetype_manual(values = group_linetypes) +
    labs(x = x_lab, y = y_lab, title = title_lab, color = group_lab, linetype = group_lab) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
}

## ============================================================
## Region main effects
## ============================================================
pS1 <- plot_main_effect(res_PC1$model, dat5, "PC1", "Region",
                        y_lab = "PC1 (growth)", x_lab = "Region",
                        title_lab = "")

pS2 <- plot_main_effect(res_germ$model, dat5, "days_to_germ", "Region",
                        y_lab = "Days to germination", x_lab = "Region",
                        title_lab = "", response_scale = TRUE)

## ============================================================
## Treatment main effects
## ============================================================
pS3 <- plot_main_effect(res_PC1$model, dat5, "PC1", "Treatment",
                        y_lab = "PC1 (growth)", x_lab = "Current treatment",
                        title_lab = "")

pS4 <- plot_main_effect(res_FD$model, dat5, "FD", "Treatment",
                        y_lab = "Flowering date", x_lab = "Current treatment",
                        title_lab = "")

pS5 <- plot_main_effect(res_sla$model, dat5, "sla", "Treatment",
                        y_lab = "Specific Leaf Area", x_lab = "Current treatment",
                        title_lab = "")

## ============================================================
## Region:Treatment interactions
## ============================================================
pS6 <- plot_interaction(res_PC1$model, dat5, "PC1", "Treatment", "Region",
                        y_lab = "PC1 (growth)", x_lab = "Current treatment",
                        group_lab = "Region",
                        title_lab = "")

pS7 <- plot_interaction(res_FD$model, dat5, "FD", "Treatment", "Region",
                        y_lab = "Flowering date", x_lab = "Current treatment",
                        group_lab = "Region",
                        title_lab = "")

## ============================================================
## S1_trt (parent treatment) main effects
## ============================================================
pS8 <- plot_main_effect(res_germ$model, dat5, "days_to_germ", "S1_trt",
                        y_lab = "Days to germination", x_lab = "Parent treatment (S1_trt)",
                        title_lab = "", response_scale = TRUE)

pS9 <- plot_main_effect(res_FD$model, dat5, "FD", "S1_trt",
                        y_lab = "Flowering date", x_lab = "Parent treatment (S1_trt)",
                        title_lab = "")

## ============================================================
## Combine into supplemental figure pages
## ============================================================

supp_fig_main_effects <- (pB|pC) / (pS1 | pS2) / (pS3 | pS4) / (pS5 | pS8) / (pS9 | plot_spacer()) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F","G","H","I")))

supp_fig_interactions <- (pS6 | pS7) +
  plot_annotation(tag_levels = list(c("S6", "S7")))

supp_fig_S1trt <- (pS8 | pS9) +
  plot_annotation(tag_levels = list(c("S8", "S9")))

supp_fig_main_effects
supp_fig_interactions
supp_fig_S1trt

#ggsave("supp_fig_main_effects.png", supp_fig_main_effects, width = 9, height = 11, dpi = 300)
#ggsave("supp_fig_interactions.png", supp_fig_interactions, width = 9, height = 4.5, dpi = 300)
#ggsave("supp_fig_S1trt.png", supp_fig_S1trt, width = 9, height = 4.5, dpi = 300)

### =========================================================
### Part 3 - Pollen
### Q3. Is it adaptive? 
### Q4. Does it evolve?
### =========================================================

### -----------------------------------------------------------
### 3A. Pollen viability (binomial GLMM)
###   Response: cbind(Viable, Total - Viable)
###   Baseline TGP / regional-difference question
###    (mirrors dat5 for PC1/germ/FD/sla)
### -----------------------------------------------------------

### -----------------------------------------------------------
### 1. Read, merge, and prepare data
### -----------------------------------------------------------

pollen <- read_csv("Data/S2_pollen.csv")
s0 <- read.csv("Data/s2_s0_ids.csv")
Treatment_history <- read_csv("Data/Trthis_IDs.csv")

dat2 <- merge(s0, Treatment_history, by = "S0_ID")
pollen <- merge(pollen, dat2, by = "ID")
pollen$Total <- pollen$Viable + pollen$Aborted

# convert factors ONCE, before splitting -- avoids duplicating this

pollen$Treatment <- as.factor(pollen$Treatment)
pollen$S0_trt    <- as.factor(pollen$S0_trt)
pollen$S1_trt    <- as.factor(pollen$S1_trt)
pollen$PrePeak   <- as.factor(pollen$PrePeak)
pollen$Region    <- as.factor(pollen$Region)
pollen$Obs       <- 1:nrow(pollen)   # unique ID per row, for OLRE

# filter to pre-drought
pollen_pre <- filter(pollen, PrePeak == "Pre")

contrasts(pollen_pre$Region)    <- contr.sum
contrasts(pollen_pre$S0_trt)    <- contr.sum
contrasts(pollen_pre$S1_trt)    <- contr.sum
contrasts(pollen_pre$Treatment) <- contr.sum

## --- Fit pipeline: additive-S1_trt justified via LRT, OLRE for
##     overdispersion, matching the region-split pipeline below ---

f_add_pre  <- cbind(Viable, Total - Viable) ~ Region * Treatment * S0_trt + S1_trt + (1 | S0_ID) + (1 | Obs)
f_full_pre <- cbind(Viable, Total - Viable) ~ Region * Treatment * S0_trt * S1_trt + (1 | S0_ID) + (1 | Obs)

m_add_pre  <- glmer(f_add_pre,  family = binomial, data = pollen_pre)
m_full_pre <- glmer(f_full_pre, family = binomial, data = pollen_pre)

lrt_pre <- anova(m_add_pre, m_full_pre)
lrt_pre   # justifies (or not) additive S1_trt for this baseline model

sim_res_pre <- simulateResiduals(m_add_pre, n = 1000)
testDispersion(sim_res_pre)
plot(sim_res_pre)

## --- Fixed-effects tests ---
Anova(m_add_pre, type = "III")

## --- Visualization: Region x S0_trt interaction (mirrors the SLA
emm_pollen_pre <- as.data.frame(emmeans(m_add_pre, ~ S0_trt | Region, type = "response"))

p_pollen_baseline <- ggplot(emm_pollen_pre, aes(x = S0_trt, y = prob, color = Region, group = Region)) +
  geom_point(size = 5) +
  geom_line(aes(linetype = Region), linewidth = 1) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  scale_color_manual(breaks = c("North", "South"), values = c("navyblue", "magenta2")) +
  scale_linetype_manual(breaks = c("North", "South"), values = c("solid", "dashed")) +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Pollen viability (proportion)",
       color = "Region", linetype = "Region") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold")
  )

p_pollen_baseline
#ggsave("pollen_baseline_RegionxS0trt.png", p_pollen_baseline, width = 5.5, height = 4.5, dpi = 300)


### =========================================================
### 2B: Region-split PrePeak follow-up
### =========================================================

pollen_n <- filter(pollen, Region == "North")
pollen_s <- filter(pollen, Region == "South")

### -----------------------------------------------------------
### 1. Set sum-to-zero contrasts (both regions)
### -----------------------------------------------------------

contrasts(pollen_n$S1_trt)    <- contr.sum
contrasts(pollen_n$PrePeak)   <- contr.sum
contrasts(pollen_n$S0_trt)    <- contr.sum
contrasts(pollen_n$Treatment) <- contr.sum

contrasts(pollen_s$S1_trt)    <- contr.sum
contrasts(pollen_s$PrePeak)   <- contr.sum
contrasts(pollen_s$S0_trt)    <- contr.sum
contrasts(pollen_s$Treatment) <- contr.sum


### -----------------------------------------------------------
### 2. Fit pipeline
### -----------------------------------------------------------

fit_pollen_model <- function(data, region_label) {
  
  f_add  <- cbind(Viable, Total - Viable) ~ PrePeak * S0_trt * Treatment + S1_trt + (1 | S0_ID)
  f_full <- cbind(Viable, Total - Viable) ~ PrePeak * S0_trt * Treatment * S1_trt + (1 | S0_ID)

  f_add_olre  <- update(f_add,  . ~ . + (1 | Obs))
  f_full_olre <- update(f_full, . ~ . + (1 | Obs))
  
  m_add_olre  <- glmer(f_add_olre,  family = binomial, data = data)
  m_full_olre <- glmer(f_full_olre, family = binomial, data = data)
  
  cat("\n---", region_label, ": isSingular (add / full) ---\n")
  print(c(add = isSingular(m_add_olre), full = isSingular(m_full_olre)))
  
  # LRT: does allowing S1_trt to interact improve fit?
  lrt <- anova(m_add_olre, m_full_olre)
  
  sim_res <- simulateResiduals(m_add_olre, n = 1000)
  
  list(
    m_add = m_add_olre,
    m_full = m_full_olre,
    lrt = lrt,
    sim_res = sim_res,
    dispersion_test = testDispersion(sim_res, plot = FALSE)
  )
}

## --- Inspect diagnostics before trusting results ---
res_pollen_N$lrt
res_pollen_N$dispersion_test          # DHARMa test; dispersion ~1, p > .05 is good
plot(res_pollen_N$sim_res)            # full DHARMa residual diagnostics

res_pollen_S$lrt
res_pollen_S$dispersion_test
plot(res_pollen_S$sim_res)

## --- Fixed-effects tests: Wald chi-square (primary for GLMMs) ---
Anova(res_pollen_N$m_add, type = "III")
Anova(res_pollen_S$m_add, type = "III")



## --- Pollen, South: PrePeak:S0_trt
simple_effects_pollen_S <- emmeans(res_pollen_S$m_add, ~ S0_trt | PrePeak, type = "response") %>%
  pairs()
simple_effects_pollen_S

emm_cross_pollen_S <- emmeans(res_pollen_S$m_add, ~ S0_trt * PrePeak, type = "response")
interaction_contrast_pollen_S <- contrast(emm_cross_pollen_S, interaction = "pairwise")
interaction_contrast_pollen_S





simple_effects_pollen_baseline <- emmeans(m_add_pre, ~ S0_trt | Region, type = "response") %>%
  pairs()
simple_effects_pollen_baseline

emm_cross_pollen_baseline <- emmeans(m_add_pre, ~ S0_trt * Region, type = "response")
interaction_contrast_pollen_baseline <- contrast(emm_cross_pollen_baseline, interaction = "pairwise")
interaction_contrast_pollen_baseline

## --- Pollen: percentage difference (D vs W), per region --

emm_pollen_means <- as.data.frame(emmeans(m_add_pre, ~ S0_trt | Region, type = "response"))
emm_pollen_means

pollen_pct_diff <- emm_pollen_means %>%
  dplyr::select(Region, S0_trt, prob) %>%
  pivot_wider(names_from = S0_trt, values_from = prob) %>%
  mutate(
    pct_point_diff = (D - W) * 100,              # percentage-POINT difference (D minus W)
    pct_relative_change = (D - W) / W * 100       # relative % change, D vs W as baseline
  )
pollen_pct_diff
### -----------------------------------------------------------
### 3. Pollen viability figures
### -----------------------------------------------------------

## --- North: S0_trt main effect ---
emm_pollen_N <- as.data.frame(emmeans(res_pollen_N$m_add, ~ S0_trt, type = "response"))

pN <- ggplot() +
  geom_point(data = emm_pollen_N, aes(x = S0_trt, y = prob),
             size = 5, color = "navyblue") +
  geom_errorbar(data = emm_pollen_N, aes(x = S0_trt, ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, color = "navyblue") +
  labs(x = "Grandparent treatment", y = "Pollen viability (proportion)",
       title = "North") +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

## --- South: S0_trt main effect ---
emm_pollen_S <- as.data.frame(emmeans(res_pollen_S$m_add, ~ S0_trt, type = "response"))

pS <- ggplot() +
  geom_point(data = emm_pollen_S, aes(x = S0_trt, y = prob),
             size = 5, color = "magenta2") +
  geom_errorbar(data = emm_pollen_S, aes(x = S0_trt, ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, color = "magenta2") +
  labs(x = "Grandparent treatment", y = "Pollen viability (proportion)",
       title = "South") +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

## --- South: PrePeak x S0_trt interaction ---
emm_pollen_S <- as.data.frame(emmeans(res_pollen_S$m_add, ~ S0_trt | PrePeak, type = "response"))

pS <- ggplot(emm_pollen_S, aes(x = S0_trt, y = prob, color = PrePeak, group = PrePeak)) +
  geom_point(size = 5) +
  geom_line(aes(linetype = PrePeak),linewidth = 1) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  scale_color_manual(values = c("darkorange", "skyblue")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Pollen viability (proportion)",
       title = "South",
       color = "Drought stage", linetype = "Drought stage") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

p_FD_North_interaction
#ggsave("FD_North_YearxS0_interaction.png", p_FD_North_interaction, width = 5.5, height = 4.5, dpi = 300)


pN
pS

#ggsave("pollen_North.png", pN, width = 5, height = 4.5, dpi = 300)
#ggsave("pollen_South.png", pS, width = 5, height = 4.5, dpi = 300)

### =========================================================
### PART 3: Max seed set (lmer, despite being count data)
### Q3. Is it adaptive?
### Q4. Does it evolve?
### =========================================================

### -----------------------------------------------------------
### PART 3a: Baseline (Pre-only, unified Region model)
### -----------------------------------------------------------

### -----------------------------------------------------------
### 1. Read, merge, and prepare data
### -----------------------------------------------------------

seed <- read.csv("Data/S0_S1_S2_seed.csv")
seed$Rep <- as.numeric(!(seed$ID %% 2)) + 1
seed$Seed_total <- as.numeric(seed$Seed_total)
seed$ID <- as.character(seed$ID)

ids <- read.csv("Data/Trthis_IDs.csv")
ids$S1_ID <- as.character(ids$S1_ID)
ids$ID <- as.character(ids$ID)

new <- inner_join(seed, ids, by = "ID")

id2 <- read.csv("Data/S2_s0_ids.csv")
seed <- inner_join(new, id2, by = "S0_ID")

seed$ID    <- as.character(seed$ID)
seed$S0_ID <- as.character(seed$S0_ID)
seed$S1_ID <- as.character(seed$S1_ID)

seed$Rep       <- as.factor(seed$Rep)
seed$S0_trt    <- as.factor(seed$S0_trt)
seed$S1_trt    <- as.factor(seed$S1_trt)
seed$PrePeak   <- as.factor(seed$PrePeak)
seed$Region    <- as.factor(seed$Region)
seed$Treatment <- as.factor(seed$Treatment)

#filter to pre drought
seed_pre <- filter(seed, PrePeak == "Pre")

n_distinct(seed_pre$S0_ID)                                  # expect 20
seed_pre %>% distinct(S0_ID, Region) %>% count(Region)       # expect 10/10
seed_pre %>% distinct(S0_ID, Region, S0_trt) %>% count(Region, S0_trt)

contrasts(seed_pre$S0_trt)    <- contr.sum
contrasts(seed_pre$S1_trt)    <- contr.sum
contrasts(seed_pre$Region)    <- contr.sum
contrasts(seed_pre$Treatment) <- contr.sum
contrasts(seed_pre$Rep)       <- contr.sum


res_seed <- fit_trait_models("Seed_total", seed_pre)

# basic residual diagnostics 
qqmath(res_seed$model)
hist(resid(res_seed$model), breaks = 40)

res_seed$lrt              # justifies (or not) additive S1_trt
res_seed$singular         # should be FALSE
res_seed$summary
res_seed$anova_typeIII


# robustness check (same as PC1/germ/FD/sla)
Anova(res_seed$model, type = "III", test.statistic = "F")   # Kenward-Roger



## --- Seed set, baseline: Region:S0_trt (Wald p=.0003, LRT p=.0105) ---
simple_effects_seed <- emmeans(res_seed$model, ~ S0_trt | Region) %>%
  pairs()
simple_effects_seed

emm_cross_seed <- emmeans(res_seed$model, ~ S0_trt * Region)
interaction_contrast_seed <- contrast(emm_cross_seed, interaction = "pairwise")
interaction_contrast_seed

## --- Seed set: percentage difference (D vs W), per region ---

emm_seed_means <- as.data.frame(emmeans(res_seed$model, ~ S0_trt | Region))
emm_seed_means

seed_pct_diff <- emm_seed_means %>%
  dplyr::select(Region, S0_trt, emmean) %>%
  pivot_wider(names_from = S0_trt, values_from = emmean) %>%
  mutate(
    raw_diff = D - W,
    pct_relative_change = (D - W) / W * 100
  )
seed_pct_diff

## --- Visualization: Region x S0_trt interaction (Region:S0_trt

emm_seed <- as.data.frame(emmeans(res_seed$model, ~ S0_trt | Region))

p_seed_baseline <- ggplot(emm_seed, aes(x = S0_trt, y = emmean, color = Region, group = Region)) +
  geom_point(size = 5) +
  geom_line(aes(linetype = Region), linewidth = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  scale_color_manual(breaks = c("North", "South"), values = c("navyblue", "magenta2")) +
  scale_linetype_manual(breaks = c("North", "South"), values = c("solid", "dashed")) +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Max seed set",
       color = "Region", linetype = "Region") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold")
  )

p_seed_baseline
#ggsave("seed_baseline_RegionxS0trt.png", p_seed_baseline, width = 5.5, height = 4.5, dpi = 300)

### -----------------------------------------------------------
### PART 3B: Region-split PrePeak follow-up
### -----------------------------------------------------------

seed_n <- filter(seed, Region == "North")
seed_s <- filter(seed, Region == "South")

contrasts(seed_n$S1_trt)    <- contr.sum
contrasts(seed_n$PrePeak)   <- contr.sum
contrasts(seed_n$S0_trt)    <- contr.sum
contrasts(seed_n$Treatment) <- contr.sum
contrasts(seed_n$Rep)       <- contr.sum

contrasts(seed_s$S1_trt)    <- contr.sum
contrasts(seed_s$PrePeak)   <- contr.sum
contrasts(seed_s$S0_trt)    <- contr.sum
contrasts(seed_s$Treatment) <- contr.sum
contrasts(seed_s$Rep)       <- contr.sum

fit_seed_region_model <- function(data, region_label) {
  f_add  <- as.formula("Seed_total ~ PrePeak * S0_trt * Treatment + S1_trt + Rep + (1 | S0_ID)")
  f_full <- as.formula("Seed_total ~ PrePeak * S0_trt * Treatment * S1_trt + Rep + (1 | S0_ID)")
  
  m_add_ml  <- lmer(f_add,  data = data, REML = FALSE)
  m_full_ml <- lmer(f_full, data = data, REML = FALSE)
  
  lrt <- anova(m_add_ml, m_full_ml)
  m_final <- lmer(f_add, data = data)   # REML = TRUE default
  
  cat("\n---", region_label, ": isSingular ---\n")
  print(isSingular(m_final))
  
  list(
    lrt = lrt,
    singular = isSingular(m_final),
    model = m_final,
    summary = summary(m_final),
    anova_typeIII = anova(m_final, ddf = "Kenward-Roger")
  )
}

res_seed_N <- fit_seed_region_model(seed_n, "North")
res_seed_S <- fit_seed_region_model(seed_s, "South")

res_seed_N$lrt
res_seed_N$singular
res_seed_N$anova_typeIII

res_seed_S$lrt
res_seed_S$singular
res_seed_S$anova_typeIII

# robustness check
Anova(res_seed_N$model, type = "III", test.statistic = "F")
Anova(res_seed_S$model, type = "III", test.statistic = "F")

### -----------------------------------------------------------
### Visualizations
### -----------------------------------------------------------

emm_seed_N <- as.data.frame(emmeans(res_seed_N$model, ~ S0_trt))
emm_seed_S <- as.data.frame(emmeans(res_seed_S$model, ~ S0_trt))

pSeedN <- ggplot() +
  geom_point(data = emm_seed_N, aes(x = S0_trt, y = emmean), size = 5, color = "navyblue") +
  geom_errorbar(data = emm_seed_N, aes(x = S0_trt, ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "navyblue") +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Max seed set",
       title = "North") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )


pSeedS <- ggplot() +
  geom_point(data = emm_seed_S, aes(x = S0_trt, y = emmean), size = 5, color = "magenta2") +
  geom_errorbar(data = emm_seed_S, aes(x = S0_trt, ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "magenta2") +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Max seed set",
       title = "South") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )


seed_NS_fig <- (pSeedN | pSeedS) +
  plot_annotation(tag_levels = list(c("C", "D")))

seed_NS_fig
#ggsave("seed_North_South.png", seed_NS_fig, width = 9, height = 4.5, dpi = 300)


pollen_NS_fig <- (pN | pS) +
  plot_annotation(tag_levels = list(c("A", "B")))

pollen_NS_fig
#ggsave("pollen_North_South.png", pollen_NS_fig, width = 9, height = 4.5, dpi = 300)





### -----------------------------------------------------------
### 3. Supplemental figures
### -----------------------------------------------------------

emm_seed_trt <- as.data.frame(emmeans(res_seed$model, ~ Treatment))

pSeedTrt <- ggplot() +
  geom_point(data = emm_seed_trt, aes(x = Treatment, y = emmean), size = 3, color = "black") +
  geom_errorbar(data = emm_seed_trt, aes(x = Treatment, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, color = "black") +
  labs(x = "Current treatment", y = "Max seed set",
       title = "S10. Current treatment effect on seed set") +
  theme_classic(base_size = 13)

pSeedTrt
#ggsave("supp_seed_treatment.png", pSeedTrt, width = 5, height = 4.5, dpi = 300)



### =========================================================
### PART 4: Region-split PrePeak follow-up for growth/phenology
###    traits (PC1, Germination, Flowering Date, SLA)
### Q4 - Does it evolve?
### =========================================================
fit_growth_region_model <- function(formula_resp, data, region_label, trait_label) {
  f_add  <- as.formula(paste(formula_resp, "~ PrePeak * S0_trt * Treatment + S1_trt + Rep + (1 | S0_ID)"))
  f_full <- as.formula(paste(formula_resp, "~ PrePeak * S0_trt * Treatment * S1_trt + Rep + (1 | S0_ID)"))
  
  m_add_ml  <- lmer(f_add,  data = data, REML = FALSE)
  m_full_ml <- lmer(f_full, data = data, REML = FALSE)
  
  lrt <- anova(m_add_ml, m_full_ml)
  m_final <- lmer(f_add, data = data)  # REML = TRUE default
  
  cat("\n---", trait_label, region_label, ": isSingular ---\n")
  print(isSingular(m_final))
  
  list(
    lrt = lrt,
    singular = isSingular(m_final),
    model = m_final,
    anova_typeIII = anova(m_final, ddf = "Kenward-Roger")
  )
}

## --- PC1 ---
res4_PC1_N <- fit_growth_region_model("PC1", dat_n, "North", "PC1")
res4_PC1_S <- fit_growth_region_model("PC1", dat_s, "South", "PC1")

## --- Germination (log-transformed, per established decision) ---
res4_germ_N <- fit_growth_region_model("log(days_to_germ)", dat_n, "North", "Germination")
res4_germ_S <- fit_growth_region_model("log(days_to_germ)", dat_s, "South", "Germination")

## --- Flowering Date (untransformed, per established decision) ---
res4_FD_N <- fit_growth_region_model("FD", dat_n, "North", "Flowering Date")
res4_FD_S <- fit_growth_region_model("FD", dat_s, "South", "Flowering Date")

dat_n_sla <- filter(dat_n, sla < 700)
dat_s_sla <- filter(dat_s, sla < 700)

res4_sla_N <- fit_growth_region_model("sla", dat_n_sla, "North", "SLA")
res4_sla_S <- fit_growth_region_model("sla", dat_s_sla, "South", "SLA")


f_sla_N_interactive <- sla ~ PrePeak * S0_trt * Treatment * S1_trt + Rep + (1 | S0_ID)
m_sla_N_interactive <- lmer(f_sla_N_interactive, data = dat_n_sla)   # REML = TRUE
res4_sla_N_interactive <- list(
  model = m_sla_N_interactive,
  singular = isSingular(m_sla_N_interactive),
  anova_typeIII = anova(m_sla_N_interactive, ddf = "Kenward-Roger")
)
res4_sla_N_interactive$singular   # confirm FALSE
res4_sla_N_interactive$anova_typeIII

## --- Inspect diagnostics before trusting any result ---
res4_PC1_N$lrt;  res4_PC1_S$lrt
res4_germ_N$lrt; res4_germ_S$lrt
res4_FD_N$lrt;   res4_FD_S$lrt
res4_sla_N$lrt;  res4_sla_S$lrt

res4_PC1_N$anova_typeIII;  res4_PC1_S$anova_typeIII
res4_germ_N$anova_typeIII; res4_germ_S$anova_typeIII
res4_FD_N$anova_typeIII;   res4_FD_S$anova_typeIII
res4_sla_N$anova_typeIII;  res4_sla_S$anova_typeIII


## Visualization Flowering date, North -- Year (PrePeak) x Grandparent ---

emm_FD_N_plot <- as.data.frame(emmeans(res4_FD_N$model, ~ S0_trt | PrePeak))

p_FD_North_interaction <- ggplot(emm_FD_N_plot, aes(x = S0_trt, y = emmean, color = PrePeak, group = PrePeak)) +
  geom_point(size = 5) +
  geom_line(aes(linetype = PrePeak), linewidth = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  scale_color_manual(values = c("darkorange", "skyblue")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_discrete(labels = c("D" = "Dry", "W" = "Wet")) + # Explicitly renames the tick labels
  labs(x = "Grandparent treatment", y = "Flowering Date",
       title = "North",
       color = "Drought stage", linetype = "Drought stage") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16) # Changes both X and Y axis tick text size
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

p_FD_North_interaction
#ggsave("FD_North_YearxS0_interaction.png", p_FD_North_interaction, width = 5, height = 4.5, dpi = 300)


## =========================================================
### PART 5: Effect-size summary across all traits, by region
###    Main effects only (PrePeak/Year, S0_trt/Grandparent,
###    S1_trt/Parent, Treatment/Current) -- interactions excluded
### Q5 - Relative effect of each. 
### =========================================================


## -----------------------------------------------------------
## 1. Build region-split, full-PrePeak datasets for the growth/
##    phenology traits (PC1, germination, FD, SLA) -- these were
##    previously only modeled on the Pre-only unified dataset (dat5);
##    this is the first time PrePeak is included for them.
## -----------------------------------------------------------

dat_n <- filter(dat, Region == "North")
dat_s <- filter(dat, Region == "South")

for (d in c("dat_n", "dat_s")) {
  df <- get(d)
  contrasts(df$PrePeak)   <- contr.sum
  contrasts(df$S0_trt)    <- contr.sum
  contrasts(df$S1_trt)    <- contr.sum
  contrasts(df$Treatment) <- contr.sum
  contrasts(df$Rep)       <- contr.sum
  assign(d, df)
}

## -----------------------------------------------------------
## 2. Helper: fit a main-effects-only model and extract partial
##    eta-squared for PrePeak, S0_trt, S1_trt, Treatment
## -----------------------------------------------------------

get_effect_sizes <- function(model, n_obs, trait_label, region_label) {
  
  focal_terms <- c("PrePeak", "S0_trt", "S1_trt", "Treatment")
  is_glmm <- inherits(model, "glmerMod")
  
  if (!is_glmm) {
    # Gaussian lmer: get proper F-test (Kenward-Roger df) via car::Anova,
    # then real partial omega-squared via effectsize::F_to_omega2()
    aov_tab <- as.data.frame(Anova(model, type = "III", test.statistic = "F"))
    aov_tab$Term <- rownames(aov_tab)
    aov_tab <- aov_tab[aov_tab$Term %in% focal_terms, ]
    
    omega2_vals <- mapply(function(f, df1, df2) {
      effectsize::F_to_omega2(f, df1, df2)$Omega2_partial
    }, aov_tab$F, aov_tab$Df, aov_tab$Df.res)
    
  } else {

    aov_tab <- as.data.frame(Anova(model, type = "III"))
    aov_tab$Term <- rownames(aov_tab)
    aov_tab <- aov_tab[aov_tab$Term %in% focal_terms, ]
    
    omega2_vals <- mapply(function(chisq, df) {
      max(0, (chisq - df) / (chisq + n_obs))   # clamp at 0, same convention as omega2
    }, aov_tab$Chisq, aov_tab$Df)
  }
  
  data.frame(
    Region = region_label,
    Trait = trait_label,
    Term = aov_tab$Term,
    Eta2 = omega2_vals   # column name kept as Eta2 for downstream code compatibility;
    # values are now partial omega-squared (see note above for pollen)
  )
}


## -----------------------------------------------------------
## 3. Fit main-effects-only models: one per trait x region
## -----------------------------------------------------------

f_main_gaussian <- function(resp) as.formula(paste(resp, "~ PrePeak + S0_trt + S1_trt + Treatment + Rep + (1 | S0_ID)"))
f_main_count    <- function(resp) as.formula(paste(resp, "~ PrePeak + S0_trt + S1_trt + Treatment + Rep + (1 | S0_ID)"))
f_main_binomial <- cbind(Viable, Total - Viable) ~ PrePeak + S0_trt + S1_trt + Treatment + (1 | S0_ID) + (1 | Obs)
f_main_seed     <- Seed_total ~ PrePeak + S0_trt + S1_trt + Treatment + Rep + (1 | S0_ID)

results_list <- list()
diagnostics_list <- list()   # residual-normality check per Gaussian model

## -----------------------------------------------------------
## Growth/phenology traits (Gaussian): 8 traits, both regions.

gaussian_traits <- c("days_to_germ", "FD", "branches", "Diameter",
                     "Height", "nodes", "biomass", "sla")
trait_labels <- c(days_to_germ = "Germ", FD = "FlowerDate", branches = "BranchNum",
                  Diameter = "Diameter", Height = "Height", nodes = "LeafNum",
                  biomass = "Biomass", sla = "SLA")

for (trait in gaussian_traits) {
  for (region in c("North", "South")) {
    d <- if (region == "North") dat_n else dat_s
    

    if (trait == "biomass") {
      d <- filter(d, biomass > 0)
    }

    if (trait == "sla") {
      d <- filter(d, sla < 700)
    }
    

    if (trait == "nodes") {
      d <- filter(d, nodes == round(nodes), nodes <= 30)
    }
    

    if (trait == "Diameter") {
      d <- filter(d, Diameter <= 10)
    }

    trait_formula_resp <- switch(trait,
                                 days_to_germ = "log(days_to_germ)",
                                 biomass = "sqrt(biomass)",
                                 trait   # default: untransformed (includes sla)
    )
    
    m <- lmer(f_main_gaussian(trait_formula_resp), data = d)
    
    resids <- resid(m)
    sw <- if (length(resids) <= 5000) shapiro.test(resids) else NULL
    
    diagnostics_list[[paste(trait, region)]] <- data.frame(
      Trait = trait_labels[trait], Region = region,
      shapiro_p = if (!is.null(sw)) sw$p.value else NA,
      resid_skew = (mean((resids - mean(resids))^3) / sd(resids)^3)  # rough skewness
    )
    
    results_list[[paste(trait, region)]] <- get_effect_sizes(m, nrow(d), trait_labels[trait], region)
  }
}

diagnostics_df <- bind_rows(diagnostics_list)
diagnostics_df


## -----------------------------------------------------------
## Flower number (count: Poisson -> negative binomial 
## -----------------------------------------------------------

for (region in c("North", "South")) {
  d <- if (region == "North") dat_n else dat_s
  
  m_pois <- glmer(f_main_count("flower_total"), family = poisson, data = d)
  sim_pois <- simulateResiduals(m_pois, n = 1000)
  disp_test <- testDispersion(sim_pois, plot = FALSE)
  
  cat("\n--- FlowerNum", region, ": Poisson dispersion test ---\n")
  print(disp_test)

  m_final <- if (disp_test$p.value < 0.05) {
    glmer.nb(f_main_count("flower_total"), data = d)
  } else {
    m_pois
  }
  
  results_list[[paste("flower_total", region)]] <- get_effect_sizes(m_final, nrow(d), "FlowerNum", region)
}

## --- Pollen viability (binomial) ---
for (region in c("North", "South")) {
  d <- if (region == "North") pollen_n else pollen_s
  m <- glmer(f_main_binomial, family = binomial, data = d)
  results_list[[paste("pollen", region)]] <- get_effect_sizes(m, nrow(d), "Pollen", region)
}

## --- Seed set (Gaussian -- confirmed appropriate via earlier Poisson/NB checks) ---
for (region in c("North", "South")) {
  d <- if (region == "North") seed_n else seed_s
  m <- lmer(f_main_seed, data = d)
  results_list[[paste("seed", region)]] <- get_effect_sizes(m, nrow(d), "Seed Set", region)
}

effect_size_df <- bind_rows(results_list)

effect_size_df$Term <- factor(effect_size_df$Term,
                              levels = c("PrePeak", "S0_trt", "S1_trt", "Treatment"),
                              labels = c("Year", "Grandparent Treatment", "Parent Treatment", "Current Treatment"))

effect_size_df$Trait <- factor(effect_size_df$Trait,
                               levels = c("Germ", "LeafNum", "BranchNum", "Diameter", "Height", "Biomass",
                                          "SLA", "FlowerNum", "FlowerDate", "Pollen", "Seed Set"))

effect_size_df


## -----------------------------------------------------------
## Effect-size table (all 11 traits, additive main-effects-only
## models) -- 
## -----------------------------------------------------------

effect_size_wide <- effect_size_df %>%
  mutate(Eta2_display = sprintf("%.3f", Eta2)) %>%
  dplyr::select(Trait, Region, Term, Eta2_display) %>%
  pivot_wider(names_from = Term, values_from = Eta2_display) %>%
  arrange(Trait, Region)

effect_size_wide

effect_size_gt <- effect_size_wide %>%
  gt() %>%
  tab_header(title = "Supplemental Table: Effect Sizes (Partial Omega-Squared) Across All 11 Traits") %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title()) %>%
  tab_footnote(
    footnote = "Partial omega-squared (see Part 4 methods note for justification over eta-squared) from main-effects-only models (Year/PrePeak + Grandparent Treatment/S0_trt + Parent Treatment/S1_trt + Current Treatment/Treatment; interactions excluded by design -- see Part 4). Gaussian traits: F_to_omega2() on Kenward-Roger F-tests. Pollen (binomial GLMM): chi-square-based bias-corrected analog, an approximation (see Part 4 note) -- not directly comparable to the Gaussian traits' values.",
    locations = cells_title()
  )

effect_size_gt
gtsave(effect_size_gt, "Supplemental_EffectSize_table.docx")


## -----------------------------------------------------------
## 4. Plot: grouped bar chart, North and South as separate panels
## -----------------------------------------------------------

term_colors <- c("Year" = "#F0A500", "Grandparent Treatment" = "#D6336C",
                 "Parent Treatment" = "#495867", "Current Treatment" = "#748CAB")

y_cap <- 0.1
y_axis_max <- 0.12   # headroom above the cap, for labels only -- bars
# themselves are capped in the data (below) and
# can never reach into this space
dodge_width <- 0.8
bar_width <- 0.75
slot_width <- dodge_width / 4

# cap bar HEIGHT in the data itself 
effect_size_df$Eta2_capped <- pmin(effect_size_df$Eta2, y_cap)

## labels for bars exceeding the cap
make_overflow_labels <- function(df) {
  df %>%
    filter(Eta2 > y_cap) %>%
    mutate(
      Trait_num = as.numeric(Trait),
      term_pos = as.numeric(Term),
      x_pos = Trait_num + (term_pos - 2.5) * slot_width,
      label = sprintf("%.2f", Eta2)
    )
}

overflow_N <- make_overflow_labels(filter(effect_size_df, Region == "North"))
overflow_S <- make_overflow_labels(filter(effect_size_df, Region == "South"))

p_effect_N <- ggplot(filter(effect_size_df, Region == "North"),
                     aes(x = Trait, y = Eta2_capped, fill = Term)) +
  geom_col(position = position_dodge(dodge_width), width = bar_width, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = y_cap, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0.06, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_text(data = overflow_N, aes(x = x_pos, y = y_cap + 0.008, label = label),
            inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = term_colors) +
  scale_y_continuous(limits = c(0, y_axis_max)) +
  labs(
    x = "Trait", 
    y = expression(Effect~size~(partial~italic(omega)^2)), 
    title = "North"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

p_effect_S <- ggplot(filter(effect_size_df, Region == "South"),
                     aes(x = Trait, y = Eta2_capped, fill = Term)) +
  geom_col(position = position_dodge(dodge_width), width = bar_width, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = y_cap, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0.06, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_text(data = overflow_S, aes(x = x_pos, y = y_cap + 0.008, label = label),
            inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = term_colors) +
  scale_y_continuous(limits = c(0, y_axis_max)) +
  labs(
    x = "Trait", 
    y = expression(Effect~size~(partial~italic(omega)^2)), 
    title = "South"
  )+
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

effect_size_fig <- (p_effect_N / p_effect_S) +
  plot_annotation(tag_levels = list(c("A", "B")))

effect_size_fig
ggsave("effect_size_by_trait_region.png", effect_size_fig, width = 8, height = 9, dpi = 300)



### =========================================================
### PART 6: Grandparent x Parent x Current treatment interaction
### =========================================================

### -----------------------------------------------------------
### 1. Build region-split pre-drought datasets
### -----------------------------------------------------------

dat5_n     <- filter(dat5, Region == "North")
dat5_s     <- filter(dat5, Region == "South")
dat5_sla_n <- filter(dat5_sla, Region == "North")
dat5_sla_s <- filter(dat5_sla, Region == "South")

pollen_pre_n <- filter(pollen_pre, Region == "North")
pollen_pre_s <- filter(pollen_pre, Region == "South")

seed_pre_n <- filter(seed_pre, Region == "North")
seed_pre_s <- filter(seed_pre, Region == "South")


### -----------------------------------------------------------
### 2. Fit pipeline: S0_trt * S1_trt * Treatment + Rep, per trait,
###    per region. 
### -----------------------------------------------------------

fit_3way_gaussian <- function(formula, data) {
  m <- lmer(formula, data = data)
  list(
    model = m,
    singular = isSingular(m),
    anova_satterthwaite = anova(m),
    anova_KR = Anova(m, type = "III", test.statistic = "F")
  )
}

## --- PC1 ---
res3_PC1_N <- fit_3way_gaussian(PC1 ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_n)
res3_PC1_S <- fit_3way_gaussian(PC1 ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_s)

## --- Germination ---
res3_germ_N <- fit_3way_gaussian(log(days_to_germ) ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_n)
res3_germ_S <- fit_3way_gaussian(log(days_to_germ) ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_s)

## --- Flowering date (untransformed, per earlier decision) ---
res3_FD_N <- fit_3way_gaussian(FD ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_n)
res3_FD_S <- fit_3way_gaussian(FD ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_s)

## --- SLA (sqrt-transformed, cleaned data) ---
res3_sla_N <- fit_3way_gaussian(sla ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_sla_n)
res3_sla_S <- fit_3way_gaussian(sla ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), dat5_sla_s)

## --- Pollen viability (binomial, OLRE for overdispersion) ---
fit_3way_pollen <- function(data) {
  m <- glmer(cbind(Viable, Total - Viable) ~ S0_trt * S1_trt * Treatment + (1 | S0_ID) + (1 | Obs),
             family = binomial, data = data)
  sim <- simulateResiduals(m, n = 1000)
  list(model = m, singular = isSingular(m), anova_wald = Anova(m, type = "III"),
       dispersion_test = testDispersion(sim, plot = FALSE))
}

res3_pollen_N <- fit_3way_pollen(pollen_pre_n)
res3_pollen_S <- fit_3way_pollen(pollen_pre_s)

## --- Seed set (Gaussian) ---
res3_seed_N <- fit_3way_gaussian(Seed_total ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), seed_pre_n)
res3_seed_S <- fit_3way_gaussian(Seed_total ~ S0_trt * S1_trt * Treatment + Rep + (1 | S0_ID), seed_pre_s)

## --- Inspect singularity + dispersion before trusting any result ---
list(PC1_N = res3_PC1_N$singular, PC1_S = res3_PC1_S$singular,
     germ_N = res3_germ_N$singular, germ_S = res3_germ_S$singular,
     FD_N = res3_FD_N$singular, FD_S = res3_FD_S$singular,
     sla_N = res3_sla_N$singular, sla_S = res3_sla_S$singular,
     pollen_N = res3_pollen_N$singular, pollen_S = res3_pollen_S$singular,
     seed_N = res3_seed_N$singular, seed_S = res3_seed_S$singular)

res3_pollen_N$dispersion_test
res3_pollen_S$dispersion_test

### -----------------------------------------------------------
### 4. Full results table: ALL terms (not just the 3-way), all 6
###    traits x 2 regions
### -----------------------------------------------------------

extract_all_terms <- function(anova_tab, trait_label, region_label) {
  df <- as.data.frame(anova_tab)
  df$Term <- rownames(df)
  df <- df[df$Term != "(Intercept)", ]
  
  # statistic column name differs by source: "Chisq" (Wald, glmer),
  # "F value" (lmerTest's native anova() -- true Satterthwaite), or
  # "F" (car::Anova with test.statistic="F" -- Kenward-Roger)
  stat_col <- intersect(c("Chisq", "F value", "F"), names(df))[1]
  p_col <- intersect(c("Pr(>Chisq)", "Pr(>F)"), names(df))[1]
  # Df column also differs: plain "Df" (car::Anova, both Wald and KR)
  # vs "NumDF"/"DenDF" (lmerTest's native anova()) -- report NumDF
  # where DenDF isn't available under a single "Df" name for display
  df_col <- if ("Df" %in% names(df)) "Df" else "NumDF"
  
  data.frame(
    Trait = trait_label,
    Region = region_label,
    Term = df$Term,
    Statistic = df[[stat_col]],
    Df = df[[df_col]],
    DenDF = if ("DenDF" %in% names(df)) round(df$DenDF, 1) else NA,
    p_value = df[[p_col]]
  )
}

full_results_table <- bind_rows(
  extract_all_terms(res3_PC1_N$anova_satterthwaite, "PC1", "North"),
  extract_all_terms(res3_PC1_S$anova_satterthwaite, "PC1", "South"),
  extract_all_terms(res3_germ_N$anova_satterthwaite, "Germination", "North"),
  extract_all_terms(res3_germ_S$anova_satterthwaite, "Germination", "South"),
  extract_all_terms(res3_FD_N$anova_satterthwaite, "FloweringDate", "North"),
  extract_all_terms(res3_FD_S$anova_satterthwaite, "FloweringDate", "South"),
  extract_all_terms(res3_sla_N$anova_satterthwaite, "SLA", "North"),
  extract_all_terms(res3_sla_S$anova_satterthwaite, "SLA", "South"),
  extract_all_terms(res3_pollen_N$anova_wald, "Pollen", "North"),
  extract_all_terms(res3_pollen_S$anova_wald, "Pollen", "South"),
  extract_all_terms(res3_seed_N$anova_satterthwaite, "SeedSet", "North"),
  extract_all_terms(res3_seed_S$anova_satterthwaite, "SeedSet", "South")
) %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ ".",
      TRUE ~ ""
    ),
    Statistic = round(Statistic, 2)
  ) %>%
  mutate(p_display = ifelse(p_value < 0.001, "<.001", as.character(round(p_value, 3))))

full_results_table

full_results_wide <- full_results_table %>%
  mutate(cell = paste0(Statistic, " (p=", p_display, sig, ")")) %>%
  dplyr::select(Trait, Region, Term, cell) %>%
  pivot_wider(names_from = Term, values_from = cell)

full_results_wide

# write_csv(full_results_table, "part5_full_results_long.csv")
# write_csv(full_results_wide, "part5_full_results_wide.csv")



### -----------visualization 

## --- Flowering date, South -- S0_trt x S1_trt interaction ---
emm_FD_S_plot <- as.data.frame(emmeans(res3_FD_S$model, ~ S0_trt | S1_trt))

p_FD_South_interaction <- ggplot(emm_FD_S_plot, aes(x = S0_trt, y = emmean, color = S1_trt, group = S1_trt)) +
  geom_point(position = position_dodge(0.2), size = 3) +
  geom_line(aes(linetype = S1_trt), position = position_dodge(0.2), linewidth = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(0.2)) +
  scale_color_manual(values = c("D" = "#8B4513", "W" = "#1E90FF")) +
  scale_linetype_manual(values = c("D" = "solid", "W" = "dotted")) +
  labs(x = "Grandparent treatment", y = "Flowering date",
       title = "South: Grandparent x Parent treatment interaction",
       color = "Parent\ntreatment", linetype = "Parent\ntreatment") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

p_FD_South_interaction
ggsave("FD_South_S0xS1_interaction.png", p_FD_South_interaction, width = 5.5, height = 4.5, dpi = 300)

### =========================================================
### PART 7: Supplemental table -- LRT justification for additive
### =========================================================

## -----------------------------------------------------------
## Helper: extract Chisq, Df, p-value from an anova() LRT object
## -----------------------------------------------------------

extract_lrt <- function(lrt_obj, trait_label, region_label) {
  df <- as.data.frame(lrt_obj)
  # the LRT row (second model) holds the Chisq/Df/p-value; the first
  # model's row is always NA on these columns
  row <- df[2, ]
  data.frame(
    Trait = trait_label,
    Region = region_label,
    Chisq = round(row$Chisq, 2),
    Df = row$Df,
    p_value = row[["Pr(>Chisq)"]]
  )
}

lrt_summary <- bind_rows(
  extract_lrt(res_PC1$lrt,  "PC1",           "Unified"),
  extract_lrt(res_germ$lrt, "Germination",   "Unified"),
  extract_lrt(res_FD$lrt,   "Flowering Date","Unified"),
  extract_lrt(res_sla$lrt,  "SLA",           "Unified"),
  extract_lrt(lrt_pre,         "Pollen",     "Baseline"),
  extract_lrt(res_pollen_N$lrt,"Pollen",     "North"),
  extract_lrt(res_pollen_S$lrt,"Pollen",     "South"),
  extract_lrt(res_seed$lrt,    "Seed Set",   "Baseline"),
  extract_lrt(res_seed_N$lrt,  "Seed Set",   "North"),
  extract_lrt(res_seed_S$lrt,  "Seed Set",   "South"),
  # Part 10: region-split, PrePeak-inclusive follow-up models --
  # completes LRT coverage for every model in the study, per methods.
  extract_lrt(res4_PC1_N$lrt,  "PC1",            "North (PrePeak)"),
  extract_lrt(res4_PC1_S$lrt,  "PC1",            "South (PrePeak)"),
  extract_lrt(res4_germ_N$lrt, "Germination",    "North (PrePeak)"),
  extract_lrt(res4_germ_S$lrt, "Germination",    "South (PrePeak)"),
  extract_lrt(res4_FD_N$lrt,   "Flowering Date", "North (PrePeak)"),
  extract_lrt(res4_FD_S$lrt,   "Flowering Date", "South (PrePeak)"),
  extract_lrt(res4_sla_N$lrt,  "SLA",            "North (PrePeak)"),
  extract_lrt(res4_sla_S$lrt,  "SLA",            "South (PrePeak)")
) %>%
  mutate(
    p_display = ifelse(p_value < 0.001, "<.001", sprintf("%.3f", p_value)),
    sig = case_when(p_value < .001 ~ "***", p_value < .01 ~ "**", p_value < .05 ~ "*", p_value < .10 ~ ".", TRUE ~ "n.s."),
    Conclusion = ifelse(p_value >= .05, "Additive S1_trt supported", "Interactive S1_trt NOT rejected"),
    # override for SLA North (PrePeak): p=.053 is conventionally >=
    # .05 (would default to "Additive supported"), but the analyst
    # judged this borderline result sufficient to prefer the
    # interactive structure for this specific cell -- overriding the
    # automatic rule with an explicit note rather than silently
    # relabeling. This is the ONLY model in the entire analysis using
    # the interactive S1_trt structure; every other cell (17 of 18
    # comparisons here) uses additive.
    Conclusion = ifelse(Trait == "SLA" & Region == "North (PrePeak)",
                        "Interactive S1_trt USED (analyst override, p=.053 judged sufficient -- see res4_sla_N_interactive)",
                        Conclusion)
  )

lrt_summary

## -----------------------------------------------------------
## Formatted gt table
## -----------------------------------------------------------

lrt_gt <- lrt_summary %>%
  dplyr::select(Trait, Region, Chisq, Df, p_display, sig, Conclusion) %>%
  gt() %>%
  tab_header(title = "Supplemental Table: Likelihood Ratio Tests for Additive vs. Fully Interactive S1_trt Structure") %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title()) %>%
  cols_label(
    p_display = "p-value",
    sig = "",
    Conclusion = "Interpretation"
  ) %>%
  tab_footnote(
    footnote = "Each row compares a model with S1_trt (parent treatment) entered additively against the same model with S1_trt fully crossed with all other fixed effects (chi-square test on twice the log-likelihood difference). Non-significant results (all rows here) support treating S1_trt as an additive covariate rather than allowing it to interact -- the structure used throughout this analysis for every applicable trait/model. *** p<.001, ** p<.01, * p<.05, . p<.10.",
    locations = cells_title()
  )

lrt_gt
#gtsave(lrt_gt, "Supplemental_LRT_table.docx")
