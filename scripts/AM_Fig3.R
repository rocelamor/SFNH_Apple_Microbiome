################################################################################
################################################################################
##
## AM_Fig3
## Bacterial abundance classes, diversity/richness responses, and matched-control richness deficit
##
## A  ASV abundance classification
## B  Partial-Spearman diversity/richness associations
## C  Negative-binomial GLMM for rare and total ASV richness
## D  Richness deficit relative to matched uninoculated controls
##
################################################################################
################################################################################


## =============================================================================
## 0. PACKAGES
## =============================================================================

library(tidyverse)
library(vegan)
library(patchwork)
library(scales)
library(grid)
library(lme4)      # glmer.nb() — negative-binomial GLMM, ships with lme4 (no glmmTMB/TMB compile needed)
library(DHARMa)


## =============================================================================
## PATHS
## =============================================================================

## Run this script with the working directory set to the repository root
## (open the repo as an RStudio project, or run `Rscript scripts/<name>.R`
## from the repo root). Place the data files in ./data ; outputs are
## written to ./output/Figures and ./output/Tables.
input_dir <- "data"

output_dir <- "output"

figures_dir <- file.path(output_dir, "Figures")
tables_dir  <- file.path(output_dir, "Tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


## =============================================================================
## 2. INPUT FILES
## =============================================================================

metadata_file <- file.path(
  input_dir,
  "metadata.csv"
)

counts_file <- file.path(
  input_dir,
  "bac_counts_rarefied14000.csv"
)

class_file <- file.path(
  input_dir,
  "bac_abundance_class.csv"
)


## =============================================================================
## 3. READ DATA
## =============================================================================

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cnt_ref <- read.csv(
  counts_file,
  row.names = 1,
  check.names = FALSE
) |>
  as.matrix()

storage.mode(
  cnt_ref
) <- "numeric"

abcl <- read.csv(
  class_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


## =============================================================================
## 4. FORMAT METADATA
## =============================================================================

meta <- meta |>
  mutate(
    ES_Number = factor(
      ES_Number,
      levels = c("ES1", "ES2", "ES3")
    ),
    
    Collection_day = factor(
      Collection_day,
      levels = c("12D", "15D", "18D", "21D")
    )
  )


## =============================================================================
## 5. APPLE / WORM SAMPLES (n = 36) — used for Panels A-D and the GLMM
## =============================================================================

W <- meta |>
  filter(
    Compartment == "As",
    Treatment == "Worm",
    SampleID %in% colnames(cnt_ref)
  ) |>
  arrange(
    ES_Number,
    Collection_day,
    Pot_Replicate
  ) |>
  mutate(
    apple_worms = worm,
    cell_id = interaction(
      ES_Number,
      Collection_day,
      drop = TRUE
    )
  )

cnt_apple <- cnt_ref[, W$SampleID, drop = FALSE]

stopifnot(
  identical(W$SampleID, colnames(cnt_apple))
)

cat("\nApple/nematode-inoculated samples:", nrow(W), "\n")
cat("Cells (should be 12):", n_distinct(W$cell_id), "\n")


## =============================================================================
## 6. ABUNDANCE CLASSES
## =============================================================================

rare_asv     <- abcl$ASV_ID[abcl$class == "Rare"]
moderate_asv <- abcl$ASV_ID[abcl$class == "Moderate"]
crt_asv      <- abcl$ASV_ID[abcl$class == "Conditionally rare/abundant"]

rare_asv_present     <- intersect(rare_asv, rownames(cnt_apple))
moderate_asv_present <- intersect(moderate_asv, rownames(cnt_apple))
crt_asv_present      <- intersect(crt_asv, rownames(cnt_apple))


## =============================================================================
## 7. DIVERSITY / RICHNESS METRICS (Worm samples only, n = 36)
## =============================================================================

W$Shannon <- vegan::diversity(t(cnt_apple), index = "shannon")

W$Total_richness <- colSums(cnt_apple > 0)

W$Rare_richness <- colSums(
  cnt_apple[rare_asv_present, , drop = FALSE] > 0
)

W$Moderate_richness <- colSums(
  cnt_apple[moderate_asv_present, , drop = FALSE] > 0
)

W$CRT_richness <- colSums(
  cnt_apple[crt_asv_present, , drop = FALSE] > 0
)


## =============================================================================
## 8. PARTIAL SPEARMAN FUNCTION (unchanged from original script)
## =============================================================================

partial_spearman <- function(y, x, covars) {
  
  keep <- complete.cases(y, x, covars)
  
  y <- y[keep]
  x <- x[keep]
  covars <- covars[keep, , drop = FALSE]
  
  mm <- model.matrix(~ ., data = as.data.frame(covars))
  
  ry <- residuals(lm(rank(y) ~ mm - 1))
  rx <- residuals(lm(rank(x) ~ mm - 1))
  
  ct <- cor.test(ry, rx, method = "pearson")
  
  c(
    rho = unname(ct$estimate),
    p   = ct$p.value,
    n   = length(y)
  )
}

cov_ed <- W[, c("ES_Number", "Collection_day")]

worm_predictor <- log1p(W$apple_worms)


## =============================================================================
## 9. COLORS + THEME — matched to current Figure 4
## =============================================================================

ES_COL <- c(
  "ES1" = "#264653",
  "ES2" = "#2a9d8f",
  "ES3" = "#e76f51"
)

theme_sfnh <- function(base_size = 20) {
  
  theme_classic(base_size = base_size, base_family = "sans") +
    
    theme(
      text = element_text(family = "sans", colour = "black"),
      axis.title = element_text(family = "sans", size = base_size + 1, colour = "black"),
      axis.text  = element_text(family = "sans", size = base_size, colour = "black"),
      axis.line  = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.7, colour = "black"),
      legend.title = element_blank(),
      legend.text  = element_text(family = "sans", size = base_size - 1),
      plot.margin  = margin(10, 14, 10, 14)
    )
}


################################################################################
################################################################################
##
## PANEL A — BACTERIAL ASV ABUNDANCE CLASSIFICATION (unchanged from original)
##
################################################################################
################################################################################

class_summary <- abcl |>
  group_by(class) |>
  summarise(n_asv = n(), .groups = "drop") |>
  mutate(pct_asv = 100 * n_asv / nrow(cnt_ref))

reads_per_asv <- tibble(
  ASV_ID = rownames(cnt_ref),
  reads  = rowSums(cnt_ref)
)

reads_by_class <- abcl |>
  select(ASV_ID, class) |>
  left_join(reads_per_asv, by = "ASV_ID") |>
  group_by(class) |>
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") |>
  mutate(pct_reads = 100 * reads / sum(cnt_ref))

fig3a_table <- class_summary |>
  left_join(
    reads_by_class |> select(class, pct_reads),
    by = "class"
  ) |>
  mutate(
    class = factor(
      class,
      levels = c("Rare", "Moderate", "Conditionally rare/abundant", "Abundant")
    )
  ) |>
  arrange(class)

write.csv(
  fig3a_table,
  file.path(tables_dir, "AM_Fig3A_Abundance_Classification.csv"),
  row.names = FALSE
)

p3a_data <- fig3a_table |>
  filter(class != "Abundant") |>
  mutate(class = droplevels(class)) |>
  pivot_longer(
    cols = c(pct_asv, pct_reads),
    names_to = "metric",
    values_to = "pct"
  ) |>
  mutate(
    metric = recode(metric, pct_asv = "ASVs", pct_reads = "Sequence reads"),
    metric = factor(metric, levels = c("ASVs", "Sequence reads")),
    label = if_else(
      metric == "ASVs",
      scales::comma(n_asv),
      paste0(sprintf("%.1f", pct), "%")
    )
  )

max_pct <- max(p3a_data$pct, na.rm = TRUE)
label_offset <- max_pct * 0.06

p3a <- ggplot(p3a_data, aes(x = class, y = pct, fill = metric)) +
  
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.68, colour = "black", linewidth = 0.45
  ) +
  
  geom_text(
    aes(y = pct + label_offset, label = label),
    position = position_dodge(width = 0.72),
    family = "sans", size = 5.0, hjust = 0
  ) +
  
  scale_fill_manual(
    values = c("ASVs" = "grey20", "Sequence reads" = "grey70"),
    name = NULL
  ) +
  
  scale_x_discrete(
    labels = c(
      "Rare" = "Rare",
      "Moderate" = "Moderate",
      "Conditionally rare/abundant" = "Conditionally\nrare/abundant"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, max_pct * 1.55),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_flip(clip = "off") +
  
  labs(x = NULL, y = "Proportion (%)") +
  
  theme_sfnh(base_size = 20) +
  
  theme(
    legend.position = c(0.82, 0.60),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 15),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(margin = margin(t = 4)),
    plot.margin = margin(8, 10, 4, 10)
  )


################################################################################
################################################################################
##
## PANEL B — PARTIAL SPEARMAN EFFECT SIZES, NOW WITH BH CORRECTION
##
################################################################################
################################################################################

run_metric <- function(metric, y) {
  
  r <- partial_spearman(y, worm_predictor, cov_ed)
  
  tibble(
    metric = metric,
    rho = as.numeric(r["rho"]),
    p   = as.numeric(r["p"]),
    n   = as.numeric(r["n"])
  )
}

effect_results <- bind_rows(
  run_metric("Shannon", W$Shannon),
  run_metric("Total richness", W$Total_richness),
  run_metric("Cond. rare/abundant", W$CRT_richness),
  run_metric("Moderate", W$Moderate_richness),
  run_metric("Rare", W$Rare_richness)
) |>
  mutate(
    ## Benjamini-Hochberg correction across the five metrics.
    q = p.adjust(p, method = "BH"),
    
    metric = factor(
      metric,
      levels = c("Shannon", "Total richness", "Cond. rare/abundant", "Moderate", "Rare")
    ),
    
    ## Stars now reflect the BH-adjusted q-value, not the raw P.
    stars = case_when(
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      TRUE      ~ ""
    )
  )

cat(
  "\n============================================================\n",
  "AM_Fig3B — EFFECT SIZES (BH-corrected)\n",
  "============================================================\n"
)

print(effect_results)

write.csv(
  effect_results,
  file.path(tables_dir, "AM_Fig3B_Effect_Sizes_BH.csv"),
  row.names = FALSE
)

p3b <- ggplot(effect_results, aes(x = metric, y = rho)) +
  
  geom_hline(yintercept = 0, linewidth = 0.6, colour = "black") +
  
  geom_col(width = 0.65, fill = "grey25", colour = "black", linewidth = 0.45) +
  
  geom_text(
    aes(label = stars, y = rho - 0.035),
    family = "sans", size = 7
  ) +
  
  annotate(
    "text", x = 0.65, y = 0.10, hjust = 0, size = 4.6, family = "sans",
    label = "*  q<0.05    **  q<0.01    ***  q<0.001",
    colour = "grey20"
  ) +
  
  scale_y_continuous(
    limits = c(-0.75, 0.15),
    breaks = seq(-0.6, 0, by = 0.2)
  ) +
  
  labs(x = NULL, y = expression(Partial~Spearman~rho)) +
  
  theme_sfnh(base_size = 20) +
  
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 17),
    legend.position = "none"
  )


################################################################################
################################################################################
##
## PANEL C — RARE ASV RICHNESS vs ACTUAL APPLE WORM ABUNDANCE
##
## Annotation statistic is unchanged (rank-based partial Spearman).
## The dashed line is now an ES + day - adjusted (ANCOVA) fitted curve,
## averaged across the 12 observed ES x day combinations, instead of an
## unadjusted OLS line that ignored ES and day entirely.
##
################################################################################
################################################################################

rare_res <- partial_spearman(W$Rare_richness, worm_predictor, cov_ed)

rare_label <- paste0(
  "rho = ", sprintf("%.3f", rare_res["rho"]),
  "\nP ", ifelse(
    rare_res["p"] < 0.001,
    "< 0.001",
    paste0("= ", sprintf("%.3f", rare_res["p"]))
  ),
  "\nn = ", as.integer(rare_res["n"])
)

## ---- ANCOVA-adjusted fitted curve for the dashed line ----------------------

ancova_fit <- lm(
  Rare_richness ~ worm_predictor + ES_Number + Collection_day,
  data = W
)

## Prediction grid across the observed abundance range, averaged over the
## 12 ES x day combinations (balanced design, so an unweighted average of
## predictions across all combinations is the population-averaged curve).

worm_grid <- seq(
  min(worm_predictor),
  max(worm_predictor),
  length.out = 200
)

es_day_grid <- expand.grid(
  ES_Number = levels(W$ES_Number),
  Collection_day = levels(W$Collection_day)
)

pred_grid_c <- expand.grid(
  worm_predictor = worm_grid,
  ES_Number = levels(W$ES_Number),
  Collection_day = levels(W$Collection_day)
) |>
  mutate(
    ES_Number = factor(ES_Number, levels = levels(W$ES_Number)),
    Collection_day = factor(Collection_day, levels = levels(W$Collection_day))
  )

pred_grid_c$fit <- predict(ancova_fit, newdata = pred_grid_c)

adjusted_curve <- pred_grid_c |>
  group_by(worm_predictor) |>
  summarise(fit = mean(fit), .groups = "drop") |>
  mutate(apple_worms = expm1(worm_predictor))

## NOTE: this raw-scatter + ANCOVA-adjusted-line panel is now a
## SUPPLEMENTARY exhibit, not the main-figure Panel C — see the "PANEL C
## (MAIN FIGURE) — GLMM DOSE-RESPONSE" block further down, which replaces
## it as Panel C. Kept here (renamed) as p3c_supp_rare_scatter because it
## is still a legitimate, simpler univariate view of the Rare-richness
## trend and may be useful for reviewers or a supplementary figure.
p3c_supp_rare_scatter <- ggplot(W, aes(x = apple_worms + 1, y = Rare_richness, colour = ES_Number)) +
  
  geom_point(size = 4.6, alpha = 0.82) +
  
  geom_line(
    data = adjusted_curve,
    aes(x = apple_worms + 1, y = fit),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.1,
    linetype = 2
  ) +
  
  annotate(
    "label", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.15,
    label = rare_label, family = "sans", size = 5.3, linewidth = 0.4, fill = "white"
  ) +
  
  scale_x_log10(
    labels = scales::label_number(),
    expand = expansion(mult = c(0.03, 0.35))
  ) +
  
  scale_colour_manual(values = ES_COL, name = NULL) +
  
  labs(
    x = "Nematodes recovered from apple\n(+1, log scale)",
    y = "Rare ASV richness"
  ) +
  
  theme_sfnh(base_size = 20) +
  
  theme(
    legend.position = c(0.97, 0.78),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_blank(),
    legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
    axis.title.x = element_text(margin = margin(t = 5)),
    plot.margin = margin(8, 10, 6, 10)
  )


################################################################################
################################################################################
##
## ADDITIONAL ANALYSIS 2
##
## MIXED-EFFECTS (GLMM) DOSE-RESPONSE MODEL
##
## Richness ~ log(nematode abundance) + (1 | ES x day cell)
##
## Negative-binomial family, matching the framework already used elsewhere
## in the manuscript for nematode-abundance counts. Random intercept per
## cell accounts for the nesting (3 pots per cell) that the partial-
## Spearman approach does not model.
##
## NOTE: fit with lme4::glmer.nb() rather than glmmTMB(), since glmmTMB
## was not available in the environment this was run in. glmer.nb() gives
## the same negative-binomial GLMM with a random intercept; the only
## practical differences are (1) it profiles theta (the NB dispersion
## parameter) via repeated calls to glmer() rather than estimating it
## jointly, which is slightly slower but numerically fine for a dataset
## this size, and (2) predict() on a merMod object doesn't return
## se.fit directly, so the CI band below is built manually from the
## fixed-effects variance-covariance matrix (population-level curve,
## random effects averaged out — same quantity glmmTMB's re.form=NA
## would have given).
##
################################################################################
################################################################################

fit_total_glmm <- glmer.nb(
  Total_richness ~ worm_predictor + (1 | cell_id),
  data = W
)

fit_rare_glmm <- glmer.nb(
  Rare_richness ~ worm_predictor + (1 | cell_id),
  data = W
)

cat("\n============================================================\n")
cat("ADDITIONAL ANALYSIS 2 — GLMM SUMMARY (Total richness)\n")
cat("============================================================\n")
print(summary(fit_total_glmm))
cat("\nSingular fit?", isSingular(fit_total_glmm), "\n")
cat("Estimated theta (NB dispersion):", getME(fit_total_glmm, "glmer.nb.theta"), "\n")

cat("\n============================================================\n")
cat("ADDITIONAL ANALYSIS 2 — GLMM SUMMARY (Rare richness)\n")
cat("============================================================\n")
print(summary(fit_rare_glmm))
cat("\nSingular fit?", isSingular(fit_rare_glmm), "\n")
cat("Estimated theta (NB dispersion):", getME(fit_rare_glmm, "glmer.nb.theta"), "\n")

## ---- Residual diagnostics — check before trusting the P-values ------------

sim_total <- simulateResiduals(fit_total_glmm)
sim_rare  <- simulateResiduals(fit_rare_glmm)

pdf(file.path(figures_dir, "AM_Fig3C_GLMM_Residual_Diagnostics.pdf"), width = 10, height = 6)
plot(sim_total, title = "Total richness GLMM")
plot(sim_rare, title = "Rare richness GLMM")
dev.off()

## ---- Effect sizes in interpretable units -----------------------------------

extract_effect <- function(model, label) {
  
  ## lme4's summary(model)$coefficients is a plain matrix (no $cond list,
  ## that's a glmmTMB-specific structure)
  est <- summary(model)$coefficients["worm_predictor", ]
  
  beta <- est["Estimate"]
  se   <- est["Std. Error"]
  pval <- est["Pr(>|z|)"]
  
  ci_low  <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  
  fold_per_10x      <- exp(beta * log(10))
  fold_per_10x_low  <- exp(ci_low  * log(10))
  fold_per_10x_high <- exp(ci_high * log(10))
  
  tibble(
    metric = label,
    beta = beta,
    se = se,
    p = pval,
    pct_change_per_10x    = (fold_per_10x - 1) * 100,
    pct_change_ci_low     = (fold_per_10x_low  - 1) * 100,
    pct_change_ci_high    = (fold_per_10x_high - 1) * 100
  )
}

glmm_effects <- bind_rows(
  extract_effect(fit_total_glmm, "Total ASV richness"),
  extract_effect(fit_rare_glmm,  "Rare ASV richness")
)

cat("\nGLMM effect sizes (% change in richness per 10-fold abundance increase):\n")
print(glmm_effects)

write.csv(
  glmm_effects,
  file.path(tables_dir, "AM_Fig3C_GLMM_Effect_Sizes.csv"),
  row.names = FALSE
)

## ---- Fitted dose-response curve over raw data ------------------------------

pred_grid2 <- tibble(
  apple_worms = exp(
    seq(log(min(W$apple_worms) + 1), log(max(W$apple_worms) + 1), length.out = 200)
  ) - 1
) |>
  mutate(worm_predictor = log1p(apple_worms))

predict_with_ci <- function(model, newdata) {
  
  ## predict.merMod (lme4) has no se.fit argument, unlike glmmTMB's
  ## predict(). Build the population-level (re.form = NA) curve and its
  ## CI manually from the fixed-effects design matrix and vcov, which is
  ## the standard delta-method approach for merMod objects and gives the
  ## same quantity glmmTMB's se.fit would have: uncertainty in the fixed
  ## effects only, random effects averaged out.
  
  X <- model.matrix(~ worm_predictor, data = newdata)
  beta_hat <- fixef(model)
  V        <- vcov(model)
  
  eta     <- as.numeric(X %*% beta_hat)
  se_eta  <- sqrt(rowSums((X %*% V) * X))
  
  newdata |>
    mutate(
      fit      = exp(eta),
      fit_low  = exp(eta - 1.96 * se_eta),
      fit_high = exp(eta + 1.96 * se_eta)
    )
}

pred_total2 <- predict_with_ci(fit_total_glmm, pred_grid2) |>
  mutate(metric = "Total ASV richness")

pred_rare2 <- predict_with_ci(fit_rare_glmm, pred_grid2) |>
  mutate(metric = "Rare ASV richness")

pred_all2 <- bind_rows(pred_total2, pred_rare2)

raw_long2 <- W |>
  select(ES_Number, apple_worms, Total_richness, Rare_richness) |>
  pivot_longer(
    cols = c(Total_richness, Rare_richness),
    names_to = "metric",
    values_to = "richness"
  ) |>
  mutate(
    metric = recode(
      metric,
      Total_richness = "Total ASV richness",
      Rare_richness   = "Rare ASV richness"
    )
  )

p_add2 <- ggplot() +
  
  geom_ribbon(
    data = pred_all2,
    aes(x = apple_worms + 1, ymin = fit_low, ymax = fit_high),
    alpha = 0.15, fill = "grey30"
  ) +
  
  geom_line(
    data = pred_all2,
    aes(x = apple_worms + 1, y = fit),
    linewidth = 0.9, colour = "black"
  ) +
  
  geom_point(
    data = raw_long2,
    aes(x = apple_worms + 1, y = richness, colour = ES_Number),
    size = 2.8, alpha = 0.8
  ) +
  
  scale_colour_manual(values = ES_COL, name = NULL) +
  
  scale_x_log10(labels = scales::label_number()) +
  
  facet_wrap(~ metric, scales = "free_y") +
  
  labs(
    x = "Nematodes recovered from apple (+1, log scale)",
    y = "ASV richness",
    title = "Additional analysis 2: GLMM dose-response (random intercept per ES x day cell)"
  ) +
  
  theme_sfnh(base_size = 16) +
  
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    plot.title = element_text(size = 14, face = "plain")
  )

ggsave(
  file.path(figures_dir, "AM_Fig3C_GLMM_Dose_Response.pdf"),
  p_add2, width = 200, height = 110, units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_Fig3C_GLMM_Dose_Response.png"),
  p_add2, width = 200, height = 110, units = "mm", dpi = 600, bg = "white"
)


################################################################################
################################################################################
##
## PANEL C (MAIN FIGURE) — GLMM DOSE-RESPONSE, RESTYLED TO REPLACE THE OLD
## RARE-RICHNESS SCATTER
##
################################################################################
################################################################################

## Figure layout: legend and per-panel labels are placed inside each panel.
## richness" title must sit INSIDE the data frame, not float outside it
## the way an external facet strip or an outside legend would. So: (1)
## facet strip text is switched off entirely (strip.text/background =
## element_blank()) and replaced with an in-panel bold title drawn as a
## geom_text layer in the top-left corner of each facet; (2) the ES1/2/3
## colour legend is set to an inset position (legend.position = c(x, y))
## that lands inside the Total-richness facet (the right-hand panel,
## which has more empty space in its upper-right corner than the Rare
## panel does) rather than "top"/"right", which would sit outside the
## panels.
##
## BUGFIX: an earlier version of this block produced a spurious THIRD,
## empty facet (with the legend stranded inside it) because "metric" was
## a plain character column in some data frames and an ad-hoc factor in
## others, with the factor's levels pulled from `levels(pred_all2$metric)`
## — which is NULL for a character column, so it silently fell back to
## alphabetical ordering instead of actually matching. facet_wrap()
## unions the "metric" values across every layer's data, so any stray
## level anywhere (even one row) adds a whole extra panel. Fixed here by
## forcing every data frame that feeds this plot onto the exact same
## explicit 2-level factor, and by pinning `nrow = 1, drop = TRUE` on
## facet_wrap() so no unused level can ever add a panel again.
##
## Also switched the in-panel text from x = -Inf / y = Inf with vjust > 1
## (fragile inside facets — depends on clipping behaviour and is easy to
## push off-panel invisibly) to real, per-facet computed coordinates, so
## the title and stats block are guaranteed to land inside the visible
## data area regardless of each facet's free y-scale range.

metric_levels <- c("Rare ASV richness", "Total ASV richness")

pred_all2 <- pred_all2 |>
  mutate(metric = factor(as.character(metric), levels = metric_levels))

raw_long2 <- raw_long2 |>
  mutate(metric = factor(as.character(metric), levels = metric_levels))

## Per-facet placement — using REAL empty space, not manufactured
## headroom. An earlier version reserved space by expanding the y-axis
## 55% above the tallest point, which fixed the overlap but distorted
## the axis (data got compressed into the bottom half of the panel,
## e.g. Total richness axis stretched to 800 when the real max is ~420).
## Reverted that: axis expansion is back to a modest default. Instead,
## the label sits in the corner that is ALREADY empty in this
## particular dataset: in both panels, the few points that reach high
## richness are all at LOW abundance (the worm+1 = 1 column and the
## single low-abundance outlier), while HIGH abundance always comes
## with LOWER richness (that's the declining trend itself). So the
## top-right corner — high abundance, high richness — has essentially
## no points in either panel. Label block is right-aligned there
## (hjust = 1, anchored near x_max), just under the facet's own y_max
## rather than above it.
##
## Also switched geom_text -> geom_label for the title (solid white
## background, thin border) as a safety margin: if a stray point still
## ends up nearby, the label's opaque background sits on top of it
## rather than text and point tangling together illegibly. This is the
## same device the very first version of this panel (the ANCOVA-line
## rare-richness scatter, still kept as a supplementary figure) already
## used successfully.
facet_ranges <- raw_long2 |>
  group_by(metric) |>
  summarise(
    y_max  = max(richness, na.rm = TRUE),
    x_min  = min(apple_worms + 1, na.rm = TRUE),
    x_max  = max(apple_worms + 1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    x_pos   = x_max,
    y_title = y_max * 0.99,
    y_stats = y_max * 0.90
  )

glmm_facet_titles <- facet_ranges |>
  mutate(title = as.character(metric))

glmm_facet_stats <- glmm_effects |>
  mutate(
    metric = factor(as.character(metric), levels = metric_levels),
    label = paste0(
      sprintf("%+.1f", pct_change_per_10x), "% per 10× abundance",
      "\n95% CI [", sprintf("%+.1f", pct_change_ci_low), ", ",
      sprintf("%+.1f", pct_change_ci_high), "]",
      "\nP ", ifelse(p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", p)))
    )
  ) |>
  left_join(facet_ranges |> select(metric, x_pos, y_stats), by = "metric")

p3c <- ggplot() +
  
  geom_ribbon(
    data = pred_all2,
    aes(x = apple_worms + 1, ymin = fit_low, ymax = fit_high),
    alpha = 0.15, fill = "grey30"
  ) +
  
  geom_line(
    data = pred_all2,
    aes(x = apple_worms + 1, y = fit),
    linewidth = 1.1, colour = "black"
  ) +
  
  geom_point(
    data = raw_long2,
    aes(x = apple_worms + 1, y = richness, colour = ES_Number),
    size = 4.0, alpha = 0.82
  ) +
  
  ## In-panel facet title (replaces the external strip) — top-right,
  ## right-aligned, solid background for safety margin
  geom_label(
    data = glmm_facet_titles,
    aes(x = x_pos, y = y_title, label = title),
    hjust = 1, vjust = 1, family = "sans", fontface = "bold", size = 5.0,
    fill = "white", label.size = 0, label.padding = unit(0.15, "lines")
  ) +
  
  ## In-panel effect-size annotation, just below the title, same corner
  geom_label(
    data = glmm_facet_stats,
    aes(x = x_pos, y = y_stats, label = label),
    hjust = 1, vjust = 1, family = "sans", size = 5, lineheight = 0.95,
    fill = "white", label.size = 0, label.padding = unit(0.15, "lines")
  ) +
  
  scale_colour_manual(values = ES_COL, name = NULL) +
  
  ## Add extra right-side space inside each facet so the ES legend can
  ## sit within the panel without overlapping the observed data.
  scale_x_log10(
    labels = scales::label_number(),
    expand = expansion(
      mult = c(0.03, 0.26)
    )
  ) +
  
  facet_wrap(~ metric, scales = "free_y", nrow = 1, drop = TRUE) +
  
  labs(
    x = "Nematodes recovered from apple\n(+1, log scale)",
    y = "ASV richness"
  ) +
  
  theme_sfnh(base_size = 20) +
  
  theme(
    ## Strip switched off — facet titles are now drawn inside the panel
    strip.background = element_blank(),
    strip.text = element_blank(),
    
    ## Inset legend, inside the Total-richness (right) facet rather than
    ## outside the plot area.
    ##
    ## Placement history on this one: top-right at a fixed high fraction
    ## put it in a big manufactured empty margin (fixed now — axis no
    ## longer artificially expanded, see facet_ranges above). Bottom-
    ## right at a fixed low fraction landed on a real point (a high-
    ## abundance ES2 sample sits at roughly y = 200 out of a ~420 max in
    ## the Total panel, i.e. ~48% up the panel). Top-right is now also
    ## where the title/stats label block lives (see geom_label calls
    ## above), so the legend can't go there either.
    ##
    ## ES legend placed in the blank right-side space inside the
    ## Total ASV facet so it does not overlap plotted data.
    legend.position = c(0.985, 0.55),
    legend.justification = c(1, 0.5),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 15),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    
    axis.title.x = element_text(margin = margin(t = 5)),
    
    ## Panel C y-axis title spacing.
    ## The large gap is handled at the patchwork-layout level below,
    ## so an extreme negative margin is no longer needed here.
    axis.title.y = element_text(
      margin = margin(r = 2)
    ),
    
    plot.margin = margin(8, 10, 6, 4)
  )


################################################################################
################################################################################
##
################################################################################
################################################################################
##
## PANEL D — MATCHED-CONTROL RICHNESS DEFICIT
##
## Richness deficit = matched control richness - nematode-pot richness
## Positive values therefore indicate lower richness than the matched control.
##
################################################################################
################################################################################

## 2. GET ALL APPLE SAMPLES
##
## Unlike the original Fig. 4 analysis, this section includes both:
##
##   36 Worm apple samples
##   12 Control apple samples
##
## =============================================================================

apple_all <- meta |>
  filter(
    Compartment == "As",
    Treatment %in% c(
      "Worm",
      "Control"
    ),
    SampleID %in% colnames(
      cnt_ref
    )
  ) |>
  arrange(
    ES_Number,
    Collection_day,
    Treatment,
    Pot_Replicate
  ) |>
  mutate(
    cell_id = interaction(
      ES_Number,
      Collection_day,
      drop = TRUE,
      sep = "_"
    )
  )


cat(
  "\n============================================================\n",
  "MATCHED-CONTROL RICHNESS DEFICIT\n",
  "============================================================\n\n"
)

cat(
  "Apple samples found:",
  nrow(
    apple_all
  ),
  "\n"
)

print(
  table(
    apple_all$Treatment
  )
)


## Expected:
## Control = 12
## Worm    = 36

stopifnot(
  sum(
    apple_all$Treatment == "Control"
  ) == 12
)

stopifnot(
  sum(
    apple_all$Treatment == "Worm"
  ) == 36
)


## =============================================================================
## 3. APPLE COUNT MATRIX
## =============================================================================

cnt_apple_all <- cnt_ref[
  ,
  apple_all$SampleID,
  drop = FALSE
]


stopifnot(
  identical(
    apple_all$SampleID,
    colnames(
      cnt_apple_all
    )
  )
)


## =============================================================================
## 4. CALCULATE TOTAL AND RARE ASV RICHNESS
##
## Uses exactly the same abundance-class definition already used in Fig. 4.
##
## =============================================================================

rare_asv_all <- intersect(
  rare_asv,
  rownames(
    cnt_apple_all
  )
)


apple_all$Total_richness <- colSums(
  cnt_apple_all > 0
)


apple_all$Rare_richness <- colSums(
  cnt_apple_all[
    rare_asv_all,
    ,
    drop = FALSE
  ] > 0
)


## =============================================================================
## 5. EXTRACT THE MATCHED CONTROL FOR EACH ES x DAY CELL
##
## There should be exactly one control apple for each of the 12 cells.
##
## =============================================================================

control_richness <- apple_all |>
  filter(
    Treatment == "Control"
  ) |>
  select(
    cell_id,
    ES_Number,
    Collection_day,
    control_SampleID = SampleID,
    control_Total_richness = Total_richness,
    control_Rare_richness = Rare_richness
  )


control_check <- control_richness |>
  count(
    cell_id
  )


print(
  control_check
)


stopifnot(
  nrow(
    control_richness
  ) == 12
)

stopifnot(
  all(
    control_check$n == 1
  )
)


## =============================================================================
## 6. MATCH EACH WORM POT TO THE CONTROL FROM THE SAME ES x DAY CELL
## =============================================================================

deficit_data <- apple_all |>
  filter(
    Treatment == "Worm"
  ) |>
  mutate(
    apple_worms = worm
  ) |>
  left_join(
    control_richness |>
      select(
        cell_id,
        control_SampleID,
        control_Total_richness,
        control_Rare_richness
      ),
    by = "cell_id"
  ) |>
  mutate(
    
    Total_richness_deficit =
      control_Total_richness -
      Total_richness,
    
    Rare_richness_deficit =
      control_Rare_richness -
      Rare_richness,
    
    log_apple_worms =
      log1p(
        apple_worms
      )
  )


stopifnot(
  nrow(
    deficit_data
  ) == 36
)

stopifnot(
  !any(
    is.na(
      deficit_data$control_Total_richness
    )
  )
)

stopifnot(
  !any(
    is.na(
      deficit_data$control_Rare_richness
    )
  )
)


## =============================================================================
## 7. INSPECT MATCHING
## =============================================================================

deficit_check <- deficit_data |>
  select(
    SampleID,
    control_SampleID,
    ES_Number,
    Collection_day,
    apple_worms,
    Total_richness,
    control_Total_richness,
    Total_richness_deficit,
    Rare_richness,
    control_Rare_richness,
    Rare_richness_deficit
  )


print(
  as_tibble(deficit_check),
  n = 36
)


write.csv(
  deficit_check,
  file.path(
    tables_dir,
    "AM_Fig3D_Matched_Control_Richness_Deficit_Data.csv"
  ),
  row.names = FALSE
)


################################################################################
## 8. ASSOCIATION WITH LOCAL APPLE NEMATODE ABUNDANCE
##
## Same partial-Spearman framework used in the original Fig. 4 analysis:
##
##   response  = richness deficit relative to matched control
##   predictor = local apple nematode abundance
##   covariates = ecological setup + collection day
##
################################################################################


deficit_covars <- deficit_data[
  ,
  c(
    "ES_Number",
    "Collection_day"
  )
]


## TOTAL ASV RICHNESS DEFICIT

total_deficit_res <- partial_spearman(
  deficit_data$Total_richness_deficit,
  deficit_data$log_apple_worms,
  deficit_covars
)


## RARE ASV RICHNESS DEFICIT

rare_deficit_res <- partial_spearman(
  deficit_data$Rare_richness_deficit,
  deficit_data$log_apple_worms,
  deficit_covars
)


## =============================================================================
## 9. RESULTS TABLE
## =============================================================================

deficit_results <- tibble(
  
  metric = c(
    "Total ASV richness deficit",
    "Rare ASV richness deficit"
  ),
  
  rho = c(
    as.numeric(
      total_deficit_res[
        "rho"
      ]
    ),
    as.numeric(
      rare_deficit_res[
        "rho"
      ]
    )
  ),
  
  p = c(
    as.numeric(
      total_deficit_res[
        "p"
      ]
    ),
    as.numeric(
      rare_deficit_res[
        "p"
      ]
    )
  ),
  
  n = c(
    as.numeric(
      total_deficit_res[
        "n"
      ]
    ),
    as.numeric(
      rare_deficit_res[
        "n"
      ]
    )
  )
)


cat(
  "\n============================================================\n",
  "RICHNESS DEFICIT RESULTS\n",
  "============================================================\n"
)

print(
  deficit_results
)


write.csv(
  deficit_results,
  file.path(
    tables_dir,
    "AM_Fig3D_Matched_Control_Richness_Deficit_Results.csv"
  ),
  row.names = FALSE
)


################################################################################
## 10. FIGURE
################################################################################

library(grid)


## -----------------------------------------------------------------------------
## Prepare plotting data
## -----------------------------------------------------------------------------

deficit_plot_data <- deficit_data |>
  select(
    SampleID,
    ES_Number,
    Collection_day,
    apple_worms,
    Total_richness_deficit,
    Rare_richness_deficit
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Total_richness_deficit,
      Rare_richness_deficit
    ),
    names_to = "metric",
    values_to = "richness_deficit"
  ) |>
  mutate(
    metric = recode(
      metric,
      Total_richness_deficit = "Total ASV richness",
      Rare_richness_deficit = "Rare ASV richness"
    ),
    metric = factor(
      metric,
      levels = c(
        "Total ASV richness",
        "Rare ASV richness"
      )
    )
  )


## -----------------------------------------------------------------------------
## Statistics labels + in-panel titles
##
## IMPORTANT:
## Greek rho is rendered using the SAME method as AM_Fig4:
##     as.expression(bquote(... rho ...))
## This avoids Unicode/font-encoding problems in PDF output.
## -----------------------------------------------------------------------------

metric_levels_d <- c(
  "Rare ASV richness",
  "Total ASV richness"
)

deficit_plot_data <- deficit_plot_data |>
  mutate(
    metric = factor(
      as.character(metric),
      levels = metric_levels_d
    )
  )


## Per-facet ranges for placing titles/statistics INSIDE the frame.
deficit_ranges <- deficit_plot_data |>
  group_by(metric) |>
  summarise(
    x_max = max(apple_worms + 1, na.rm = TRUE),
    y_min = min(richness_deficit, na.rm = TRUE),
    y_max = max(richness_deficit, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    y_range = y_max - y_min,
    x_pos   = x_max,
    x_pos   = x_max * 1.25,
    y_title = y_max - 0.04 * y_range,
    y_stats = y_max - 0.75 * y_range
  )


## In-panel titles
deficit_titles <- deficit_ranges |>
  transmute(
    metric,
    x = x_pos,
    y = y_title,
    label = as.character(metric)
  )


## Helper used for plotmath P-value labels
## Same formatting logic used in Fig. 4.
format_p <- function(p) {
  if (p < 0.001) {
    "< 0.001"
  } else {
    paste0("= ", sprintf("%.3f", p))
  }
}


## Extract statistics
rho_total <- deficit_results$rho[
  deficit_results$metric == "Total ASV richness deficit"
]

p_total <- deficit_results$p[
  deficit_results$metric == "Total ASV richness deficit"
]

rho_rare <- deficit_results$rho[
  deficit_results$metric == "Rare ASV richness deficit"
]

p_rare <- deficit_results$p[
  deficit_results$metric == "Rare ASV richness deficit"
]


## EXACT SAME rho-rendering strategy used in AM_Fig4
label_total <- as.expression(
  bquote(
    atop(
      "Partial" ~ rho == .(sprintf("%+.3f", rho_total)),
      italic(P) ~ .(format_p(p_total))
    )
  )
)

label_rare <- as.expression(
  bquote(
    atop(
      "Partial" ~ rho == .(sprintf("%+.3f", rho_rare)),
      italic(P) ~ .(format_p(p_rare))
    )
  )
)


## Annotation positions, one row per facet
ann_total <- deficit_ranges |>
  filter(metric == "Total ASV richness") |>
  transmute(
    metric,
    x = x_pos,
    y = y_stats
  )

ann_rare <- deficit_ranges |>
  filter(metric == "Rare ASV richness") |>
  transmute(
    metric,
    x = x_pos,
    y = y_stats
  )


## -----------------------------------------------------------------------------
## Plot
## -----------------------------------------------------------------------------

p_deficit <- ggplot(
  deficit_plot_data,
  aes(
    x = apple_worms + 1,
    y = richness_deficit,
    colour = ES_Number
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.7,
    colour = "grey40"
  ) +
  
  geom_point(
    size = 4.0,
    alpha = 0.82
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "black",
    linewidth = 1.0,
    linetype = 2
  ) +
  
  ## In-panel facet titles
  geom_text(
    data = deficit_titles,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    family = "sans",
    fontface = "bold",
    size = 5.0,
    colour = "black"
  ) +
  
  ## Rare-richness statistics
  geom_text(
    data = ann_rare,
    aes(
      x = x,
      y = y
    ),
    inherit.aes = FALSE,
    label = label_rare,
    hjust = 1,
    vjust = 1,
    family = "sans",
    size = 4.6,
    lineheight = 0.95,
    colour = "black"
  ) +
  
  ## Total-richness statistics
  geom_text(
    data = ann_total,
    aes(
      x = x,
      y = y
    ),
    inherit.aes = FALSE,
    label = label_total,
    hjust = 1,
    vjust = 1,
    family = "sans",
    size = 4.6,
    lineheight = 0.95,
    colour = "black"
  ) +
  
  facet_wrap(
    ~ metric,
    nrow = 1,
    scales = "free_y",
    drop = TRUE
  ) +
  
  scale_x_log10(
    labels = scales::label_number(),
    expand = expansion(mult = c(0.03, 0.22))
  ) +
  
  scale_colour_manual(
    values = ES_COL,
    name = NULL
  ) +
  
  labs(
    x = "Nematodes recovered from apple\n(+1, log scale)",
    y = "Richness deficit relative to matched control"
  ) +
  
  theme_sfnh(
    base_size = 20
  ) +
  
  theme(
    ## No external facet-strip boxes
    strip.background = element_blank(),
    strip.text = element_blank(),
    
    ## Main Fig. 3 already has the ES legend in Panel C,
    ## so suppress the duplicate legend in D.
    legend.position = "none",
    
    axis.title.x = element_text(
      margin = margin(t = 5)
    ),
    
    plot.margin = margin(8, 10, 6, 4)
  )


p3d_base <- p_deficit
p3d <- p3d_base


## Standalone version keeps an ES legend for inspection/export
p3d_standalone <- p_deficit +
  theme(
    legend.position = "top",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )


p3d_standalone




################################################################################
################################################################################
##
## COMBINE + EXPORT MAIN AM_Fig3 (A, B, C, D)
##
## A + B on the first row, C on the second row, D on the third row.
## C and D each contain two internal facets, so stacking them preserves
## readability for the first render. Layout can be adjusted after inspection.
##
################################################################################
################################################################################

top_row_fig3 <- (
  patchwork::free(
    p3a,
    side = "b"
  ) |
    p3b
)

fig3 <- (
  top_row_fig3 /
    patchwork::free(
      p3c,
      side = "l"
    ) /
    patchwork::free(
      p3d,
      side = "l"
    )
) +
  plot_layout(
    heights = c(1, 1.05, 1.05)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 24
    )
  )

fig3

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig3_Richness_Response_MAIN_ABCD.pdf"
  ),
  plot = fig3,
  width = 340,
  height = 455,
  units = "mm"
)

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig3_Richness_Response_MAIN_ABCD.png"
  ),
  plot = fig3,
  width = 340,
  height = 455,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(figures_dir, "AM_Fig3A_Abundance_Classification.pdf"),
  p3a, width = 180, height = 140, units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_Fig3B_Effect_Sizes_BH.pdf"),
  p3b, width = 180, height = 140, units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_Fig3C_GLMM_Dose_Response.pdf"),
  p3c, width = 260, height = 140, units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_Fig3D_Matched_Control_Richness_Deficit.pdf"),
  p3d_standalone, width = 260, height = 150, units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_Fig3D_Matched_Control_Richness_Deficit.png"),
  p3d_standalone, width = 260, height = 150, units = "mm", dpi = 600, bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(figures_dir, "AM_Fig3_ABCD_sessionInfo.txt")
)

cat(
  "\n============================================================\n",
  "AM_Fig3 (A–D) COMPLETE\n",
  "============================================================\n\n",
  "All output saved under:\n",
  output_dir,
  "\n"
)

