################################################################################
################################################################################
## AM_FigS4
## FUNGAL ITS NON-RESPONSE AS A SPECIFICITY CONTROL
################################################################################
################################################################################


## =============================================================================
## 0. PACKAGES
## =============================================================================

library(tidyverse)
library(vegan)
library(patchwork)
library(scales)
library(lme4)
library(DHARMa)


## =============================================================================
## 1. PATHS
## =============================================================================

## Run this script with the working directory set to the repository root
## (open the repo as an RStudio project, or run `Rscript scripts/<name>.R`
## from the repo root). Place the data files in ./data ; outputs are
## written to ./output/Figures and ./output/Tables.
input_dir <- "data"

output_dir <- "output"

figures_dir <- file.path(output_dir, "Figures")
tables_dir  <- file.path(output_dir, "Tables")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)


## =============================================================================
## 2. INPUT FILES
## =============================================================================

its_metadata_file <- file.path(input_dir, "its_metadata_apple.csv")
its_counts_file   <- file.path(input_dir, "its_counts_rarefied14000_apple.csv")


## =============================================================================
## 3. READ DATA
## =============================================================================

its_meta <- read.csv(its_metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

fungal_cnt_ref <- read.csv(
  its_counts_file,
  row.names = 1,
  check.names = FALSE
) |>
  as.matrix()

storage.mode(fungal_cnt_ref) <- "numeric"

cat("\nFungal ITS table dimensions:", dim(fungal_cnt_ref)[1], "ASVs x",
    dim(fungal_cnt_ref)[2], "samples\n")

## Confirm the table is rarefied to a single common depth (14,000, per the
## filename and independently verified: min = max = 14,000 across all 37
## samples).
depth_check <- colSums(fungal_cnt_ref)
cat("Fungal sample-total range after rarefaction (should be constant):",
    range(depth_check), "\n")
stopifnot(
  "Fungal count table does not appear to be rarefied to a single depth" =
    diff(range(depth_check)) == 0
)


## =============================================================================
## 4. FORMAT METADATA
## =============================================================================

its_meta <- its_meta |>
  mutate(
    ES_Number = factor(ES_Number, levels = c("ES1", "ES2", "ES3")),
    Collection_day = factor(Collection_day, levels = c("12D", "15D", "18D", "21D")),
    cell_id = interaction(ES_Number, Collection_day, drop = TRUE)
  )


## =============================================================================
## 5. APPLE / NEMATODE SAMPLES — fungal ITS worm-inoculated apples
##
## its_metadata_apple.csv is self-contained (no join to the bacterial
## metadata.csv needed or possible — SampleID conventions differ). Nematode
## abundance is read directly from Ave_worm_pot.
## =============================================================================

WF <- its_meta |>
  filter(
    Treatment == "Worm",
    SampleID %in% colnames(fungal_cnt_ref)
  ) |>
  arrange(ES_Number, Collection_day, Pot_Replicate) |>
  mutate(
    apple_worms = Ave_worm_pot
  )

cnt_apple_fungal <- fungal_cnt_ref[, WF$SampleID, drop = FALSE]

stopifnot(identical(WF$SampleID, colnames(cnt_apple_fungal)))

cat("\nFungal apple/nematode samples matched:", nrow(WF), "\n")
cat("Cells represented:", n_distinct(WF$cell_id), "of 12\n")

## Confirmed against the actual data: n = 28 (not 36), 12 cells represented
## but unbalanced (1-3 worm-inoculated pots per cell; ES3/21D has only 1).
## This is expected — flag loudly if it ever changes (e.g. after adding
## more ITS libraries), since the interpretation below assumes this n.
if (nrow(WF) != 28) {
  warning(
    "Fungal worm-inoculated sample count (", nrow(WF), ") differs from the ",
    "n = 28 confirmed when this script was written. If ITS libraries were ",
    "added or removed, re-check cell balance and update the Methods/caption ",
    "n accordingly."
  )
}


## =============================================================================
## 6. DIVERSITY / RICHNESS METRICS (fungal ITS)
## =============================================================================

WF$Fungal_Shannon        <- vegan::diversity(t(cnt_apple_fungal), index = "shannon")
WF$Fungal_Total_richness <- colSums(cnt_apple_fungal > 0)

## Pre-transform the x-variable by hand rather than relying on scale_x_log10().
## Reason: annotate(x = -Inf, ...) combined with scale_x_log10() sends -Inf
## through log10() BEFORE ggplot's Inf-handling special case applies, which
## returns NaN (not -Inf) and silently drops the label row (confirmed via the
## "In log(x, base) : NaNs produced" / "Removed 1 row ... (geom_label())"
## warnings). Plotting on an already-log10-transformed column with a LINEAR
## scale avoids this entirely, since -Inf/Inf are untouched on a linear scale.
WF$log_worms <- log10(WF$apple_worms + 1)

worm_axis_breaks <- scales::breaks_log()(range(WF$apple_worms + 1, na.rm = TRUE))
worm_axis_breaks <- worm_axis_breaks[worm_axis_breaks > 0]
worm_axis_labels <- scales::label_number()(worm_axis_breaks)
worm_axis_pos    <- log10(worm_axis_breaks)


## =============================================================================
## 7. PARTIAL SPEARMAN FUNCTION (identical to the bacterial pipeline)
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

cov_ed_fungal   <- WF[, c("ES_Number", "Collection_day")]
worm_predictor_fungal <- log1p(WF$apple_worms)


## =============================================================================
## 8. COLORS + THEME — matched to Figures 3 and 4
## =============================================================================

ES_COL <- c("ES1" = "#264653", "ES2" = "#2a9d8f", "ES3" = "#e76f51")

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
## PANEL A — FUNGAL SHANNON DIVERSITY vs NEMATODE ABUNDANCE
## PANEL B — FUNGAL TOTAL RICHNESS vs NEMATODE ABUNDANCE
##
## Same partial-Spearman test and figure style as the bacterial Fig. 3
## panels, applied to the fungal ITS metrics. The expected (specificity-
## control) result is a NON-significant association for both.
##
################################################################################
################################################################################

fungal_shannon_res <- partial_spearman(WF$Fungal_Shannon, worm_predictor_fungal, cov_ed_fungal)
fungal_total_res   <- partial_spearman(WF$Fungal_Total_richness, worm_predictor_fungal, cov_ed_fungal)

fungal_results <- tibble(
  metric = c("Fungal Shannon diversity", "Fungal ASV richness"),
  rho = c(as.numeric(fungal_shannon_res["rho"]), as.numeric(fungal_total_res["rho"])),
  p   = c(as.numeric(fungal_shannon_res["p"]),   as.numeric(fungal_total_res["p"])),
  n   = c(as.numeric(fungal_shannon_res["n"]),   as.numeric(fungal_total_res["n"]))
)

cat("\n============================================================\n",
    "AM_FigS4 — FUNGAL ITS SPECIFICITY CONTROL\n",
    "============================================================\n")
print(fungal_results)

write.csv(
  fungal_results,
  file.path(tables_dir, "AM_FigS4_Fungal_Specificity_Control.csv"),
  row.names = FALSE
)

make_label <- function(res) {
  
  p_part <- if (res["p"] < 0.001) {
    "italic(P) < 0.001"
  } else {
    paste0("italic(P) == ", sprintf("%.3f", res["p"]))
  }
  
  paste0(
    "rho == ", sprintf("%+.3f", res["rho"]),
    " * \";\" ~~ ", p_part
  )
}

# Statistics are placed at the upper-right of each panel.
# Using Inf keeps the placement stable if the data range changes.
stats_x <- Inf
stats_y_A <- Inf
stats_y_B <- Inf


AM_FigS4A <- ggplot(WF, aes(x = log_worms, y = Fungal_Shannon, colour = ES_Number)) +
  geom_point(size = 4.2, alpha = 0.82) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1.0, linetype = 2) +
  annotate(
    "text",
    x = stats_x,
    y = stats_y_A,
    hjust = 1.05,
    vjust = 1.10,
    label = make_label(fungal_shannon_res),
    parse = TRUE,
    family = "sans",
    size = 4,
    colour = "black",
    lineheight = 0.95
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(
    breaks = worm_axis_pos, labels = worm_axis_labels,
    expand = expansion(mult = c(0.05, 0.28))
  ) +
  scale_colour_manual(values = ES_COL, name = NULL) +
  labs(x = "Nematodes recovered from apple\n(+1, log scale)", y = "Fungal Shannon diversity") +
  theme_sfnh(base_size = 18) +
  theme(
    legend.position = "none"
  )

AM_FigS4B <- ggplot(WF, aes(x = log_worms, y = Fungal_Total_richness, colour = ES_Number)) +
  geom_point(size = 4.2, alpha = 0.82) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1.0, linetype = 2) +
  annotate(
    "text",
    x = stats_x,
    y = stats_y_B,
    hjust = 1.05,
    vjust = 1.10,
    label = make_label(fungal_total_res),
    parse =TRUE,
    family = "sans",
    size = 4,
    colour = "black",
    lineheight = 0.95
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(
    breaks = worm_axis_pos, labels = worm_axis_labels,
    expand = expansion(mult = c(0.05, 0.28))
  ) +
  scale_colour_manual(values = ES_COL, name = NULL) +
  labs(x = "Nematodes recovered from apple\n(+1, log scale)", y = "Fungal ASV richness") +
  theme_sfnh(base_size = 18) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 8)
  )


################################################################################
## COMBINE + EXPORT
################################################################################

AM_FigS4 <- (AM_FigS4A | AM_FigS4B) +
  plot_layout(widths = c(1, 1.18)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "sans", face = "bold", size = 22))

AM_FigS4

ggsave(
  file.path(figures_dir, "AM_FigS4_Fungal_Specificity_Control.pdf"),
  AM_FigS4, width = 290, height = 130, units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_FigS4_Fungal_Specificity_Control.png"),
  AM_FigS4, width = 290, height = 130, units = "mm", dpi = 600, bg = "white"
)

cat("\n============================================================\n",
    "AM_FigS4 COMPLETE\n",
    "============================================================\n\n",
    "All output saved under:\n", output_dir, "\n")

