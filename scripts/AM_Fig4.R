################################################################################
# AM_Fig4
# Apple bacterial community displacement and dispersion
#
# A  Bray-Curtis distance to matched uninoculated control
# B  Between-cell replicate dispersion
# C  Residual within-cell replicate dispersion
#
# Each cell = ecological setup x collection day
#           = 3 nematode-inoculated pots + 1 matched uninoculated control
################################################################################


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(patchwork)
  library(scales)
})


# =============================================================================
# PATHS
# =============================================================================

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


# ==============================================================================
# 2. READ DATA
# ==============================================================================

meta <- read.csv(
  file.path(input_dir, "metadata.csv"),
  stringsAsFactors = FALSE
)

cnt <- read.csv(
  file.path(input_dir, "bac_counts_rarefied14000.csv"),
  row.names = 1,
  check.names = FALSE
) |>
  as.matrix()

storage.mode(cnt) <- "numeric"


# ==============================================================================
# 3. FORMAT METADATA
# ==============================================================================

meta$ES_Number <- factor(
  meta$ES_Number,
  levels = c("ES1", "ES2", "ES3")
)

meta$Collection_day <- factor(
  meta$Collection_day,
  levels = c("12D", "15D", "18D", "21D")
)


# ==============================================================================
# 4. APPLE SAMPLES AND BRAY-CURTIS DISTANCES
# ==============================================================================

apple_meta <- meta |>
  filter(
    Compartment == "As",
    SampleID %in% colnames(cnt)
  )

apple_ids <- apple_meta$SampleID

cnt_apple <- cnt[
  ,
  apple_ids,
  drop = FALSE
]

rel_apple <- sweep(
  cnt_apple,
  2,
  colSums(cnt_apple),
  "/"
) * 100

D_all <- as.matrix(
  vegdist(
    t(rel_apple),
    method = "bray"
  )
)


# ==============================================================================
# 5. WORM AND CONTROL METADATA
# ==============================================================================

W <- apple_meta |>
  filter(Treatment == "Worm") |>
  mutate(
    apple_worms = worm,
    worm_log1p = log1p(worm),
    cell = interaction(
      ES_Number,
      Collection_day,
      drop = TRUE
    )
  )

C <- apple_meta |>
  filter(Treatment == "Control") |>
  mutate(
    cell = interaction(
      ES_Number,
      Collection_day,
      drop = TRUE
    )
  )

cat("\nApple/nematode-inoculated samples:", nrow(W), "\n")
cat("Apple/Control samples:", nrow(C), "\n")
cat("ES x day cells:", n_distinct(W$cell), "\n\n")


# ==============================================================================
# 6. SHARED FUNCTIONS
# ==============================================================================

partial_spearman <- function(y, x, covars) {
  
  keep <- complete.cases(y, x, covars)
  
  y <- y[keep]
  x <- x[keep]
  covars <- covars[keep, , drop = FALSE]
  
  mm <- model.matrix(
    ~ .,
    data = as.data.frame(covars)
  )
  
  ry <- residuals(
    lm(rank(y) ~ mm - 1)
  )
  
  rx <- residuals(
    lm(rank(x) ~ mm - 1)
  )
  
  ct <- cor.test(
    ry,
    rx,
    method = "pearson"
  )
  
  c(
    rho = unname(ct$estimate),
    p = ct$p.value,
    n = length(y)
  )
}


center_within <- function(x, group) {
  
  r <- rank(
    x,
    na.last = "keep"
  )
  
  r - ave(
    r,
    group,
    FUN = function(z) {
      mean(z, na.rm = TRUE)
    }
  )
}


perm_within <- function(
    y,
    x,
    strata,
    nperm = 9999,
    seed = 17
) {
  
  set.seed(seed)
  
  obs <- abs(
    cor(
      y,
      x,
      method = "pearson",
      use = "complete.obs"
    )
  )
  
  perm_r <- numeric(nperm)
  
  for (b in seq_len(nperm)) {
    
    xp <- x
    
    for (g in unique(strata)) {
      
      idx <- which(strata == g)
      
      xp[idx] <- sample(
        x[idx]
      )
    }
    
    perm_r[b] <- abs(
      cor(
        y,
        xp,
        method = "pearson",
        use = "complete.obs"
      )
    )
  }
  
  (
    sum(perm_r >= obs) + 1
  ) / (
    nperm + 1
  )
}


format_p <- function(p) {
  
  if (p < 0.001) {
    
    "< 0.001"
    
  } else {
    
    paste0(
      "= ",
      sprintf("%.3f", p)
    )
  }
}


# ==============================================================================
# 7. COLORS AND THEME
#    Colors retained from the running old Figure 6 script.
# ==============================================================================

ES_COL <- c(
  ES1 = "#264653",
  ES2 = "#2a9d8f",
  ES3 = "#e76f51"
)

theme_sfnh <- function(base_size = 18) {
  
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    
    theme(
      text = element_text(
        family = "sans",
        colour = "black"
      ),
      
      axis.title = element_text(
        family = "sans",
        colour = "black"
      ),
      
      axis.text = element_text(
        family = "sans",
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      legend.title = element_blank(),
      
      # No legend box.
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      
      plot.margin = margin(
        12,
        15,
        12,
        15
      )
    )
}


# ==============================================================================
# 8. AM_Fig4A
#    BRAY-CURTIS DISTANCE TO MATCHED CONTROL
#    Old Figure 6A
# ==============================================================================

dist_ctrl <- map_dfr(
  W$SampleID,
  function(w) {
    
    es <- W$ES_Number[
      W$SampleID == w
    ]
    
    day <- W$Collection_day[
      W$SampleID == w
    ]
    
    ctrl <- C$SampleID[
      C$ES_Number == es &
        C$Collection_day == day
    ]
    
    if (length(ctrl) != 1) {
      return(NULL)
    }
    
    tibble(
      SampleID = w,
      ES_Number = es,
      Collection_day = day,
      apple_worms = W$apple_worms[
        W$SampleID == w
      ],
      worm_log1p = W$worm_log1p[
        W$SampleID == w
      ],
      dist_ctrl = D_all[
        ctrl,
        w
      ]
    )
  }
)

res_A <- partial_spearman(
  y = dist_ctrl$dist_ctrl,
  x = dist_ctrl$worm_log1p,
  covars = dist_ctrl[
    ,
    c(
      "ES_Number",
      "Collection_day"
    )
  ]
)

rho_A <- unname(res_A["rho"])
P_A   <- unname(res_A["p"])

# Panel A
label_A <- as.expression(
  bquote(
    atop(
      "Partial" ~ rho == .(sprintf("%+.3f", rho_A)),
      italic(P) ~ .(format_p(P_A))
    )
  )
)

p4a <- ggplot(
  dist_ctrl,
  aes(
    x = apple_worms + 1,
    y = dist_ctrl,
    colour = ES_Number
  )
) +
  
  geom_point(
    size = 3.2,
    alpha = 0.85
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 0.7,
    linetype = 2
  ) +
  
  scale_x_log10(
    breaks = c(
      1,
      10,
      100,
      1000,
      10000
    ),
    labels = c(
      "0",
      "10",
      "100",
      "1,000",
      "10,000"
    ),
    expand = expansion(
      mult = c(0.03, 0.20)
    )
  ) +
  
  scale_colour_manual(
    values = ES_COL,
    name = NULL
  ) +
  
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = label_A,
    hjust = 1.05,
    vjust = 1.15,
    family = "sans",
    size = 4.0,
    linewidth = 0,
    label.padding = unit(
      0.20,
      "lines"
    ),
    fill = "white"
  ) +
  
  labs(
    x = "Nematodes recovered from apple",
    y = "Bray-Curtis distance\nto matched control"
  ) +
  
  theme_sfnh() +
  
  theme(
    legend.position = "none")


p4a
# ==============================================================================
# 9. DISPERSION DATA
#    Mean Bray-Curtis distance from each nematode-inoculated pot to the other two nematode-inoculated pots
#    within the same ES x collection-day cell.
# ==============================================================================

disp <- map_dfr(
  W$SampleID,
  function(w) {
    
    es <- W$ES_Number[
      W$SampleID == w
    ]
    
    day <- W$Collection_day[
      W$SampleID == w
    ]
    
    others <- W$SampleID[
      W$ES_Number == es &
        W$Collection_day == day &
        W$SampleID != w
    ]
    
    if (length(others) == 0) {
      return(NULL)
    }
    
    tibble(
      SampleID = w,
      ES_Number = es,
      Collection_day = day,
      
      cell = interaction(
        es,
        day,
        drop = TRUE
      ),
      
      apple_worms = W$apple_worms[
        W$SampleID == w
      ],
      
      mean_pair = mean(
        D_all[
          w,
          others
        ]
      )
    )
  }
)


# ==============================================================================
# 10. AM_Fig4B
#     BETWEEN-CELL DISPERSION
#     Old Figure 6D
# ==============================================================================

cell_disp <- disp |>
  group_by(
    ES_Number,
    Collection_day,
    cell
  ) |>
  summarise(
    cell_dispersion = mean(
      mean_pair
    ),
    median_worms = median(
      apple_worms
    ),
    .groups = "drop"
  )

res_B <- cor.test(
  cell_disp$median_worms,
  cell_disp$cell_dispersion,
  method = "spearman",
  exact = FALSE
)

rho_B <- unname(
  res_B$estimate
)

P_B <- res_B$p.value

label_B <- as.expression(
  bquote(
    atop(
      rho == .(sprintf("%+.3f", rho_B)),
      italic(P) ~ .(format_p(P_B))
    )
  )
)

p4b <- ggplot(
  cell_disp,
  aes(
    x = median_worms + 1,
    y = cell_dispersion,
    colour = ES_Number
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 0.7,
    linetype = 2
  ) +
  
  scale_x_log10(
    breaks = c(
      1,
      10,
      100,
      1000,
      10000
    ),
    labels = c(
      "0",
      "10",
      "100",
      "1,000",
      "10,000"
    ),
    expand = expansion(
      mult = c(0.03, 0.20)
    )
  ) +
  
  scale_colour_manual(
    values = ES_COL,
    name = NULL
  ) +
  
  annotate(
    "label",
    x = Inf,
    y = -Inf,
    label = label_B,
    hjust = 1.05,
    vjust = -0.15,
    family = "sans",
    size = 4.0,
    linewidth = 0,
    label.padding = unit(
      0.20,
      "lines"
    ),
    fill = "white"
  ) +
  
  labs(
    x = "Median nematodes recovered from apple",
    y = "Mean within-cell\nBray-Curtis dispersion"
  ) +
  
  theme_sfnh() +
  
  theme(
    legend.position = "none"
  )

p4b
# ==============================================================================
# 11. AM_Fig4C
#     RESIDUAL WITHIN-CELL DISPERSION
#
#     Replaces the old raw within-cell dispersion panel.
#
#     Rationale:
#     The old Panel C could become positive simply because samples with
#     different nematode abundances occupy different positions along one common
#     nematode-associated compositional gradient. Here we first remove that
#     fitted directional association from community composition, then ask
#     whether residual community dispersion still increases with nematode
#     abundance within ES x collection-day cells.
#
#     Residualization is performed on the first 6 Bray-Curtis PCoA axes.
#     Each axis is fitted exactly as in the existing residual-dispersion script:
#         PCoA axis ~ ES_Number + Collection_day + log1p(local nematode abundance)
#     The residual distances are then evaluated using the same within-cell
#     rank-centering and restricted permutation test used previously.
# ==============================================================================

# ------------------------------------------------------------------------------
# 11A. Reproduce the OLD raw within-cell result for comparison only.
#      This is NOT plotted in the revised AM_Fig4.
# ------------------------------------------------------------------------------

disp_raw_within <- disp |>
  mutate(
    worm_rank_centered = center_within(
      apple_worms,
      cell
    ),
    dispersion_rank_centered = center_within(
      mean_pair,
      cell
    )
  )

r_C_raw <- cor(
  disp_raw_within$worm_rank_centered,
  disp_raw_within$dispersion_rank_centered,
  method = "pearson",
  use = "complete.obs"
)

P_C_raw <- perm_within(
  y = disp_raw_within$dispersion_rank_centered,
  x = disp_raw_within$worm_rank_centered,
  strata = disp_raw_within$cell,
  nperm = 9999,
  seed = 17
)


# ------------------------------------------------------------------------------
# 11B. Bray-Curtis PCoA for the 36 nematode-inoculated apple samples.
# ------------------------------------------------------------------------------

worm_ids <- as.character(W$SampleID)

D_worm <- as.dist(
  D_all[
    worm_ids,
    worm_ids,
    drop = FALSE
  ]
)

pcoa_resid <- cmdscale(
  D_worm,
  k = 6,
  eig = TRUE
)

axis_scores <- as.data.frame(
  pcoa_resid$points
)

colnames(axis_scores) <- paste0(
  "PCoA",
  seq_len(ncol(axis_scores))
)

axis_scores$SampleID <- rownames(axis_scores)

resid_dat <- axis_scores |>
  as_tibble() |>
  left_join(
    W |>
      select(
        SampleID,
        ES_Number,
        Collection_day,
        cell,
        apple_worms,
        worm_log1p
      ),
    by = "SampleID"
  )

stopifnot(nrow(resid_dat) == 36)
stopifnot(!anyNA(resid_dat$apple_worms))


# ------------------------------------------------------------------------------
# 11C. Remove the fitted directional nematode-abundance effect from each axis.
#      This reproduces the existing residual-dispersion test exactly.
# ------------------------------------------------------------------------------

n_axes <- 6

resid_mat <- matrix(
  NA_real_,
  nrow = nrow(resid_dat),
  ncol = n_axes,
  dimnames = list(
    resid_dat$SampleID,
    paste0("Resid", seq_len(n_axes))
  )
)

for (k in seq_len(n_axes)) {
  form_k <- as.formula(
    paste0(
      "PCoA",
      k,
      " ~ ES_Number + Collection_day + worm_log1p"
    )
  )
  
  fit_k <- lm(
    form_k,
    data = resid_dat
  )
  
  resid_mat[, k] <- residuals(fit_k)
}


# ------------------------------------------------------------------------------
# 11D. Euclidean distances in residual composition space.
# ------------------------------------------------------------------------------

D_resid <- as.matrix(
  dist(
    resid_mat,
    method = "euclidean"
  )
)


# ------------------------------------------------------------------------------
# 11E. Residual within-cell dispersion for each pot.
#      For each pot: mean residual distance to its two cell-mates.
# ------------------------------------------------------------------------------

disp_resid <- map_dfr(
  seq_len(nrow(resid_dat)),
  function(i) {
    
    this_id <- resid_dat$SampleID[i]
    this_cell <- resid_dat$cell[i]
    
    same_cell_idx <- which(
      resid_dat$cell == this_cell &
        seq_len(nrow(resid_dat)) != i
    )
    
    if (length(same_cell_idx) != 2) {
      stop(
        "Expected exactly 2 nematode-inoculated cell-mates for sample ",
        this_id,
        ", but found ",
        length(same_cell_idx),
        "."
      )
    }
    
    tibble(
      SampleID = this_id,
      ES_Number = resid_dat$ES_Number[i],
      Collection_day = resid_dat$Collection_day[i],
      cell = this_cell,
      apple_worms = resid_dat$apple_worms[i],
      residual_dispersion = mean(
        D_resid[
          this_id,
          resid_dat$SampleID[same_cell_idx]
        ]
      )
    )
  }
)


# ------------------------------------------------------------------------------
# 11F. Same within-cell rank-centering and restricted permutation test used for
#      the old Panel C, now applied to RESIDUAL community dispersion.
# ------------------------------------------------------------------------------

disp_resid_within <- disp_resid |>
  mutate(
    worm_rank_centered = center_within(
      apple_worms,
      cell
    ),
    residual_dispersion_rank_centered = center_within(
      residual_dispersion,
      cell
    )
  )

r_C <- cor(
  disp_resid_within$worm_rank_centered,
  disp_resid_within$residual_dispersion_rank_centered,
  method = "pearson",
  use = "complete.obs"
)

P_C <- perm_within(
  y = disp_resid_within$residual_dispersion_rank_centered,
  x = disp_resid_within$worm_rank_centered,
  strata = disp_resid_within$cell,
  nperm = 9999,
  seed = 17
)

label_C <- as.expression(
  bquote(
    atop(
      "Residual within-cell" ~ italic(r) == .(sprintf("%+.3f", r_C)),
      "Permutation" ~ italic(P) ~ .(format_p(P_C))
    )
  )
)


# ------------------------------------------------------------------------------
# 11G. Revised Panel C.
# ------------------------------------------------------------------------------

p4c <- ggplot(
  disp_resid_within,
  aes(
    x = worm_rank_centered,
    y = residual_dispersion_rank_centered,
    colour = ES_Number
  )
) +
  
  geom_point(
    size = 3.2,
    alpha = 0.85
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 0.7,
    linetype = 2
  ) +
  
  scale_colour_manual(
    values = ES_COL,
    name = NULL
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = label_C,
    hjust = 1.05,
    vjust = -0.15,
    family = "sans",
    size = 4.0
  )+
  
  labs(
    x = "Nematode abundance\n(within-cell rank-centered)",
    y = "Residual community dispersion\n(within-cell rank-centered)"
  ) +
  
  theme_sfnh() +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = 13),
    
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    
    legend.box.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 8
    ),
    
    legend.margin = margin(0, 0, 0, 0)
  )

p4c

# ==============================================================================
# 12. PRINT STATISTICAL RESULTS
# ==============================================================================

cat(
  "\n============================================================\n",
  "AM_Fig4 RESULTS\n",
  "============================================================\n\n"
)

cat(
  "Panel A - distance to matched control\n"
)

cat(
  sprintf(
    "rho = %+.3f, P = %.5f, n = %d\n\n",
    rho_A,
    P_A,
    nrow(dist_ctrl)
  )
)

cat(
  "Panel B - between-cell dispersion\n"
)

cat(
  sprintf(
    "rho = %+.3f, P = %.5f, n = %d cells\n\n",
    rho_B,
    P_B,
    nrow(cell_disp)
  )
)

cat(
  "Old raw Panel C - within-cell dispersion (diagnostic only)\n"
)

cat(
  sprintf(
    "r = %+.3f, permutation P = %.5f, n = %d pots\n\n",
    r_C_raw,
    P_C_raw,
    nrow(disp_raw_within)
  )
)

cat(
  "Panel C - residual within-cell dispersion\n"
)

cat(
  sprintf(
    "r = %+.3f, permutation P = %.5f, n = %d pots\n\n",
    r_C,
    P_C,
    nrow(disp_resid_within)
  )
)

if (P_C < 0.05 && sign(r_C) == sign(r_C_raw)) {
  cat(
    "Residual dispersion remains positively associated with nematode abundance\n",
    "after removing the fitted directional abundance-associated shift.\n",
    "This supports dispersion beyond a shared directional gradient, but does\n",
    "not by itself establish stochastic community assembly.\n\n"
  )
} else {
  cat(
    "The residual dispersion association weakens, disappears, or changes sign\n",
    "after removing the fitted directional abundance-associated shift.\n",
    "The raw dispersion pattern should therefore not be interpreted as evidence\n",
    "of divergence beyond a shared directional gradient.\n\n"
  )
}


# ==============================================================================
# 13. SAVE TABLES
# ==============================================================================

write.csv(
  dist_ctrl,
  file.path(
    tables_dir,
    "AM_Fig4A_Distance_to_Control.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_disp,
  file.path(
    tables_dir,
    "AM_Fig4B_Between_Cell_Dispersion.csv"
  ),
  row.names = FALSE
)

write.csv(
  disp_resid_within,
  file.path(
    tables_dir,
    "AM_Fig4C_Residual_Within_Cell_Dispersion.csv"
  ),
  row.names = FALSE
)

# Save residual-axis values for reproducibility.
resid_axis_out <- cbind(
  resid_dat |>
    select(
      SampleID,
      ES_Number,
      Collection_day,
      cell,
      apple_worms,
      worm_log1p
    ),
  as.data.frame(resid_mat)
)

write.csv(
  resid_axis_out,
  file.path(
    tables_dir,
    "AM_Fig4C_Residual_PCoA_Axes.csv"
  ),
  row.names = FALSE
)

write.csv(
  tibble(
    analysis = c(
      "old_raw_within_cell_diagnostic",
      "residual_within_cell_main"
    ),
    r = c(
      r_C_raw,
      r_C
    ),
    permutation_P = c(
      P_C_raw,
      P_C
    ),
    n = c(
      nrow(disp_raw_within),
      nrow(disp_resid_within)
    )
  ),
  file.path(
    tables_dir,
    "AM_Fig4C_Raw_vs_Residual_Statistics.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. COMBINE AM_Fig4
#
# Three retained panels only.
# A and B are retained; C is the new residual-dispersion analysis.
#
# Do NOT collect guides, because the legends need to stay inside the data frame.
# ==============================================================================

fig4 <- (
  p4a |
    p4b |
    p4c
) +
  plot_layout(
    widths = c(
      1,
      1,
      1
    )
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

fig4


# ==============================================================================
# 15. EXPORT
# ==============================================================================

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig4_Community_Dispersion_MAIN_revised.pdf"
  ),
  plot = fig4,
  width = 390,
  height = 135,
  units = "mm"
)

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig4_Community_Dispersion_MAIN_revised.png"
  ),
  plot = fig4,
  width = 390,
  height = 135,
  units = "mm",
  dpi = 600,
  bg = "white"
)


# Individual panels

ggsave(
  file.path(
    figures_dir,
    "AM_Fig4A_Distance_to_Control.pdf"
  ),
  p4a,
  width = 150,
  height = 120,
  units = "mm"
)

ggsave(
  file.path(
    figures_dir,
    "AM_Fig4B_Between_Cell_Dispersion.pdf"
  ),
  p4b,
  width = 150,
  height = 120,
  units = "mm"
)

ggsave(
  file.path(
    figures_dir,
    "AM_Fig4C_Residual_Within_Cell_Dispersion.pdf"
  ),
  p4c,
  width = 150,
  height = 120,
  units = "mm"
)


# ==============================================================================
# 16. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    tables_dir,
    "AM_Fig4_sessionInfo.txt"
  )
)

cat(
  "\nAM_Fig4 complete.\n",
  "Saved under:\n",
  output_dir,
  "\n"
)

