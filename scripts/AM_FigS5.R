################################################################################
# AM_FigS5
# Directional alignment and high-vs-low convergence
#
# Extracted from former Figure 5 panels B and C.
#
# New panel labels:
#   A  Directional alignment of community displacement
#   B  High-vs-low nematode-abundance convergence
#
# Each cell = ecological setup x collection day
#           = 3 Worm pots + 1 matched Control pot
################################################################################


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(patchwork)
})


# ==============================================================================
# 1. PATHS
# ==============================================================================

## Run this script with the working directory set to the repository root
## (open the repo as an RStudio project, or run `Rscript scripts/<name>.R`
## from the repo root). Place the data files in ./data ; outputs are
## written to ./output/Figures and ./output/Tables.
input_dir <- "data"

output_dir <- "output"

figures_dir <- file.path(
  output_dir,
  "Figures"
)

tables_dir <- file.path(
  output_dir,
  "Tables"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. READ DATA
# ==============================================================================

meta <- read.csv(
  file.path(
    input_dir,
    "metadata.csv"
  ),
  stringsAsFactors = FALSE
)

cnt <- read.csv(
  file.path(
    input_dir,
    "bac_counts_rarefied14000.csv"
  ),
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
  levels = c(
    "ES1",
    "ES2",
    "ES3"
  )
)

meta$Collection_day <- factor(
  meta$Collection_day,
  levels = c(
    "12D",
    "15D",
    "18D",
    "21D"
  )
)

meta$Treatment <- factor(
  meta$Treatment,
  levels = c(
    "Control",
    "Worm"
  )
)


# ==============================================================================
# 4. PREPARE APPLE BACTERIAL DATA
# ==============================================================================

apple_meta <- meta %>%
  filter(
    Compartment == "As",
    SampleID %in% colnames(cnt)
  ) %>%
  arrange(
    match(
      SampleID,
      colnames(cnt)
    )
  )

apple_ids <- apple_meta$SampleID

cnt_apple <- cnt[
  ,
  apple_ids,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(cnt_apple),
    apple_meta$SampleID
  )
)

# Relative abundance (%)
rel_apple <- sweep(
  cnt_apple,
  2,
  colSums(cnt_apple),
  "/"
) * 100

# Bray-Curtis among all 48 apple samples:
# 36 Worm + 12 Control
D_all <- as.matrix(
  vegdist(
    t(rel_apple),
    method = "bray"
  )
)


# ==============================================================================
# 5. WORM AND CONTROL METADATA
# ==============================================================================

W <- apple_meta %>%
  filter(
    Treatment == "Worm"
  ) %>%
  mutate(
    apple_worms = worm,
    worm_log1p = log1p(worm),
    cell = interaction(
      ES_Number,
      Collection_day,
      drop = TRUE
    )
  )

C <- apple_meta %>%
  filter(
    Treatment == "Control"
  ) %>%
  mutate(
    cell = interaction(
      ES_Number,
      Collection_day,
      drop = TRUE
    )
  )

cat(
  "\nApple/Worm samples:",
  nrow(W),
  "\n"
)

cat(
  "Apple/Control samples:",
  nrow(C),
  "\n"
)

cat(
  "ES x collection-day cells:",
  n_distinct(W$cell),
  "\n\n"
)

stopifnot(
  nrow(W) == 36,
  nrow(C) == 12,
  n_distinct(W$cell) == 12
)


# ==============================================================================
# 6. SHARED HELPERS
# ==============================================================================

format_p <- function(p) {
  
  if (p < 0.001) {
    
    "< 0.001"
    
  } else {
    
    paste0(
      "= ",
      sprintf(
        "%.3f",
        p
      )
    )
  }
}


# ==============================================================================
# 7. PUBLICATION THEME
# ==============================================================================

theme_sfnh <- function(
    base_size = 18
) {
  
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    
    theme(
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      legend.position = "none",
      
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


################################################################################
################################################################################
#
# AM_FigS5A
# DIRECTIONAL ALIGNMENT OF COMMUNITY DISPLACEMENT
#
################################################################################
################################################################################


# ==============================================================================
# 8. SIX-AXIS PCoA
# ==============================================================================

# Same dimensionality used in the original analysis
K <- 6

pc <- cmdscale(
  as.dist(
    D_all
  ),
  k = K,
  eig = TRUE
)$points


# ==============================================================================
# 9. BUILD CONTROL -> HIGHEST-NEMATODE DISPLACEMENT VECTOR FOR EACH CELL
# ==============================================================================

cells <- W %>%
  distinct(
    ES_Number,
    Collection_day
  )

vecs <- list()

high_pots <- character()

low_pots <- character()

cell_ids <- character()


for (
  i in seq_len(
    nrow(
      cells
    )
  )
) {
  
  es <- cells$ES_Number[i]
  
  day <- cells$Collection_day[i]
  
  
  ctrl <- C$SampleID[
    C$ES_Number == es &
      C$Collection_day == day
  ]
  
  
  cell_worms <- W %>%
    filter(
      ES_Number == es,
      Collection_day == day
    )
  
  
  stopifnot(
    length(ctrl) == 1,
    nrow(cell_worms) == 3
  )
  
  
  high <- cell_worms$SampleID[
    which.max(
      cell_worms$apple_worms
    )
  ]
  
  
  low <- cell_worms$SampleID[
    which.min(
      cell_worms$apple_worms
    )
  ]
  
  
  v <- pc[
    high,
    ,
    drop = TRUE
  ] -
    pc[
      ctrl,
      ,
      drop = TRUE
    ]
  
  
  vector_length <- sqrt(
    sum(
      v^2
    )
  )
  
  stopifnot(
    is.finite(vector_length),
    vector_length > 0
  )
  
  
  v <- v / vector_length
  
  
  this_cell <- paste(
    es,
    day,
    sep = "_"
  )
  
  
  vecs[[this_cell]] <- v
  
  high_pots <- c(
    high_pots,
    high
  )
  
  low_pots <- c(
    low_pots,
    low
  )
  
  cell_ids <- c(
    cell_ids,
    this_cell
  )
}


V <- do.call(
  rbind,
  vecs
)

stopifnot(
  nrow(V) == 12,
  ncol(V) == K
)


# ==============================================================================
# 10. OBSERVED COSINE SIMILARITY
# ==============================================================================

cosmat <- V %*% t(V)

cos_vals <- cosmat[
  upper.tri(
    cosmat
  )
]

obs_mean_cos <- mean(
  cos_vals
)

stopifnot(
  length(cos_vals) == 66
)


# ==============================================================================
# 11. RANDOM-ORIENTATION NULL MODEL
# ==============================================================================

set.seed(1)

null_means <- replicate(
  5000,
  {
    
    rv <- matrix(
      rnorm(
        nrow(V) * K
      ),
      nrow = nrow(V),
      ncol = K
    )
    
    rv <- rv /
      sqrt(
        rowSums(
          rv^2
        )
      )
    
    rc <- rv %*% t(rv)
    
    mean(
      rc[
        upper.tri(
          rc
        )
      ]
    )
  }
)


# 95% interval of the simulated null distribution of mean cosine similarity
null_ci <- quantile(
  null_means,
  probs = c(
    0.025,
    0.975
  ),
  names = FALSE
)

ci_lo <- null_ci[1]

ci_hi <- null_ci[2]


# Two-sided simulation-based P value
P_A <- 2 * min(
  mean(
    null_means <= obs_mean_cos
  ),
  mean(
    null_means >= obs_mean_cos
  )
)

P_A <- min(
  P_A,
  1
)


cat(
  "\n============================================================\n",
  "AM_FigS5A - DIRECTIONAL ALIGNMENT\n",
  "============================================================\n"
)

cat(
  sprintf(
    "Mean cosine similarity = %+.6f\n",
    obs_mean_cos
  )
)

cat(
  sprintf(
    "95%% null interval      = [%+.6f, %+.6f]\n",
    ci_lo,
    ci_hi
  )
)

cat(
  sprintf(
    "Two-sided simulated P  = %.6f\n\n",
    P_A
  )
)


# ==============================================================================
# 12. PANEL A PLOT
# ==============================================================================

cos_df <- tibble(
  cosine = cos_vals
)

label_A <- paste0(
  "Mean cosine = ",
  sprintf("%+.3f", obs_mean_cos),
  "; P ",
  format_p(P_A),
  "\n95% null CI [",
  sprintf("%.3f", ci_lo),
  ", ",
  sprintf("%.3f", ci_hi),
  "]"
)


AM_FigS5A <- ggplot(
  cos_df,
  aes(
    x = cosine
  )
) +
  
  annotate(
    "rect",
    xmin = ci_lo,
    xmax = ci_hi,
    ymin = -Inf,
    ymax = Inf,
    fill = "lightblue",
    alpha = 0.30
  ) +
  
  geom_histogram(
    bins = 15,
    fill = "grey70",
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_vline(
    xintercept = obs_mean_cos,
    colour = "red",
    linetype = "dashed",
    linewidth = 1.0
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_A,
    hjust = 1.03,
    vjust = 1.10,
    family = "sans",
    size = 3.2,
    colour = "black",
    lineheight = 0.95
  ) +
  
  scale_x_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = c(
      -1,
      -0.5,
      0,
      0.5,
      1
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  labs(
    x = "Cosine similarity between displacement vectors",
    y = "Number of cell pairs"
  ) +
  
  theme_sfnh()


AM_FigS5A


################################################################################
################################################################################
#
# AM_FigS5B
# CONVERGENCE OF HIGHEST- VS LOWEST-NEMATODE COMMUNITIES
#
################################################################################
################################################################################


# ==============================================================================
# 13. SELECT HIGHEST AND LOWEST NEMATODE POT WITHIN EACH CELL
# ==============================================================================

W <- W %>%
  mutate(
    rank_in_cell = ave(
      apple_worms,
      cell,
      FUN = function(z) {
        rank(
          z,
          ties.method = "first"
        )
      }
    )
  )


panel_B_selection <- W %>%
  select(
    SampleID,
    ES_Number,
    Collection_day,
    apple_worms,
    cell,
    rank_in_cell
  ) %>%
  filter(
    rank_in_cell %in% c(
      1,
      3
    )
  ) %>%
  mutate(
    Panel_group = if_else(
      rank_in_cell == 3,
      "High",
      "Low"
    )
  ) %>%
  arrange(
    cell,
    Panel_group
  )


panel_B_pairs <- panel_B_selection %>%
  select(
    cell,
    Panel_group,
    SampleID
  ) %>%
  tidyr::pivot_wider(
    names_from = Panel_group,
    values_from = SampleID
  ) %>%
  arrange(
    cell
  )


stopifnot(
  nrow(panel_B_pairs) == 12,
  !any(
    is.na(
      panel_B_pairs$High
    )
  ),
  !any(
    is.na(
      panel_B_pairs$Low
    )
  )
)


hi <- panel_B_pairs$High

lo <- panel_B_pairs$Low


# ==============================================================================
# 14. PAIRWISE BRAY-CURTIS DISTANCES AMONG SELECTED COMMUNITIES
# ==============================================================================

get_pairwise_distances <- function(ids) {
  
  d <- D_all[
    ids,
    ids,
    drop = FALSE
  ]
  
  d[
    upper.tri(
      d
    )
  ]
}


high_dist <- get_pairwise_distances(
  hi
)

low_dist <- get_pairwise_distances(
  lo
)


stopifnot(
  length(hi) == 12,
  length(lo) == 12,
  length(high_dist) == 66,
  length(low_dist) == 66
)


mean_high <- mean(
  high_dist
)

mean_low <- mean(
  low_dist
)

# Negative = highest-nematode communities are more similar
# Positive = highest-nematode communities are more dissimilar
delta_obs <- mean_high - mean_low


# ==============================================================================
# 15. EXACT WITHIN-CELL PERMUTATION TEST
#
# Each of the 12 cells contributes one matched High/Low pair.
# All 2^12 = 4096 possible within-cell label assignments are evaluated.
# ==============================================================================

n_cells <- nrow(
  panel_B_pairs
)

n_exact <- 2^n_cells


swap_matrix <- sapply(
  0:(n_exact - 1),
  function(z) {
    as.logical(
      intToBits(z)[
        seq_len(
          n_cells
        )
      ]
    )
  }
)


perm_delta <- vapply(
  seq_len(
    ncol(
      swap_matrix
    )
  ),
  function(k) {
    
    swap <- swap_matrix[
      ,
      k
    ]
    
    
    perm_hi <- ifelse(
      swap,
      panel_B_pairs$Low,
      panel_B_pairs$High
    )
    
    
    perm_lo <- ifelse(
      swap,
      panel_B_pairs$High,
      panel_B_pairs$Low
    )
    
    
    mean(
      get_pairwise_distances(
        perm_hi
      )
    ) -
      mean(
        get_pairwise_distances(
          perm_lo
        )
      )
  },
  numeric(1)
)


P_B <- mean(
  abs(
    perm_delta
  ) >=
    abs(
      delta_obs
    )
)


cat(
  "\n============================================================\n",
  "AM_FigS5B - EXACT WITHIN-CELL PERMUTATION TEST\n",
  "============================================================\n"
)

cat(
  sprintf(
    "Mean Bray-Curtis, highest-nematode pots = %.6f\n",
    mean_high
  )
)

cat(
  sprintf(
    "Mean Bray-Curtis, lowest-nematode pots  = %.6f\n",
    mean_low
  )
)

cat(
  sprintf(
    "Observed difference (High - Low)       = %+.6f\n",
    delta_obs
  )
)

cat(
  sprintf(
    "Exact permutations                     = %d\n",
    n_exact
  )
)

cat(
  sprintf(
    "Two-sided exact permutation P          = %.6f\n\n",
    P_B
  )
)


# ==============================================================================
# 16. PANEL B PLOT
# ==============================================================================

conv_df <- bind_rows(
  
  tibble(
    group = "High worm",
    distance = high_dist
  ),
  
  tibble(
    group = "Low worm",
    distance = low_dist
  )
) %>%
  
  mutate(
    group = factor(
      group,
      levels = c(
        "Low worm",
        "High worm"
      )
    )
  )


label_B <- paste0(
  "High mean = ",
  sprintf(
    "%.3f",
    mean_high
  ),
  "\nLow mean = ",
  sprintf(
    "%.3f",
    mean_low
  ),
  "\nPermutation P ",
  format_p(
    P_B
  )
)

y_label_B <- max(conv_df$distance, na.rm = TRUE) -
  0.70 * diff(range(conv_df$distance, na.rm = TRUE))

AM_FigS5B <- ggplot(
  conv_df,
  aes(
    x = group,
    y = distance,
    fill = group
  )
) +
  
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    colour = "black",
    linewidth = 0.6
  ) +
  
  geom_jitter(
    width = 0.10,
    alpha = 0.45,
    size = 2
  ) +
  
  scale_fill_manual(
    values = c(
      "Low worm" = "grey80",
      "High worm" = "grey35"
    ),
    guide = "none"
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_label_B,
    label = label_B,
    hjust = 0.5,
    vjust = 1,
    family = "sans",
    size = 3.0,
    colour = "black",
    lineheight = 0.95
  ) +
  
  scale_x_discrete(
    labels = c(
      "Low worm" = "Lowest-density\npot per cell",
      "High worm" = "Highest-density\npot per cell"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Between-cell Bray-Curtis distance"
  ) +
  
  theme_sfnh()


AM_FigS5B


################################################################################
################################################################################
#
# 17. SAVE INDIVIDUAL PANELS FIRST
#
################################################################################
################################################################################


# AM_FigS5A

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS5A_Directional_Alignment.pdf"
  ),
  plot = AM_FigS5A,
  width = 150,
  height = 120,
  units = "mm"
)

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS5A_Directional_Alignment.png"
  ),
  plot = AM_FigS5A,
  width = 150,
  height = 120,
  units = "mm",
  dpi = 600,
  bg = "white"
)


# AM_FigS5B

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS5B_High_vs_Low_Convergence.pdf"
  ),
  plot = AM_FigS5B,
  width = 150,
  height = 120,
  units = "mm"
)

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS5B_High_vs_Low_Convergence.png"
  ),
  plot = AM_FigS5B,
  width = 150,
  height = 120,
  units = "mm",
  dpi = 600,
  bg = "white"
)


################################################################################
################################################################################
#
# 18. SAVE TABLES
#
################################################################################
################################################################################


write.csv(
  tibble(
    cosine_similarity = cos_vals
  ),
  file.path(
    tables_dir,
    "AM_FigS5A_Observed_Cosine_Similarities.csv"
  ),
  row.names = FALSE
)

write.csv(
  tibble(
    null_mean_cosine = null_means
  ),
  file.path(
    tables_dir,
    "AM_FigS5A_Random_Vector_Null_Distribution.csv"
  ),
  row.names = FALSE
)

write.csv(
  panel_B_selection,
  file.path(
    tables_dir,
    "AM_FigS5B_Selected_High_Low_Pots.csv"
  ),
  row.names = FALSE
)

write.csv(
  tibble(
    permuted_difference = perm_delta
  ),
  file.path(
    tables_dir,
    "AM_FigS5B_Exact_Within_Cell_Permutation.csv"
  ),
  row.names = FALSE
)


################################################################################
################################################################################
#
# 19. COMBINE AM_FigS5
#
################################################################################
################################################################################


AM_FigS5 <- (
  AM_FigS5A |
    AM_FigS5B
) +
  
  plot_layout(
    widths = c(
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
      size = 22
    )
  )


AM_FigS5


################################################################################
################################################################################
#
# 20. SAVE COMBINED FIGURE
#
################################################################################
################################################################################


ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS5.pdf"
  ),
  plot = AM_FigS5,
  width = 310,
  height = 125,
  units = "mm"
)

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS5.png"
  ),
  plot = AM_FigS5,
  width = 310,
  height = 125,
  units = "mm",
  dpi = 600,
  bg = "white"
)


################################################################################
################################################################################
#
# 21. SESSION INFORMATION
#
################################################################################
################################################################################


capture.output(
  sessionInfo(),
  file = file.path(
    tables_dir,
    "AM_FigS5_sessionInfo.txt"
  )
)


cat(
  "\n============================================================\n",
  "AM_FigS5 COMPLETE\n",
  "============================================================\n\n"
)

cat(
  "Saved individual panels and combined figure to:\n",
  figures_dir,
  "\n\n"
)

cat(
  "Combined files:\n",
  file.path(
    figures_dir,
    "AM_FigS5.pdf"
  ),
  "\n",
  file.path(
    figures_dir,
    "AM_FigS5.png"
  ),
  "\n"
)

