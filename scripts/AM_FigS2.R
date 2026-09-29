################################################################################
################################################################################
##
## AM_FigS2
## Bacterial community structure across the SFNH system and within the apple compartment
##
## A  Apple versus soil bacterial community composition
## B  Apple bacterial community composition versus local nematode abundance
##
################################################################################
################################################################################


## =============================================================================
## 1. PACKAGES
## =============================================================================

library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)
library(scales)
library(ggtext)


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
## 3. INPUT FILES
## =============================================================================

metadata_file <- file.path(
  input_dir,
  "metadata.csv"
)

bacteria_file <- file.path(
  input_dir,
  "bac_counts_rarefied14000.csv"
)


## =============================================================================
## 4. READ DATA
## =============================================================================

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

bacteria <- read_csv(
  bacteria_file,
  show_col_types = FALSE
)


## =============================================================================
## 5. FORMAT METADATA
## =============================================================================

metadata <- metadata %>%
  mutate(
    
    ES_Number = factor(
      ES_Number,
      levels = c(
        "ES1",
        "ES2",
        "ES3"
      )
    ),
    
    day_num = as.numeric(
      stringr::str_extract(
        as.character(Collection_day),
        "\\d+"
      )
    ),
    
    Collection_day = factor(
      day_num,
      levels = c(
        12,
        15,
        18,
        21
      ),
      labels = c(
        "12D",
        "15D",
        "18D",
        "21D"
      )
    ),
    
    Compartment = factor(
      Compartment,
      levels = c(
        "As",
        "Ss",
        "Sd"
      )
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Worm"
      )
    )
  )


## =============================================================================
## 6. COLOURS
## =============================================================================

ES_COL <- c(
  "ES1" = "#0072B2",
  "ES2" = "#D55E00",
  "ES3" = "#009E73"
)

REGION_COL <- c(
  "Apple" = "#CC79A7",
  "Soil" = "#0072B2"
)


## =============================================================================
## 7. PUBLICATION THEME
## =============================================================================

theme_sfnh <- function(
    base_size = 14
) {
  
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
        size = base_size + 1,
        colour = "black"
      ),
      
      axis.text = element_text(
        size = base_size,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.6,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.6,
        colour = "black"
      ),
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = base_size - 1
      ),
      
      legend.position = "top",
      
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      
      strip.text = element_text(
        size = base_size + 1,
        face = "bold"
      ),
      
      plot.margin = margin(
        8,
        10,
        8,
        10
      )
    )
}


## =============================================================================
## 8. PREPARE BACTERIAL COUNT MATRIX
## =============================================================================

## First column contains ASV IDs

asv_id_col <- colnames(
  bacteria
)[1]


## ASVs = rows, samples = columns

bac_matrix <- bacteria %>%
  column_to_rownames(
    var = asv_id_col
  ) %>%
  as.matrix()


storage.mode(
  bac_matrix
) <- "numeric"


## Transpose:
## samples = rows
## ASVs = columns

bac_samples <- t(
  bac_matrix
)


cat(
  "\nRaw bacterial matrix:\n"
)

print(
  dim(
    bac_samples
  )
)


## =============================================================================
## 9. MATCH COUNT TABLE TO METADATA
## =============================================================================

common_samples <- intersect(
  rownames(
    bac_samples
  ),
  metadata$SampleID
)


cat(
  "\nSamples shared between bacterial table and metadata:",
  length(
    common_samples
  ),
  "\n"
)


bac_samples <- bac_samples[
  common_samples,
  ,
  drop = FALSE
]


metadata_bac <- metadata %>%
  
  filter(
    SampleID %in%
      common_samples
  ) %>%
  
  arrange(
    match(
      SampleID,
      rownames(
        bac_samples
      )
    )
  )


stopifnot(
  identical(
    metadata_bac$SampleID,
    rownames(
      bac_samples
    )
  )
)


cat(
  "Final bacterial samples:",
  nrow(
    bac_samples
  ),
  "\n"
)

cat(
  "Final ASVs:",
  ncol(
    bac_samples
  ),
  "\n"
)


################################################################################
################################################################################
##
## AM_FigS2A
##
## OVERALL BACTERIAL COMMUNITY STRUCTURE
##
## IMPORTANT:
## Soil surface (Ss) and deep soil (Sd) are NOT summed into synthetic samples.
## They remain independent bacterial-community samples and are assigned
## the same plotting/statistical region label: Soil.
##
################################################################################
################################################################################


## =============================================================================
## 10. BRAY-CURTIS DISTANCE
## =============================================================================

bray_all <- vegdist(
  bac_samples,
  method = "bray"
)


## =============================================================================
## 11. PCoA
## =============================================================================

pcoa_all <- cmdscale(
  bray_all,
  k = 2,
  eig = TRUE
)


## =============================================================================
## 12. VARIANCE REPRESENTED BY PCoA AXES
## =============================================================================

eig_all <- pcoa_all$eig


var_all <- eig_all[
  eig_all > 0
] /
  
  sum(
    eig_all[
      eig_all > 0
    ]
  ) *
  
  100


PC1_all <- round(
  var_all[
    1
  ],
  1
)


PC2_all <- round(
  var_all[
    2
  ],
  1
)


## =============================================================================
## 13. PCoA DATAFRAME
## =============================================================================

pcoa_all_df <- as.data.frame(
  pcoa_all$points
)


colnames(
  pcoa_all_df
) <- c(
  "PCoA1",
  "PCoA2"
)


pcoa_all_df$SampleID <- rownames(
  pcoa_all_df
)


pcoa_all_df <- pcoa_all_df %>%
  
  left_join(
    metadata_bac,
    by = "SampleID"
  ) %>%
  
  mutate(
    
    Region = case_when(
      as.character(Compartment) == "As" ~ "Apple",
      as.character(Compartment) %in% c("Ss", "Sd") ~ "Soil",
      TRUE ~ NA_character_
    ),
    
    Region = factor(
      Region,
      levels = c(
        "Apple",
        "Soil"
      )
    )
  )


## =============================================================================
## 14. MARGINAL PERMANOVA — PANEL A
##
## Question:
##
## Does region (Apple vs Soil) explain bacterial community composition
## after accounting for:
##
##   ES_Number
##   Collection_day
##   Treatment
##
## by = "margin" means each term is evaluated while the other terms
## remain in the model.
##
## =============================================================================

perm_data_A <- metadata_bac %>%
  mutate(
    
    Region = case_when(
      as.character(Compartment) == "As" ~ "Apple",
      as.character(Compartment) %in% c("Ss", "Sd") ~ "Soil",
      TRUE ~ NA_character_
    ),
    
    Region = factor(
      Region,
      levels = c(
        "Apple",
        "Soil"
      )
    )
  )


stopifnot(
  !any(is.na(perm_data_A$Region)),
  identical(
    perm_data_A$SampleID,
    labels(
      bray_all
    )
  )
)

cat(
  "\nAM_FigS2A region counts:\n"
)

print(
  table(
    perm_data_A$Region
  )
)


set.seed(
  1
)


perm_A <- adonis2(
  bray_all ~
    ES_Number +
    Collection_day +
    Treatment +
    Region,
  data = perm_data_A,
  permutations = 9999,
  by = "margin"
)


cat(
  "\n============================================================\n",
  "AM_FigS2A — MARGINAL PERMANOVA\n",
  "============================================================\n"
)


print(
  perm_A
)


## Extract Apple-vs-Soil region result

R2_A <- perm_A[
  "Region",
  "R2"
]


P_A <- perm_A[
  "Region",
  "Pr(>F)"
]


cat(
  sprintf(
    "\nRegion (Apple vs Soil): R2 = %.3f, P = %.6f\n",
    R2_A,
    P_A
  )
)


## Save full PERMANOVA table

perm_A_table <- as.data.frame(
  perm_A
) %>%
  rownames_to_column(
    "term"
  )


write_csv(
  perm_A_table,
  file.path(
    tables_dir,
    "AM_FigS2A_PERMANOVA.csv"
  )
)


## =============================================================================
## 15. PANEL A STATISTIC LABEL
## =============================================================================

perm_label_A <- paste0(
  "PERMANOVA<br>",
  "Region: R<sup>2</sup> = ", sprintf("%.3f", R2_A), "<br>",
  "<i>P</i> ",
  ifelse(
    P_A < 0.001,
    "&lt; 0.001",
    paste0("= ", sprintf("%.3f", P_A))
  )
)

## =============================================================================
## 16. AM_FigS2A
## =============================================================================

p2a <- ggplot(
  pcoa_all_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Region
  )
) +
  
  geom_point(
    size = 3.4,
    alpha = 0.75
  ) +
  
  ggtext::geom_richtext(
    data = data.frame(
      x = Inf,
      y = Inf,
      label = perm_label_A
    ),
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    family = "sans",
    size = 4.0,
    colour = "black",
    fill = NA,
    label.colour = NA,
    linewidth = 0,
    label.padding = unit(0, "lines"),
    lineheight = 1.15
  ) +
  
  scale_colour_manual(
    values = REGION_COL,
    name = "Region"
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.08
      )
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.20
      )
    )
  ) +
  
  labs(
    x = paste0(
      "PCoA1 (",
      PC1_all,
      "%)"
    ),
    
    y = paste0(
      "PCoA2 (",
      PC2_all,
      "%)"
    )
  ) +
  
  theme_sfnh(
    base_size = 14
  ) +
  
  theme(
    legend.position = "bottom"
  )


p2a


################################################################################
################################################################################
##
## AM_FigS2B
##
## APPLE BACTERIAL COMMUNITY STRUCTURE
##
## nematode-inoculated pots only
##
################################################################################
################################################################################


## =============================================================================
## 17. APPLE / NEMATODE-INOCULATED SAMPLES
##
## We use only nematode-inoculated pots here because local apple nematode abundance is
## the biological predictor.
##
## n should be 36.
## =============================================================================

apple_worm_meta <- metadata_bac %>%
  
  filter(
    Compartment == "As",
    Treatment == "Worm"
  ) %>%
  
  arrange(
    match(
      SampleID,
      rownames(
        bac_samples
      )
    )
  ) %>%
  
  mutate(
    apple_worms = worm,
    
    worm_log1p = log1p(
      worm
    )
  )


apple_worm_ids <- apple_worm_meta$SampleID


bac_apple_worm <- bac_samples[
  apple_worm_ids,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    apple_worm_meta$SampleID,
    rownames(
      bac_apple_worm
    )
  )
)


cat(
  "\nApple/nematode-inoculated samples for AM_FigS2B:",
  nrow(
    apple_worm_meta
  ),
  "\n"
)


## =============================================================================
## 18. BRAY-CURTIS — APPLE/WORM
## =============================================================================

bray_apple_worm <- vegdist(
  bac_apple_worm,
  method = "bray"
)


## =============================================================================
## 19. PCoA — APPLE/WORM
## =============================================================================

pcoa_apple_worm <- cmdscale(
  bray_apple_worm,
  k = 2,
  eig = TRUE
)


## =============================================================================
## 20. VARIANCE REPRESENTED BY PCoA AXES
## =============================================================================

eig_apple <- pcoa_apple_worm$eig


var_apple <- eig_apple[
  eig_apple > 0
] /
  
  sum(
    eig_apple[
      eig_apple > 0
    ]
  ) *
  
  100


PC1_apple <- round(
  var_apple[
    1
  ],
  1
)


PC2_apple <- round(
  var_apple[
    2
  ],
  1
)


## =============================================================================
## 21. APPLE PCoA DATAFRAME
## =============================================================================

pcoa_apple_df <- as.data.frame(
  pcoa_apple_worm$points
)


colnames(
  pcoa_apple_df
) <- c(
  "PCoA1",
  "PCoA2"
)


pcoa_apple_df$SampleID <- rownames(
  pcoa_apple_df
)


pcoa_apple_df <- pcoa_apple_df %>%
  
  left_join(
    apple_worm_meta,
    by = "SampleID"
  )


## =============================================================================
## 22. MARGINAL PERMANOVA — PANEL B
##
## Question:
##
## Is bacterial community composition associated with actual local
## apple worm abundance after accounting for:
##
##   ES_Number
##   Collection_day
##
## This is the beta-diversity analogue of our partial-analysis logic.
##
## =============================================================================

stopifnot(
  identical(
    apple_worm_meta$SampleID,
    labels(
      bray_apple_worm
    )
  )
)


set.seed(
  1
)


perm_B <- adonis2(
  bray_apple_worm ~
    ES_Number +
    Collection_day +
    worm_log1p,
  data = apple_worm_meta,
  permutations = 9999,
  by = "margin"
)


cat(
  "\n============================================================\n",
  "AM_FigS2B — MARGINAL PERMANOVA\n",
  "============================================================\n"
)


print(
  perm_B
)


## Extract worm-abundance term only

R2_B <- perm_B[
  "worm_log1p",
  "R2"
]


P_B <- perm_B[
  "worm_log1p",
  "Pr(>F)"
]


cat(
  sprintf(
    "\nLocal apple nematode abundance: R2 = %.3f, P = %.6f\n",
    R2_B,
    P_B
  )
)


## Save complete PERMANOVA table

perm_B_table <- as.data.frame(
  perm_B
) %>%
  rownames_to_column(
    "term"
  )


write_csv(
  perm_B_table,
  file.path(
    tables_dir,
    "AM_FigS2B_PERMANOVA.csv"
  )
)


## =============================================================================
## 23. PANEL B STATISTIC LABEL
## =============================================================================
perm_label_B <- paste0(
  "PERMANOVA<br>",
  "Worms in apple: R<sup>2</sup> = ", sprintf("%.3f", R2_B), "<br>",
  "<i>P</i> ",
  ifelse(
    P_B < 0.001,
    "&lt; 0.001",
    paste0("= ", sprintf("%.3f", P_B))
  )
)


## =============================================================================
## 24. AM_FigS2B
## =============================================================================

p2b <- ggplot(
  pcoa_apple_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = ES_Number,
    shape = Collection_day
  )
) +
  
  geom_point(
    size = 3.4,
    alpha = 0.82,
    stroke = 0.7
  ) +
  ggtext::geom_richtext(
    data = data.frame(
      x = Inf,
      y = Inf,
      label = perm_label_B
    ),
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    family = "sans",
    size = 4.0,
    colour = "black",
    fill = NA,
    label.colour = NA,
    linewidth = 0,
    label.padding = unit(0, "lines"),
    lineheight = 1.15
  ) +
  
  scale_colour_manual(
    values = ES_COL,
    name = "Ecological setup"
  ) +
  
  scale_shape_manual(
    values = c(
      "12D" = 16,
      "15D" = 17,
      "18D" = 15,
      "21D" = 18
    ),
    name = "Collection day"
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.08
      )
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.20
      )
    )
  ) +
  
  labs(
    x = paste0(
      "PCoA1 (",
      PC1_apple,
      "%)"
    ),
    
    y = paste0(
      "PCoA2 (",
      PC2_apple,
      "%)"
    )
  ) +
  
  theme_sfnh(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.box = "vertical"
  )


p2b


################################################################################
################################################################################
##
## 25. COMBINE FIGURE S2
##
################################################################################
################################################################################

fig2 <- (
  p2a |
    p2b
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
      size = 18
    )
  )


fig2


################################################################################
################################################################################
##
## 26. SAVE INDIVIDUAL PANELS
##
################################################################################
################################################################################

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS2A_Overall_Bacterial_PCoA.pdf"
  ),
  plot = p2a,
  width = 140,
  height = 115,
  units = "mm"
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS2A_Overall_Bacterial_PCoA.png"
  ),
  plot = p2a,
  width = 140,
  height = 115,
  units = "mm",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS2B_Apple_Bacterial_PCoA.pdf"
  ),
  plot = p2b,
  width = 140,
  height = 115,
  units = "mm"
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS2B_Apple_Bacterial_PCoA.png"
  ),
  plot = p2b,
  width = 140,
  height = 115,
  units = "mm",
  dpi = 600,
  bg = "white"
)


################################################################################
################################################################################
##
## 27. SAVE FINAL COMBINED FIGURE
##
################################################################################
################################################################################

ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS2_Bacterial_Community_Structure.pdf"
  ),
  plot = fig2,
  width = 330,
  height = 130,
  units = "mm"
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_FigS2_Bacterial_Community_Structure.png"
  ),
  plot = fig2,
  width = 300,
  height = 125,
  units = "mm",
  dpi = 600,
  bg = "white"
)


################################################################################
################################################################################
##
## 28. SAVE SESSION INFORMATION
##
################################################################################
################################################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "AM_FigS2_sessionInfo.txt"
  )
)


cat(
  "\n============================================================\n",
  "FIGURE S2 COMPLETE\n",
  "============================================================\n\n"
)


cat(
  "Panel A:\n",
  sprintf(
    "Region (Apple vs Soil) | ES + day + treatment: R2 = %.3f, P = %.6f\n\n",
    R2_A,
    P_A
  )
)


cat(
  "Panel B:\n",
  sprintf(
    "Actual apple worms | ES + day: R2 = %.3f, P = %.6f\n\n",
    R2_B,
    P_B
  )
)


cat(
  "Saved to:\n",
  figures_dir,
  "\n"
)

