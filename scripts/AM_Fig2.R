################################################################################
################################################################################
##
## AM_Fig2
## Local C. elegans abundance and developmental composition in decomposing apple
##
################################################################################
################################################################################


## =============================================================================
## 1. PACKAGES
## =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)


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
## 3. IMPORT METADATA
## =============================================================================

metadata_file <- file.path(
  input_dir,
  "metadata.csv"
)

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)


## Check imported data

cat("Metadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata column names:\n")
print(colnames(metadata))


## =============================================================================
## 4. CHECK REQUIRED COLUMNS
## =============================================================================

required_columns <- c(
  "SampleID",
  "ES_Number",
  "Collection_day",
  "Pot_Replicate",
  "Compartment",
  "Treatment",
  "worm",
  "L1_L2",
  "L3_L4",
  "Adult"
)

missing_columns <- setdiff(
  required_columns,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


## =============================================================================
## 5. MANUSCRIPT COLOR PALETTE
## =============================================================================

## Ecological succession colors
## Keep these IDENTICAL throughout the paper

ES_COL <- c(
  "ES1" = "#0072B2",
  "ES2" = "#D55E00",
  "ES3" = "#009E73"
)


## Developmental-stage colors

STAGE_COL <- c(
  "L1-L2" = "#56B4E9",
  "L3-L4" = "#E69F00",
  "Adult" = "#CC79A7"
)


## =============================================================================
## 6. MANUSCRIPT THEME
## =============================================================================

theme_sfnh <- function(base_size = 9) {
  
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    
    theme(
      
      text = element_text(
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
        linewidth = 0.5,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.5,
        colour = "black"
      ),
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = base_size
      ),
      
      legend.position = "top",
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = base_size + 1,
        face = "bold"
      ),
      
      plot.margin = margin(
        5, 5, 5, 5
      )
    )
}


## =============================================================================
## 7. PREPARE APPLE + NEMATODE-INOCULATED SAMPLES
## =============================================================================

apple_worm <- metadata %>%
  
  filter(
    Compartment == "As",
    Treatment == "Worm"
  ) %>%
  
  mutate(
    
    ES_Number = factor(
      ES_Number,
      levels = c(
        "ES1",
        "ES2",
        "ES3"
      )
    ),
    
    ## Extract the numeric day even if values are
    ## 12, "12", or "12D"
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
    )
  ) %>%
  
  arrange(
    ES_Number,
    day_num,
    Pot_Replicate
  )

## =============================================================================
## 8. CHECK EXPERIMENTAL DESIGN
## =============================================================================

cat("\n============================================\n")
cat("APPLE/NEMATODE-INOCULATED DATA CHECK\n")
cat("============================================\n\n")

cat("Number of Apple/nematode-inoculated samples:\n")
print(nrow(apple_worm))

cat("\nSamples per ES x collection day:\n")

print(
  table(
    apple_worm$ES_Number,
    apple_worm$Collection_day
  )
)


## We expect:
## 36 Apple/nematode-inoculated samples
## 3 replicates per ES x day


if (nrow(apple_worm) != 36) {
  
  warning(
    paste(
      "Expected 36 Apple/nematode-inoculated samples but found",
      nrow(apple_worm)
    )
  )
}


## =============================================================================
## 9. SAVE APPLE/NEMATODE-INOCULATED DATA USED FOR AM_Fig2
## =============================================================================

write_csv(
  apple_worm,
  file.path(
    tables_dir,
    "AM_Fig2_Apple_Nematode_Metadata.csv"
  )
)


################################################################################
################################################################################
##
## AM_Fig2A
## LOCAL APPLE NEMATODE ABUNDANCE
##
################################################################################
################################################################################


## =============================================================================
## 10. SUMMARIZE LOCAL APPLE NEMATODE ABUNDANCE
## =============================================================================

worm_summary <- apple_worm %>%
  
  group_by(
    ES_Number,
    Collection_day
  ) %>%
  
  summarise(
    median_worm = median(
      worm,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## Check summary

print(worm_summary)


## =============================================================================
## 11. AM_Fig2A — LOCAL APPLE NEMATODE ABUNDANCE
## =============================================================================

p2a <- ggplot() +
  
  ## Median trajectory
  geom_line(
    data = worm_summary,
    aes(
      x = Collection_day,
      y = median_worm + 1,
      colour = ES_Number,
      group = ES_Number
    ),
    linewidth = 0.8
  ) +
  
  ## Hollow median point
  geom_point(
    data = worm_summary,
    aes(
      x = Collection_day,
      y = median_worm + 1,
      colour = ES_Number
    ),
    shape = 21,
    fill = "white",
    size = 3.2,
    stroke = 1
  ) +
  
  ## Individual biological replicates
  geom_point(
    data = apple_worm,
    aes(
      x = Collection_day,
      y = worm + 1,
      colour = ES_Number
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2,
    alpha = 0.70
  ) +
  
  scale_colour_manual(
    values = ES_COL
  ) +
  
  scale_y_log10(
    breaks = c(
      1,
      10,
      100,
      1000,
      10000,
      20000
    ),
    labels = c(
      "0",
      "10",
      "100",
      "1,000",
      "10,000",
      "20,000"
    ),
    expand = expansion(
      mult = c(
        0.03,
        0.08
      )
    )
  ) +
  
  labs(
    x = "Collection day",
    y = "Nematodes recovered from apple\n(log10[count + 1])"
  ) +
  
  theme_sfnh()


p2a


## =============================================================================
## 11. SAVE AM_Fig2A
## =============================================================================

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig2A_Local_Apple_Nematode_Abundance.pdf"
  ),
  plot = p2a,
  width = 85,
  height = 75,
  units = "mm",
  device = pdf
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig2A_Local_Apple_Nematode_Abundance.png"
  ),
  plot = p2a,
  width = 85,
  height = 75,
  units = "mm",
  dpi = 600,
  bg = "white"
)


################################################################################
################################################################################
##
## AM_Fig2B
## DEVELOPMENTAL-STAGE COMPOSITION IN APPLE
##
################################################################################
################################################################################


## =============================================================================
## 12. CALCULATE DEVELOPMENTAL-STAGE COMPOSITION
## =============================================================================

apple_stage <- apple_worm %>%
  
  mutate(
    
    stage_total =
      L1_L2 +
      L3_L4 +
      Adult,
    
    p_L1L2 = if_else(
      stage_total > 0,
      L1_L2 / stage_total * 100,
      NA_real_
    ),
    
    p_L3L4 = if_else(
      stage_total > 0,
      L3_L4 / stage_total * 100,
      NA_real_
    ),
    
    p_Adult = if_else(
      stage_total > 0,
      Adult / stage_total * 100,
      NA_real_
    )
  ) %>%
  
  select(
    SampleID,
    ES_Number,
    Collection_day,
    Pot_Replicate,
    p_L1L2,
    p_L3L4,
    p_Adult
  ) %>%
  
  pivot_longer(
    
    cols = c(
      p_L1L2,
      p_L3L4,
      p_Adult
    ),
    
    names_to = "Stage",
    values_to = "Percent"
  ) %>%
  
  mutate(
    
    Stage = recode(
      Stage,
      "p_L1L2" = "L1-L2",
      "p_L3L4" = "L3-L4",
      "p_Adult" = "Adult"
    ),
    
    Stage = factor(
      Stage,
      levels = c(
        "L1-L2",
        "L3-L4",
        "Adult"
      )
    )
  )


## =============================================================================
## 13. SUMMARIZE STAGE COMPOSITION
## =============================================================================

stage_summary <- apple_worm %>%
  
  group_by(
    ES_Number,
    Collection_day
  ) %>%
  
  summarise(
    L1_L2 = sum(L1_L2, na.rm = TRUE),
    L3_L4 = sum(L3_L4, na.rm = TRUE),
    Adult = sum(Adult, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  mutate(
    Total = L1_L2 + L3_L4 + Adult
  ) %>%
  
  pivot_longer(
    cols = c(
      L1_L2,
      L3_L4,
      Adult
    ),
    names_to = "Stage",
    values_to = "Count"
  ) %>%
  
  mutate(
    Percent = Count / Total * 100,
    
    Stage = recode(
      Stage,
      "L1_L2" = "L1-L2",
      "L3_L4" = "L3-L4",
      "Adult" = "Adult"
    ),
    
    Stage = factor(
      Stage,
      levels = c(
        "Adult",
        "L3-L4",
        "L1-L2"
      )
    )
  )


## Save summary table

write_csv(
  stage_summary,
  file.path(
    tables_dir,
    "AM_Fig2B_Apple_Developmental_Stage_Summary.csv"
  )
)


## =============================================================================
## 14. PLOT DEVELOPMENTAL-STAGE COMPOSITION
## =============================================================================
## Flag ES x day combinations with no nematodes

no_worm_cells <- stage_summary %>%
  filter(Total == 0) %>%
  distinct(
    ES_Number,
    Collection_day
  )


## Keep only valid developmental-stage percentages

stage_plot <- stage_summary %>%
  filter(
    Total > 0,
    is.finite(Percent)
  )


## Plot

p2b <- ggplot(
  stage_plot,
  aes(
    x = Collection_day,
    y = Percent,
    fill = Stage
  )
) +
  
  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  
  ## Mark ES x day cells where no nematodes were recovered
  geom_text(
    data = no_worm_cells,
    aes(
      x = Collection_day,
      y = 50,
      label = "No nematodes\nrecovered"
    ),
    inherit.aes = FALSE,
    size = 2.7,
    lineheight = 0.9
  ) +
  
  facet_wrap(
    ~ ES_Number,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = STAGE_COL,
    breaks = c(
      "L1-L2",
      "L3-L4",
      "Adult"
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      100,
      25
    ),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  
  labs(
    x = "Collection day",
    y = "Developmental stage composition (%)"
  ) +
  
  theme_sfnh() +
  
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.6
    ),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      margin = margin(
        t = 4,
        r = 4,
        b = 4,
        l = 4
      )
    )
  )

p2b

p2b <- p2b +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.6
    ),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      margin = margin(
        t = 4,
        r = 4,
        b = 4,
        l = 4
      )
    )
  )

p2b
## =============================================================================
## 15. SAVE AM_Fig2B
## =============================================================================

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig2B_Apple_Developmental_Stages.pdf"
  ),
  plot = p2b,
  width = 170,
  height = 70,
  units = "mm"
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig2B_Apple_Developmental_Stages.png"
  ),
  plot = p2b,
  width = 170,
  height = 80,
  units = "mm",
  dpi = 600,
  bg = "white"
)


################################################################################
################################################################################
##
## COMBINED AM_Fig2A + 1C
##
##
################################################################################
################################################################################


## =============================================================================
## 16. COMBINE PANELS
## =============================================================================

fig2 <- p2a / p2b +
  
  plot_layout(
    heights = c(
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
      size = 12
    )
  )


## Display

fig2


## =============================================================================
## 17. SAVE COMBINED FIGURE
## =============================================================================

ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig2_Local_Nematode_Abundance_Development.pdf"
  ),
  plot = fig2,
  width = 170,
  height = 150,
  units = "mm",
  device = pdf
)


ggsave(
  filename = file.path(
    figures_dir,
    "AM_Fig2_Local_Nematode_Abundance_Development.png"
  ),
  plot = fig2,
  width = 170,
  height = 150,
  units = "mm",
  dpi = 600,
  bg = "white"
)


## =============================================================================
## 18. SESSION INFORMATION
## =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "AM_Fig2_sessionInfo.txt"
  )
)


cat("\n============================================\n")
cat("AM_Fig2 COMPLETE\n")
cat("============================================\n")

cat("\nFiles saved to:\n")
cat(figures_dir, "\n")

cat("\nTables saved to:\n")
cat(tables_dir, "\n")

