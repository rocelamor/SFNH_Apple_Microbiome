################################################################################
################################################################################
##
## AM_FigS1
## Total nematode abundance and whole-pot developmental composition
##
## A  Total nematode abundance per pot (apple + soil)
## B  Whole-pot developmental-stage composition
##
################################################################################
################################################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(vegan)
  library(FSA)
})

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


## Read metadata
meta <- read.csv(
  file.path(input_dir, "metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta$ES_Number <- factor(
  meta$ES_Number,
  levels = c("ES1", "ES2", "ES3")
)

meta$Collection_day <- factor(
  meta$Collection_day,
  levels = c("12D", "15D", "18D", "21D")
)

meta <- meta %>%
  mutate(
    pot_id = paste(
      ES_Number,
      Collection_day,
      Pot_Replicate,
      sep = "_"
    )
  )

## =============================================================================
## POT-LEVEL TOTAL WORM POPULATION
## =============================================================================

pot_summary <- meta %>%
  group_by(
    pot_id,
    ES_Number,
    Collection_day,
    Pot_Replicate,
    Treatment
  ) %>%
  summarise(
    apple_worms = sum(
      worm[Compartment == "As"],
      na.rm = TRUE
    ),
    soil_worms = sum(
      worm[Compartment %in% c("Ss", "Sd")],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    total_worms = apple_worms + soil_worms
  )

W_pot <- pot_summary %>%
  filter(Treatment == "Worm")


## =============================================================================
## NEGATIVE-BINOMIAL GLM FOR TOTAL NEMATODE ABUNDANCE
## =============================================================================

# Full model: ecological setup, collection day, and their interaction.
m_total_full <- MASS::glm.nb(
  total_worms ~ ES_Number * Collection_day,
  data = W_pot
)

# Additive model used for the likelihood-ratio test of the interaction.
m_total_additive <- MASS::glm.nb(
  total_worms ~ ES_Number + Collection_day,
  data = W_pot
)

interaction_test <- anova(
  m_total_additive,
  m_total_full,
  test = "LRT"
)

cat("\n=== ES x collection-day interaction: negative-binomial GLM ===\n")
print(interaction_test)

write_csv(
  as.data.frame(interaction_test) %>% rownames_to_column("model"),
  file.path(tables_dir, "AM_FigS1_ESxDay_LRT.csv")
)

# Simulation-based residual diagnostics.
set.seed(123)
total_dharma <- DHARMa::simulateResiduals(
  fittedModel = m_total_full,
  n = 1000
)

pdf(file.path(figures_dir, "AM_FigS1_NB_GLM_Residual_Diagnostics.pdf"))
plot(total_dharma)
dev.off()

diagnostic_tests <- tibble(
  test = c("dispersion", "zero_inflation", "uniformity"),
  p_value = c(
    DHARMa::testDispersion(total_dharma)$p.value,
    DHARMa::testZeroInflation(total_dharma)$p.value,
    DHARMa::testUniformity(total_dharma)$p.value
  )
)

write_csv(
  diagnostic_tests,
  file.path(tables_dir, "AM_FigS1_NB_GLM_Diagnostics.csv")
)

# Pairwise ecological-setup comparisons within collection day.
emm_ES_within_day <- emmeans::emmeans(
  m_total_full,
  ~ ES_Number | Collection_day,
  type = "response"
)

ES_within_day_pairs <- pairs(
  emm_ES_within_day,
  adjust = "tukey"
)

write_csv(
  as.data.frame(emm_ES_within_day),
  file.path(tables_dir, "AM_FigS1_Estimated_Marginal_Means.csv")
)

write_csv(
  as.data.frame(ES_within_day_pairs),
  file.path(tables_dir, "AM_FigS1_ES_Within_Day_Tukey.csv")
)

## =============================================================================
## RESULTS TABLES
## =============================================================================

max_total <- W_pot %>%
  arrange(desc(total_worms)) %>%
  slice(1)

total_by_es_day <- W_pot %>%
  group_by(
    ES_Number,
    Collection_day
  ) %>%
  summarise(
    n = n(),
    mean_total = mean(total_worms, na.rm = TRUE),
    median_total = median(total_worms, na.rm = TRUE),
    min_total = min(total_worms, na.rm = TRUE),
    max_total = max(total_worms, na.rm = TRUE),
    .groups = "drop"
  )

day21_summary <- total_by_es_day %>%
  filter(Collection_day == "21D")

zero_total <- W_pot %>%
  filter(total_worms == 0) %>%
  arrange(
    ES_Number,
    Collection_day,
    Pot_Replicate
  )

## =============================================================================
## WHOLE-POT DEVELOPMENTAL COMPOSITION
## =============================================================================

total_stage <- meta %>%
  filter(Treatment == "Worm") %>%
  group_by(
    pot_id,
    ES_Number,
    Collection_day,
    Pot_Replicate
  ) %>%
  summarise(
    L1_L2 = sum(L1_L2, na.rm = TRUE),
    L3_L4 = sum(L3_L4, na.rm = TRUE),
    Adult = sum(Adult, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_staged = L1_L2 + L3_L4 + Adult,
    p_L1_L2 = if_else(
      total_staged > 0,
      100 * L1_L2 / total_staged,
      NA_real_
    ),
    p_L3_L4 = if_else(
      total_staged > 0,
      100 * L3_L4 / total_staged,
      NA_real_
    ),
    p_Adult = if_else(
      total_staged > 0,
      100 * Adult / total_staged,
      NA_real_
    )
  )

stage_summary <- total_stage %>%
  select(
    ES_Number,
    Collection_day,
    p_L1_L2,
    p_L3_L4,
    p_Adult
  ) %>%
  pivot_longer(
    starts_with("p_"),
    names_to = "stage",
    values_to = "pct"
  ) %>%
  mutate(
    stage = recode(
      stage,
      p_L1_L2 = "L1-L2",
      p_L3_L4 = "L3-L4",
      p_Adult = "Adult"
    ),
    stage = factor(
      stage,
      levels = c("L1-L2", "L3-L4", "Adult")
    )
  ) %>%
  group_by(
    ES_Number,
    Collection_day,
    stage
  ) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    .groups = "drop"
  )

## Console summaries
cat("\n=== Highest individual total-nematode pot ===\n")
print(max_total)

cat("\n=== ES x day total-nematode summary ===\n")
print(total_by_es_day, n = Inf)

cat("\n=== Day 21 summary ===\n")
print(day21_summary)

cat("\n=== Recorded total-nematode zeros ===\n")
print(zero_total, n = Inf)

cat("\n=== Whole-pot developmental-stage means ===\n")
print(stage_summary, n = Inf)

## Save tables
write_csv(
  W_pot,
  file.path(
    tables_dir,
    "AM_FigS1_Total_Nematode_Pot_Data.csv"
  )
)

write_csv(
  total_by_es_day,
  file.path(
    tables_dir,
    "AM_FigS1_Total_Nematode_ES_Day_Summary.csv"
  )
)

write_csv(
  stage_summary,
  file.path(
    tables_dir,
    "AM_FigS1_Stage_Composition_Summary.csv"
  )
)

## =============================================================================
## PLOT STYLE
## =============================================================================

ES_COL <- c(
  ES1 = "#264653",
  ES2 = "#2a9d8f",
  ES3 = "#e76f51"
)

STAGE_COL <- c(
  "L1-L2" = "#d9d9d9",
  "L3-L4" = "#8c8c8c",
  "Adult" = "#262626"
)

theme_sfnh <- function(base_size = 12) {
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      text = element_text(
        family = "sans",
        colour = "black"
      ),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.4
      ),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        colour = "black"
      ),
      legend.title = element_blank()
    )
}

## =============================================================================
## PANEL A — TOTAL WORMS PER POT
## =============================================================================

p_total_worm <- ggplot(
  W_pot,
  aes(
    x = as.numeric(
      str_remove(
        as.character(Collection_day),
        "D"
      )
    ),
    y = total_worms + 1,
    colour = ES_Number
  )
) +
  geom_point(
    alpha = 0.60,
    size = 2.5,
    position = position_jitter(
      width = 0.15,
      height = 0
    )
  ) +
  stat_summary(
    fun = median,
    geom = "line",
    linewidth = 0.9,
    aes(group = ES_Number)
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 3.1,
    shape = 21,
    fill = "white",
    stroke = 0.9
  ) +
  scale_y_log10(
    labels = scales::label_number(
      big.mark = ",",
      accuracy = 1
    )
  ) +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21),
    labels = c("12D", "15D", "18D", "21D")
  ) +
  scale_colour_manual(values = ES_COL) +
  labs(
    x = "Incubation day",
    y = expression(
        "Nematodes per pot (+1, log scale)"
    )
  ) +
  theme_sfnh() +
  theme(
    legend.position = "top"
  )

## =============================================================================
## PANEL B — WHOLE-POT DEVELOPMENTAL COMPOSITION
## =============================================================================

p_stage_total <- ggplot(
  stage_summary,
  aes(
    x = Collection_day,
    y = mean_pct,
    fill = stage
  )
) +
  geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  facet_wrap(
    ~ ES_Number,
    nrow = 1
  ) +
  scale_fill_manual(values = STAGE_COL) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Incubation day",
    y = "Population (%)"
  ) +
  theme_sfnh() +
  theme(
    legend.position = "top"
  )

## =============================================================================
## COMBINED SUPPLEMENTARY FIGURE
## =============================================================================

fig_total_worm <- (
  p_total_worm /
    p_stage_total
) +
  plot_layout(
    heights = c(1, 1.1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "top",
    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 16
    )
  )

fig_total_worm

## =============================================================================
## SAVE
## =============================================================================

ggsave(
  file.path(
    figures_dir,
    "AM_FigS1_Total_Nematode_Population.pdf"
  ),
  fig_total_worm,
  width = 170,
  height = 180,
  units = "mm"
)

ggsave(
  file.path(
    figures_dir,
    "AM_FigS1_Total_Nematode_Population.png"
  ),
  fig_total_worm,
  width = 170,
  height = 180,
  units = "mm",
  dpi = 600,
  bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "AM_FigS1_sessionInfo.txt"
  )
)

cat(
  "\nSaved:\n",
  "AM_FigS1_Total_Nematode_Population.pdf\n",
  "AM_FigS1_Total_Nematode_Population.png\n"
)

## =============================================================================

