################################################################################
################################################################################
##
## AM_FigS3
## Sequential-removal sensitivity analysis of rare ASV richness
##
################################################################################
################################################################################


## =============================================================================
## PACKAGES
## =============================================================================

library(tidyverse)
library(vegan)
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
## INPUT FILES
## =============================================================================

metadata_file <- file.path(input_dir, "metadata.csv")
counts_file   <- file.path(input_dir, "bac_counts_rarefied14000.csv")
class_file    <- file.path(input_dir, "bac_abundance_class.csv")


## =============================================================================
## READ AND FORMAT DATA
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

storage.mode(cnt_ref) <- "numeric"

abcl <- read.csv(
  class_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta <- meta |>
  mutate(
    ES_Number = factor(ES_Number, levels = c("ES1", "ES2", "ES3")),
    Collection_day = factor(
      Collection_day,
      levels = c("12D", "15D", "18D", "21D")
    )
  )


## =============================================================================
## APPLE SAMPLES FROM NEMATODE-INOCULATED POTS
## =============================================================================

W <- meta |>
  filter(
    Compartment == "As",
    Treatment == "Worm",
    SampleID %in% colnames(cnt_ref)
  ) |>
  arrange(ES_Number, Collection_day, Pot_Replicate) |>
  mutate(apple_worms = worm)

stopifnot(nrow(W) == 36)

cnt_apple <- cnt_ref[, W$SampleID, drop = FALSE]
stopifnot(identical(W$SampleID, colnames(cnt_apple)))


## =============================================================================
## RARE-ASV SET
## =============================================================================

rare_asv <- abcl$ASV_ID[abcl$class == "Rare"]
rare_asv <- intersect(rare_asv, rownames(cnt_apple))


## =============================================================================
## PARTIAL SPEARMAN FUNCTION
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
    p = ct$p.value,
    n = length(y)
  )
}

cov_ed <- W[, c("ES_Number", "Collection_day")]
worm_predictor <- log1p(W$apple_worms)


## =============================================================================
## RELATIVE ABUNDANCE AND TESTABLE ASVs
## =============================================================================

rel_apple <- sweep(
  cnt_apple,
  2,
  colSums(cnt_apple),
  "/"
) * 100

prevalence <- rowSums(cnt_apple > 0)

testable_asv <- rownames(cnt_apple)[
  prevalence >= 0.25 * nrow(W)
]


## =============================================================================
## RANK ASVs BY POSITIVE ASSOCIATION WITH LOCAL NEMATODE ABUNDANCE
## =============================================================================

exp_rho <- vapply(
  testable_asv,
  function(asv) {
    as.numeric(
      partial_spearman(
        rel_apple[asv, ],
        worm_predictor,
        cov_ed
      )["rho"]
    )
  },
  numeric(1)
)

names(exp_rho) <- testable_asv
exp_rho <- sort(exp_rho, decreasing = TRUE)

write_csv(
  tibble(
    ASV_ID = names(exp_rho),
    partial_rho = as.numeric(exp_rho),
    rank = seq_along(exp_rho)
  ),
  file.path(tables_dir, "AM_FigS3_Ranked_Positive_ASVs.csv")
)


## =============================================================================
## SEQUENTIAL REMOVAL AND RE-RAREFACTION
## =============================================================================

removal_n <- c(0, 2, 5, 10, 15, 20, 30)
rarefaction_seed <- 7

run_expander_curve <- function(seed = rarefaction_seed) {
  map_dfr(
    removal_n,
    function(nk) {
      drop_asv <- if (nk == 0) {
        character(0)
      } else {
        names(exp_rho)[seq_len(nk)]
      }

      remaining_asv <- setdiff(rownames(cnt_apple), drop_asv)
      sub_counts <- cnt_apple[remaining_asv, , drop = FALSE]
      rare_remaining <- intersect(rare_asv, remaining_asv)
      depth <- min(colSums(sub_counts))

      if (nk == 0) {
        # cnt_apple is already rarefied to 14,000 reads per sample.
        richness <- colSums(
          sub_counts[rare_remaining, , drop = FALSE] > 0
        )
      } else {
        # Fixed seed makes the stochastic re-rarefaction reproducible.
        set.seed(seed)
        rarefied <- t(
          vegan::rrarefy(
            t(sub_counts),
            sample = depth
          )
        )

        richness <- colSums(
          rarefied[
            intersect(rare_asv, rownames(rarefied)),
            ,
            drop = FALSE
          ] > 0
        )
      }

      stat <- partial_spearman(
        richness,
        worm_predictor,
        cov_ed
      )

      high12 <- order(
        worm_predictor,
        decreasing = TRUE
      )[seq_len(min(12, length(worm_predictor)))]

      share_hi <- if (length(drop_asv) == 0) {
        0
      } else {
        mean(
          colSums(
            rel_apple[drop_asv, , drop = FALSE]
          )[high12]
        )
      }

      tibble(
        n_removed = nk,
        depth = depth,
        rho = as.numeric(stat["rho"]),
        p = as.numeric(stat["p"]),
        share_hi = share_hi
      )
    }
  )
}

curve_apple <- run_expander_curve()

print(curve_apple)

write_csv(
  curve_apple,
  file.path(tables_dir, "AM_FigS3_Sequential_Removal.csv")
)


## =============================================================================
## PLOT
## =============================================================================

curve_plot <- curve_apple |>
  mutate(
    sig = if_else(p < 0.05, "P < 0.05", "P >= 0.05")
  )

rho_range <- c(
  min(curve_plot$rho, 0) - 0.05,
  max(curve_plot$rho, 0) + 0.05
)

share_range <- c(
  0,
  max(curve_plot$share_hi, 1)
)

scale_share <- function(x) {
  rho_range[1] +
    (x - share_range[1]) / diff(share_range) * diff(rho_range)
}

unscale_share <- function(y) {
  share_range[1] +
    (y - rho_range[1]) / diff(rho_range) * diff(share_range)
}

p_s3 <- ggplot(
  curve_plot,
  aes(x = n_removed, y = rho)
) +
  geom_line(
    aes(y = scale_share(share_hi)),
    colour = "grey65",
    linetype = 2,
    linewidth = 1.0
  ) +
  geom_point(
    aes(y = scale_share(share_hi)),
    colour = "grey50",
    fill = "grey85",
    shape = 22,
    size = 4.2,
    stroke = 0.8
  ) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_line(colour = "black", linewidth = 1.1) +
  geom_point(
    aes(fill = sig),
    shape = 21,
    colour = "black",
    size = 4.8,
    stroke = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "P < 0.05" = "black",
      "P >= 0.05" = "white"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    name = "Partial Spearman rho",
    sec.axis = sec_axis(
      ~ unscale_share(.),
      name = "Removed sequence reads in high-nematode pots (%)"
    )
  ) +
  scale_x_continuous(breaks = removal_n) +
  coord_cartesian(ylim = rho_range) +
  labs(x = "Positively associated ASVs removed (n)") +
  theme_classic(base_size = 20, base_family = "sans") +
  theme(
    text = element_text(colour = "black"),
    legend.position = "top",
    axis.title.y.right = element_text(colour = "grey40"),
    axis.text.y.right = element_text(colour = "grey40")
  )

p_s3

ggsave(
  file.path(figures_dir, "AM_FigS3_Sequential_Removal.pdf"),
  p_s3,
  width = 160,
  height = 115,
  units = "mm"
)

ggsave(
  file.path(figures_dir, "AM_FigS3_Sequential_Removal.png"),
  p_s3,
  width = 160,
  height = 115,
  units = "mm",
  dpi = 600,
  bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "AM_FigS3_sessionInfo.txt")
)

