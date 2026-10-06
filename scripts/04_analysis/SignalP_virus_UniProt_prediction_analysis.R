# Analysis of SignalP 6.0 predictions for viral proteins
# with signal peptides annotated in UniProt.
#
# The script:
# - analyses proteins classified as OTHER by SignalP 6.0;
# - compares signal peptide lengths between predicted and non-predicted proteins;
# - evaluates prediction probabilities;
# - identifies high-confidence SignalP predictions.

# Necessary libraries
library(dplyr)
library(readr)

#### Results from SignalP 6.0, UniProt data for Viruses (signal peptide annotated)
prediction_results_viruses_signalp6_uniprot <- read.delim(
  "results/signalp/viruses_uniprot/prediction_results_viruses_uniprot_may_2026_1579.txt",
  sep = "\t",
  header = TRUE,
  skip = 1,
  comment.char = "",
  quote = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
) 

names(prediction_results_viruses_signalp6_uniprot)[1] <- "ID"

### Table of 359 non-predicted proteins by SignalP6 (UniProt, viruses annotated)

no_signalp_nodublons_viruses_entries_metadata <- read_csv(
  "data/processed/uniprot/no_signalp_nodublons_viruses_entries_metadata.csv"
)

no_signalp_prediction_merge <- no_signalp_nodublons_viruses_entries_metadata %>%
  left_join(prediction_results_viruses_signalp6_uniprot, by = c("ProteinId_input" = "ID"))

# UniProt metadata for viral proteins with annotated signal peptides
uniprotkb_viruses_taxonomy_id_10239_AND_ft_metadata <- read_csv(
  "data/processed/uniprot/uniprotkb_taxonomy_viruses_AND_ft_sign_before_update_metadata.csv"
)

max(no_signalp_prediction_merge$`SP(Sec/SPI)`)


png(
  filename = "results/figures/359_n0_predictions_signalp6.png",
  width = 8.87,
  height = 6.65,
  units = "in",
  res = 300
)

boxplot(
  no_signalp_prediction_merge[, c(
    "SP(Sec/SPI)",
    "LIPO(Sec/SPII)",
    "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)",
    "PILIN(Sec/SPIII)"
  )],
  names = c("Sec/SPI", "Sec/SPII", "Tat/SPI", "Tat/SPII", "Sec/SPIII"),
  ylab = "Probability",
  main = "SignalP6 prediction for 359 OTHER SP in Virus (UniProt data)",
  ylim = c(0, 0.55),
  las = 2,
  col = c(
    "pink",
    "skyblue",
    "lightgreen",
    "salmon",
    "plum"
  )
)

dev.off()

sum(
  no_signalp_prediction_merge$`SP(Sec/SPI)` >= 0.2 &
    no_signalp_prediction_merge$`SP(Sec/SPI)` <= 0.5,
  na.rm = TRUE
)
# Create probability intervals of 0.1
sec_spi_counts <- no_signalp_prediction_merge %>%
  mutate(
    probability_interval = cut(
      `SP(Sec/SPI)`,
      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c(
        "0–0.1",
        "0.1–0.2",
        "0.2–0.3",
        "0.3–0.4",
        "0.4–0.5",
        "≥ 0.5"
      )
    )
  ) %>%
  count(probability_interval, .drop = FALSE)

sec_spi_counts

png(
  filename = "results/figures/359_no_predictions_sec_spi_probability_fr.png",
  width = 8.87,
  height = 6.65,
  units = "in",
  res = 300
)

bar_colors <- c(
  "white",
  "white",
  "pink",
  "pink",
  "pink",
  "white"
)

bp <- barplot(
  sec_spi_counts$n,
  names.arg = sec_spi_counts$probability_interval,
  xlab = "Probabilité Sec/SPI",
  ylab = "Nombre de séquences",
  main = "Distribution des probabilités Sec/SPI\npour les 359 séquences sans peptide signal prédit",
  col = bar_colors,
  ylim = c(0, max(sec_spi_counts$n) * 1.15)
)

# Add counts above bars
text(
  x = bp,
  y = sec_spi_counts$n,
  labels = sec_spi_counts$n,
  pos = 3,
  cex = 0.9
)

dev.off()

### French version

png(
  filename = "results/figures/359_n0_predictions_signalp6_fr.png",
  width = 8.87,
  height = 6.65,
  units = "in",
  res = 300
)

boxplot(
  no_signalp_prediction_merge[, c(
    "SP(Sec/SPI)",
    "LIPO(Sec/SPII)",
    "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)",
    "PILIN(Sec/SPIII)"
  )],
  names = c("Sec/SPI", "Sec/SPII", "Tat/SPI", "Tat/SPII", "Sec/SPIII"),
  ylab = "Probabilité",
  main = "Prédiction de SignalP 6.0 pour 359 OTHER SP chez les virus (données UniProt)",
  ylim = c(0, 0.55),
  las = 2,
  col = c(
    "pink",
    "skyblue",
    "lightgreen",
    "salmon",
    "plum"
  )
)

dev.off()


### length for 355 non-predicted proteins by SignalP6 (UniProt, viruses annotated)

boxplot(no_signalp_prediction_merge$Pos_sp_end_uniprot_1, na.rm = TRUE) 

max(no_signalp_prediction_merge$Pos_sp_end_uniprot_1, na.rm = TRUE) 

prediction_results_viruses_signalp6_uniprot_SP <- prediction_results_viruses_signalp6_uniprot %>%
  filter(prediction_results_viruses_signalp6_uniprot$Prediction == "SP")

yes_signalp_prediction_merge <- prediction_results_viruses_signalp6_uniprot_SP %>%
  left_join(uniprotkb_viruses_taxonomy_id_10239_AND_ft_metadata, by = c("ID" = "ProteinId_input"))

### length for 1171 predicted proteins by SignalP6 (UniProt, viruses annotated)
boxplot(yes_signalp_prediction_merge$Pos_sp_end_uniprot_1, na.rm = TRUE) 

x_no <- na.omit(no_signalp_prediction_merge$Pos_sp_end_uniprot_1)
x_yes <- na.omit(yes_signalp_prediction_merge$Pos_sp_end_uniprot_1)

length(x_no)
length(x_yes)

test_result <- wilcox.test(
  x_no,
  x_yes,
  paired = FALSE,
  exact = FALSE
)

test_result

p <- test_result$p.value
p_label <- paste0("p = ", format.pval(p, digits = 3, eps = 0.001))


median(x_no)
median(x_yes)

wilcox_eff <- abs(
  qnorm(test_result$p.value / 2) /
    sqrt(length(x_no) + length(x_yes))
)

wilcox_eff

### Signal peptide lengths differed significantly between the groups (Wilcoxon rank-sum test), with a moderate effect size, r = 0.31.


png(
  filename = "results/figures/yes_no_predictions_signalp6_virus_uniprot.png",
  width = 7.66,
  height = 6.65,
  units = "in",
  res = 300
)

positions <- c(1, 1.5)


boxplot(
  x_no,
  x_yes,
  at = positions,
  boxwex = 0.25,
  names = c("No (n = 359)", "Yes (n = 1146)"),
  ylab = "Length of signal peptide (UniProt)",
  xlab = "Prediction by SignalP 6.0",
  main = "Comparison of SP length (UniProt dataset, viruses)",
  col = c("pink", "skyblue"),
  border = "grey30",
  ylim = c(0, 155),
  xlim = c(0.7, 1.8),
  las = 1
)

y <- 142
h <- 3

segments(positions[1], y, positions[2], y)
segments(positions[1], y, positions[1], y - h)
segments(positions[2], y, positions[2], y - h)

# p-value
text(
  x = mean(positions),
  y = y + 4,
  labels = p_label,
  cex = 0.9
)

dev.off()

#### French version

png(
  filename = "results/figures/yes_no_predictions_signalp6_virus_uniprot_fr.png",
  width = 7.66,
  height = 6.65,
  units = "in",
  res = 300
)

positions <- c(1, 1.5)


boxplot(
  x_no,
  x_yes,
  at = positions,
  boxwex = 0.25,
  names = c("Non (n = 359)", "Oui (n = 1146)"),
  ylab = "Longueur du peptide signal (UniProt)",
  xlab = "Prédiction par SignalP 6.0",
  main = "Comparaison de la longueur des peptides signaux 
  (jeu de données UniProt, virus)",
  col = c("pink", "skyblue"),
  border = "grey30",
  ylim = c(0, 155),
  xlim = c(0.7, 1.8),
  las = 1
)

y <- 142
h <- 3

segments(positions[1], y, positions[2], y)
segments(positions[1], y, positions[1], y - h)
segments(positions[2], y, positions[2], y - h)

# p-value
text(
  x = mean(positions),
  y = y + 4,
  labels = p_label,
  cex = 0.9
)

dev.off()
#####



high_prediction_score_sp_viruses <- prediction_results_viruses_signalp6_uniprot_SP %>%
  filter(prediction_results_viruses_signalp6_uniprot_SP$`SP(Sec/SPI)` > 0.999)

summary(high_prediction_score_sp_viruses$`SP(Sec/SPI)`)


signalp_viruses_uniprot_results <- read.table(
  "results/signalp/viruses_uniprot/output_signalp_viruses_uniprot.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(signalp_viruses_uniprot_results) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

signalp_viruses_uniprot_results_high_prediction_merge <- high_prediction_score_sp_viruses %>%
  left_join(signalp_viruses_uniprot_results, by = c("ID" = "ProteinId"))













