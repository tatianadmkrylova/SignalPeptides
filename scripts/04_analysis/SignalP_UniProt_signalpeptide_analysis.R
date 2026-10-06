# Analysis of signal peptide annotations and SignalP 6.0 predictions
#
# This script compares SignalP 6.0 predictions with signal peptide
# annotations from UniProt and signalpeptide.de for viruses, bacteria,
# and mammals.
#
# The analysis includes:
# - comparison of predicted and annotated signal peptide lengths;
# - comparison between UniProt and signalpeptide.de annotations;
# - calculation of SignalP 6.0 performance metrics for the viral dataset;
# - generation of figures used in the project analysis.

# Necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(ggpointdensity)
library(viridis)

####################### Viruses ################################################

####################### UniProt ################################################

######## Necessary tables (data) ###############################################

# Table from SignalP 6.0 results (predicted signal peptides)
results_viruses_signalp6_uniprot_1579 <- read.table(
  "results/signalp/viruses_uniprot/output_signalp_viruses_uniprot.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(results_viruses_signalp6_uniprot_1579) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

results_viruses_signalp6_uniprot_1579 <- results_viruses_signalp6_uniprot_1579[0:6]

# Table from SignalP 6.0 results for all protein sequences
prediction_results_viruses_signalp6_uniprot_1579 <- read.delim(
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

colnames(prediction_results_viruses_signalp6_uniprot_1579)[1] <- "ID"

# Table with UniProt metadata
uniprotkb_viruses_metadata_1579 <- read_csv(
  "data/processed/uniprot/uniprotkb_taxonomy_viruses_AND_ft_sign_before_update_metadata.csv"
)


### Merged table (UniProt metadata with SignalP results)

prediction_results_viruses_signalp6_uniprot_1579_metadata <-
  results_viruses_signalp6_uniprot_1579 %>%
  left_join(
    uniprotkb_viruses_metadata_1579,
    by = c("ProteinId" = "ProteinId_input")
  )


### Comparison of signal peptide lengths: SignalP 6.0 vs UniProt

plot_signal_len_compare_virus_uni_signalp <- ggplot(
  prediction_results_viruses_signalp6_uniprot_1579_metadata,
  aes(x = end, y = Pos_sp_end_uniprot_1, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP 6.0") +
  ylab("Length in UniProt") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Virus (N = 1202)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_signal_len_compare_virus_uni_signalp


ggsave(
  filename = "results/figures/plot_signal_len_compare_virus_uni_signalp.png",
  plot = plot_signal_len_compare_virus_uni_signalp,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

### French version

### Comparison of signal peptide lengths: SignalP 6.0 vs UniProt

plot_signal_len_compare_virus_uni_signalp_fr <- ggplot(
  prediction_results_viruses_signalp6_uniprot_1579_metadata,
  aes(x = end, y = Pos_sp_end_uniprot_1, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Longueur dans SignalP 6.0") +
  ylab("Longueur dans UniProt") +
  scale_color_viridis() +
  labs(title = "Longueur des peptides signaux chez les virus 
       (N = 1202)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_signal_len_compare_virus_uni_signalp_fr


ggsave(
  filename = "results/figures/plot_signal_len_compare_virus_uni_signalp_fr.png",
  plot = plot_signal_len_compare_virus_uni_signalp_fr,
  width = 5.83,
  height = 4.55,
  dpi = 300
)
######################### No predicted SP ######################################

prediction_results_viruses_signalp6_uniprot_1579_all_prob <-
  prediction_results_viruses_signalp6_uniprot_1579 %>%
  left_join(
    results_viruses_signalp6_uniprot_1579,
    by = c("ID" = "ProteinId")
  )

prediction_results_viruses_signalp6_uniprot_1579_all_prob_uniprot_merge <-
  prediction_results_viruses_signalp6_uniprot_1579_all_prob %>%
  left_join(
    uniprotkb_viruses_metadata_1579,
    by = c("ID" = "ProteinId_input")
  )

# Protein IDs without SignalP prediction
out_prediction_viruses_sp_uniprot <-
  prediction_results_viruses_signalp6_uniprot_1579_all_prob_uniprot_merge %>%
  filter(Prediction == "OTHER")

# Protein IDs with SignalP prediction
with_prediction_viruses_sp_uniprot <-
  prediction_results_viruses_signalp6_uniprot_1579_all_prob_uniprot_merge %>%
  filter(Prediction != "OTHER")

TP <- length(with_prediction_viruses_sp_uniprot$ID)
FN <- length(out_prediction_viruses_sp_uniprot$ID)


######################### Negative dataset #####################################

# SignalP 6.0 predictions for all proteins in the negative dataset
prediction_negatif_set <- read.delim(
  "results/signalp/negative_dataset/prediction_results.txt",
  sep = "\t",
  header = TRUE,
  skip = 1,
  comment.char = "",
  quote = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# SignalP 6.0 positive predictions in the negative dataset
results_prediction_negatif_set <- read.table(
  "results/signalp/negative_dataset/output.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

FP <- length(results_prediction_negatif_set$V1)
TN <- length(prediction_negatif_set$`# ID`) - FP

TP <- as.numeric(TP)
TN <- as.numeric(TN)
FP <- as.numeric(FP)
FN <- as.numeric(FN)


######################### Performance metrics ##################################

precision <- TP / (TP + FP)

recall <- TP / (TP + FN)

specificity <- TN / (TN + FP)

accuracy <- (TP + TN) / (TP + TN + FP + FN)

F1_score <- 2 * precision * recall /
  (precision + recall)

MCC <- (TP * TN - FP * FN) /
  sqrt(
    (TP + FP) *
      (TP + FN) *
      (TN + FP) *
      (TN + FN)
  )

balanced_accuracy <- (recall + specificity) / 2

results <- data.frame(
  Metric = c(
    "Precision",
    "Recall",
    "Specificity",
    "Accuracy",
    "F1-score",
    "MCC",
    "Balanced accuracy"
  ),
  Value = c(
    precision,
    recall,
    specificity,
    accuracy,
    F1_score,
    MCC,
    balanced_accuracy
  )
)

results


######################### Confusion matrix #####################################

conf_df <- data.frame(
  Actual = c("SP", "SP", "No SP", "No SP"),
  Predicted = c("SP", "No SP", "SP", "No SP"),
  Count = c(TP, FN, FP, TN)
)

conf_df

conf_df$Actual <- factor(
  conf_df$Actual,
  levels = c("SP", "No SP")
)

conf_df$Predicted <- factor(
  conf_df$Predicted,
  levels = c("SP", "No SP")
)

conf_matrix_vir <- ggplot(
  conf_df,
  aes(x = Predicted, y = Actual, fill = Count)
) +
  geom_tile(color = "black", linewidth = 0.8) +
  geom_text(aes(label = Count), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Predicted class",
    y = "Actual class"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )

conf_matrix_vir

ggsave(
  filename = "results/figures/conf_matrix_vir.png",
  plot = conf_matrix_vir,
  width = 4.79,
  height = 4.55,
  dpi = 300
)

### French version

conf_matrix_vir_fr <- ggplot(
  conf_df,
  aes(x = Predicted, y = Actual, fill = Count)
) +
  geom_tile(color = "black", linewidth = 0.8) +
  geom_text(aes(label = Count), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Classe prédite",
    y = "Classe réelle"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )

conf_matrix_vir_fr

ggsave(
  filename = "results/figures/conf_matrix_vir_fr.png",
  plot = conf_matrix_vir_fr,
  width = 4.79,
  height = 4.55,
  dpi = 300
)

####################### Signal Peptide Web database ############################

######## Necessary tables (data) ###############################################

# Table from SignalP 6.0 results (predicted signal peptides)
results_viruses_signalp6_spw <- read.table(
  "results/signalp/viruses_signalpeptide_de/output_signalp_viruses_spw.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(results_viruses_signalp6_spw) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

results_viruses_signalp6_spw <- results_viruses_signalp6_spw[0:6]


# Table with UniProt metadata
uniprotkb_viruses_metadata_spw <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_viruses_all_metadata.csv"
)

# Table from signalpeptide.de
viruses_metadata_site_spw <- read_csv(
  "data/raw/signalpeptide_de/signalpeptide_viruses.csv"
)


### Merged table (signalpeptide.de data with SignalP results)

results_viruses_signalp6_spw_site_merge <- results_viruses_signalp6_spw %>%
  left_join(
    viruses_metadata_site_spw,
    by = c("ProteinId" = "Accession Number")
  )


### Comparison of signal peptide lengths: SignalP 6.0 vs signalpeptide.de

plot_results_viruses_signalp6_spw_site_merge <- ggplot(
  results_viruses_signalp6_spw_site_merge,
  aes(x = end, y = Length, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP 6.0") +
  ylab("Length in SPW") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Virus (N = 5553)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_results_viruses_signalp6_spw_site_merge


ggsave(
  filename = "results/figures/plot_results_viruses_signalp6_spw_site_merge.png",
  plot = plot_results_viruses_signalp6_spw_site_merge,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

#### French version

### Comparison of signal peptide lengths: SignalP 6.0 vs signalpeptide.de

plot_results_viruses_signalp6_spw_site_merge_fr <- ggplot(
  results_viruses_signalp6_spw_site_merge,
  aes(x = end, y = Length, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Longueur dans SignalP 6.0") +
  ylab("Longueur dans SPW") +
  scale_color_viridis() +
  labs(title = "Longueur des peptides signaux chez les virus 
       (N = 5553)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_results_viruses_signalp6_spw_site_merge_fr


ggsave(
  filename = "results/figures/plot_results_viruses_signalp6_spw_site_merge_fr.png",
  plot = plot_results_viruses_signalp6_spw_site_merge_fr,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

######################## Merge UniProt vs SPW ##################################

uniprotkb_viruses_metadata_1579_viruses_metadata_site_spw_merge <-
  uniprotkb_viruses_metadata_1579 %>%
  left_join(
    viruses_metadata_site_spw,
    by = c("ProteinId_input" = "Accession Number")
  )


intersection_viruses <- uniprotkb_viruses_metadata_1579 %>%
  inner_join(
    viruses_metadata_site_spw,
    by = c("ProteinId_input" = "Accession Number")
  )

nrow(intersection_viruses)


plot_data_viruses <-
  uniprotkb_viruses_metadata_1579_viruses_metadata_site_spw_merge %>%
  filter(
    !is.na(Pos_sp_end_uniprot_1),
    !is.na(Length)
  )


plot_uniprotkb_viruses_metadata_1579_viruses_metadata_site_spw_merge <- ggplot(
  plot_data_viruses,
  aes(
    x = Pos_sp_end_uniprot_1,
    y = Length
  )
) +
  theme_classic() +
  xlab("Length in UniProt") +
  ylab("Length in SPW") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Virus (N=1078)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_uniprotkb_viruses_metadata_1579_viruses_metadata_site_spw_merge


ggsave(
  filename = "results/figures/plot_uniprotkb_viruses_metadata_1579_viruses_metadata_site_spw_merge.png",
  plot = plot_uniprotkb_viruses_metadata_1579_viruses_metadata_site_spw_merge,
  width = 5.83,
  height = 4.55,
  dpi = 300
)


# Number and percentage of identical signal peptide lengths

plot_data_viruses %>%
  summarise(
    N = n(),
    identical = sum(Pos_sp_end_uniprot_1 == Length),
    percentage_identical =
      mean(Pos_sp_end_uniprot_1 == Length) * 100
  )



####################### Bacteria ################################################

####################### UniProt ################################################

######## Necessary tables (data) ###############################################

# Table from SignalP 6.0 results (predicted signal peptides)
results_bacteria_signalp6_uniprot_824 <- read.table(
  "results/signalp/bacteria_uniprot/output_signalp_bacteria_uniprot.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(results_bacteria_signalp6_uniprot_824) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

results_bacteria_signalp6_uniprot_824 <-
  results_bacteria_signalp6_uniprot_824[0:6]


# Table from SignalP 6.0 results for all protein sequences
prediction_results_bacteria_signalp6_uniprot_824 <- read.delim(
  "results/signalp/bacteria_uniprot/prediction_results_bacteria_uniprot.txt",
  sep = "\t",
  header = TRUE,
  skip = 1,
  comment.char = "",
  quote = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(prediction_results_bacteria_signalp6_uniprot_824)[1] <- "ID"


# Table with UniProt metadata
uniprotkb_bacteria_metadata_824 <- read_csv(
  "data/processed/uniprot/uniprot_bacteria_signal_confirmed_exp_before_update_metadata.csv"
)


### Merged table (UniProt metadata with SignalP results)

prediction_results_bacteria_signalp6_uniprot_824_metadata <-
  results_bacteria_signalp6_uniprot_824 %>%
  left_join(
    uniprotkb_bacteria_metadata_824,
    by = c("ProteinId" = "ProteinId_input")
  )


### Comparison of signal peptide lengths: SignalP 6.0 vs UniProt

plot_signal_len_compare_bacteria_uni_signalp <- ggplot(
  prediction_results_bacteria_signalp6_uniprot_824_metadata,
  aes(x = end, y = Pos_sp_end_uniprot_1, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP 6.0") +
  ylab("Length in UniProt") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Bacteria (N = 824)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_signal_len_compare_bacteria_uni_signalp


ggsave(
  filename = "results/figures/plot_signal_len_compare_bacteria_uni_signalp.png",
  plot = plot_signal_len_compare_bacteria_uni_signalp,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

##### French version

### Comparison of signal peptide lengths: SignalP 6.0 vs UniProt

plot_signal_len_compare_bacteria_uni_signalp_fr <- ggplot(
  prediction_results_bacteria_signalp6_uniprot_824_metadata,
  aes(x = end, y = Pos_sp_end_uniprot_1, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Longueur dans SignalP 6.0") +
  ylab("Longueur dans UniProt") +
  scale_color_viridis() +
  labs(title = "Longueur des peptides signaux chez les bactéries
       (N = 806)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_signal_len_compare_bacteria_uni_signalp_fr

ggsave(
  filename = "results/figures/plot_signal_len_compare_bacteria_uni_signalp_fr.png",
  plot = plot_signal_len_compare_bacteria_uni_signalp_fr,
  width = 5.83,
  height = 4.55,
  dpi = 300
)


####################### Signal Peptide Web database ############################

######## Necessary tables (data) ###############################################

# Table from SignalP 6.0 results (predicted signal peptides)
results_bacteria_signalp6_spw <- read.table(
  "results/signalp/bacteria_signalpeptide_de/output_bacteria_confirmed.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(results_bacteria_signalp6_spw) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

results_bacteria_signalp6_spw <-
  results_bacteria_signalp6_spw[0:6]


# Table with UniProt metadata
uniprotkb_bacteria_metadata_spw <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_bacteria_uniprot_meta.csv"
)


# Table from signalpeptide.de
bacteria_metadata_site_spw <- read_csv(
  "data/raw/signalpeptide_de/signalpeptide_bacteria.csv"
)

# Keep only confirmed signal peptides
bacteria_metadata_site_spw_cfrm <- bacteria_metadata_site_spw %>%
  filter(`SP Status` == "confirmed")


### Merged table (signalpeptide.de data with SignalP results)

results_bacteria_signalp6_spw_site_merge <- results_bacteria_signalp6_spw %>%
  left_join(
    bacteria_metadata_site_spw_cfrm,
    by = c("ProteinId" = "Accession Number")
  )

results_bacteria_signalp6_spw_site_merge$Length <-
  as.numeric(results_bacteria_signalp6_spw_site_merge$Length)


### Comparison of signal peptide lengths: SignalP 6.0 vs signalpeptide.de

plot_results_bacteria_signalp6_spw_site_merge <- ggplot(
  results_bacteria_signalp6_spw_site_merge,
  aes(x = end, y = Length, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP 6.0") +
  ylab("Length in SPW") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Bacteria (N = 1161)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_results_bacteria_signalp6_spw_site_merge


ggsave(
  filename = "results/figures/plot_results_bacteria_signalp6_spw_site_merge.png",
  plot = plot_results_bacteria_signalp6_spw_site_merge,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

### French version

plot_results_bacteria_signalp6_spw_site_merge_fr <- ggplot(
  results_bacteria_signalp6_spw_site_merge,
  aes(x = end, y = Length, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Longueur dans SignalP 6.0") +
  ylab("Longueur dans SPW") +
  scale_color_viridis() +
  labs(title = "Longueur des peptides signaux chez les bactéries 
       (N = 1118)")+
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_results_bacteria_signalp6_spw_site_merge_fr

ggsave(
  filename = "results/figures/plot_results_bacteria_signalp6_spw_site_merge_fr.png",
  plot = plot_results_bacteria_signalp6_spw_site_merge_fr,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

######################## Merge UniProt vs SPW ##################################

uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge <-
  uniprotkb_bacteria_metadata_824 %>%
  left_join(
    bacteria_metadata_site_spw_cfrm,
    by = c("ProteinId_input" = "Accession Number")
  )

uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge$Length <-
  as.numeric(
    uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge$Length
  )


intersection_bacteria <- uniprotkb_bacteria_metadata_824 %>%
  inner_join(
    bacteria_metadata_site_spw_cfrm,
    by = c("ProteinId_input" = "Accession Number")
  )

nrow(intersection_bacteria)


plot_data_bacteria <-
  uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge %>%
  filter(
    !is.na(Pos_sp_end_uniprot_1),
    !is.na(Length)
  )


plot_uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge <- ggplot(
  plot_data_bacteria,
  aes(
    x = Pos_sp_end_uniprot_1,
    y = Length
  )
) +
  theme_classic() +
  xlab("Length in UniProt") +
  ylab("Length in SPW") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Bacteria (N=643)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge


ggsave(
  filename = "results/figures/plot_uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge.png",
  plot = plot_uniprotkb_bacteria_metadata_824_bacteria_metadata_site_spw_cfr_merge,
  width = 5.83,
  height = 4.55,
  dpi = 300
)


# Number and percentage of identical signal peptide lengths

plot_data_bacteria %>%
  summarise(
    N = n(),
    identical = sum(Pos_sp_end_uniprot_1 == Length),
    percentage_identical =
      mean(Pos_sp_end_uniprot_1 == Length) * 100
  )

####################### Mammalia ################################################

####################### UniProt ################################################

######## Necessary tables (data) ###############################################

# Table from SignalP 6.0 results (predicted signal peptides)
results_mammalia_signalp6_uniprot_1463 <- read.table(
  "results/signalp/mammalia_uniprot/output_signalp_mammalia_uniprot.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(results_mammalia_signalp6_uniprot_1463) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

results_mammalia_signalp6_uniprot_1463 <-
  results_mammalia_signalp6_uniprot_1463[0:6]


# Table from SignalP 6.0 results for all protein sequences
prediction_results_mammalia_signalp6_uniprot_1463 <- read.delim(
  "results/signalp/mammalia_uniprot/prediction_results_mammalia_uniprot.txt",
  sep = "\t",
  header = TRUE,
  skip = 1,
  comment.char = "",
  quote = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(prediction_results_mammalia_signalp6_uniprot_1463)[1] <- "ID"


# Table with UniProt metadata
uniprotkb_mammalia_metadata_1463 <- read_csv(
  "data/processed/uniprot/uniprot_mammalia_signal_confirmed_exp_before_update_metadata.csv"
)


### Merged table (UniProt metadata with SignalP results)

prediction_results_mammalia_signalp6_uniprot_1463_metadata <-
  results_mammalia_signalp6_uniprot_1463 %>%
  left_join(
    uniprotkb_mammalia_metadata_1463,
    by = c("ProteinId" = "ProteinId_input")
  )


### Comparison of signal peptide lengths: SignalP 6.0 vs UniProt

plot_signal_len_compare_mammalia_uni_signalp <- ggplot(
  prediction_results_mammalia_signalp6_uniprot_1463_metadata,
  aes(x = end, y = Pos_sp_end_uniprot_1, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP 6.0") +
  ylab("Length in UniProt") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Mammalia (N = 1463)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_signal_len_compare_mammalia_uni_signalp


ggsave(
  filename = "results/figures/plot_signal_len_compare_mammalia_uni_signalp.png",
  plot = plot_signal_len_compare_mammalia_uni_signalp,
  width = 5.83,
  height = 4.55,
  dpi = 300
)


#### French version

### Comparison of signal peptide lengths: SignalP 6.0 vs UniProt

plot_signal_len_compare_mammalia_uni_signalp_fr <- ggplot(
  prediction_results_mammalia_signalp6_uniprot_1463_metadata,
  aes(x = end, y = Pos_sp_end_uniprot_1, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Longueur dans SignalP 6.0") +
  ylab("Longueur dans UniProt") +
  scale_color_viridis() +
  labs(title = "Longueur des peptides signaux chez les mammifères 
       (N = 1436)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_signal_len_compare_mammalia_uni_signalp_fr


ggsave(
  filename = "results/figures/plot_signal_len_compare_mammalia_uni_signalp_fr.png",
  plot = plot_signal_len_compare_mammalia_uni_signalp_fr,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

####################### Signal Peptide Web database ############################

######## Necessary tables (data) ###############################################

# Table from SignalP 6.0 results (predicted signal peptides)
results_mammalia_signalp6_spw <- read.table(
  "results/signalp/mammalia_signalpeptide_de/output_mammalia_confirmed.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(results_mammalia_signalp6_spw) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

results_mammalia_signalp6_spw <-
  results_mammalia_signalp6_spw[0:6]

results_mammalia_signalp6_spw$ProteinId <-
  sub("\\|.*$", "", results_mammalia_signalp6_spw$ProteinId)


# Table with UniProt metadata
uniprotkb_mammalia_metadata_spw <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_mammalia_uniprot_meta.csv"
)


# Table from signalpeptide.de
mammalia_metadata_site_spw <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_mammalia_confirmed.csv"
)


### Merged table (signalpeptide.de data with SignalP results)

results_mammalia_signalp6_spw_site_merge <- results_mammalia_signalp6_spw %>%
  left_join(
    mammalia_metadata_site_spw,
    by = c("ProteinId" = "ProteinId")
  )


### Comparison of signal peptide lengths: SignalP 6.0 vs signalpeptide.de

plot_results_mammalia_signalp6_spw_site_merge <- ggplot(
  results_mammalia_signalp6_spw_site_merge,
  aes(x = end, y = Length, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP 6.0") +
  ylab("Length in SPW") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Mammalia (N = 2061)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_results_mammalia_signalp6_spw_site_merge


ggsave(
  filename = "results/figures/plot_results_mammalia_signalp6_spw_site_merge.png",
  plot = plot_results_mammalia_signalp6_spw_site_merge,
  width = 5.83,
  height = 4.55,
  dpi = 300
)

#### French version
### Comparison of signal peptide lengths: SignalP 6.0 vs signalpeptide.de

plot_results_mammalia_signalp6_spw_site_merge_fr <- ggplot(
  results_mammalia_signalp6_spw_site_merge,
  aes(x = end, y = Length, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Longueur dans SignalP 6.0") +
  ylab("Longueur dans SPW") +
  scale_color_viridis() +
  labs(title = "Longueur des peptides signaux chez les mammifères 
       (N = 2061)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_results_mammalia_signalp6_spw_site_merge_fr


ggsave(
  filename = "results/figures/plot_results_mammalia_signalp6_spw_site_merge_fr.png",
  plot = plot_results_mammalia_signalp6_spw_site_merge_fr,
  width = 5.83,
  height = 4.55,
  dpi = 300
)


######################## Merge UniProt vs SPW ##################################

uniprotkb_mammalia_metadata_1463_mammalia_metadata_site_spw_merge <-
  uniprotkb_mammalia_metadata_1463 %>%
  left_join(
    mammalia_metadata_site_spw,
    by = c("ProteinId_input" = "ProteinId")
  )


intersection_mammalia <- uniprotkb_mammalia_metadata_1463 %>%
  inner_join(
    mammalia_metadata_site_spw,
    by = c("ProteinId_input" = "ProteinId")
  )

nrow(intersection_mammalia)


plot_data_mammalia <-
  uniprotkb_mammalia_metadata_1463_mammalia_metadata_site_spw_merge %>%
  filter(
    !is.na(Pos_sp_end_uniprot_1),
    !is.na(Length)
  )


plot_uniprotkb_mammalia_metadata_1463_mammalia_metadata_site_spw_merge <- ggplot(
  plot_data_mammalia,
  aes(
    x = Pos_sp_end_uniprot_1,
    y = Length
  )
) +
  theme_classic() +
  xlab("Length in UniProt") +
  ylab("Length in SPW") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Mammalia (N=1261)") +
  geom_pointdensity(show.legend = TRUE) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE,
    color = "grey30",
    linetype = "dashed",
    linewidth = 0.5,
    fill = "grey70",
    alpha = 0.25
  ) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  coord_fixed(ratio = 1)

plot_uniprotkb_mammalia_metadata_1463_mammalia_metadata_site_spw_merge


ggsave(
  filename = "results/figures/plot_uniprotkb_mammalia_metadata_1463_mammalia_metadata_site_spw_merge.png",
  plot = plot_uniprotkb_mammalia_metadata_1463_mammalia_metadata_site_spw_merge,
  width = 5.83,
  height = 4.55,
  dpi = 300
)


# Number and percentage of identical signal peptide lengths

plot_data_mammalia %>%
  summarise(
    N = n(),
    identical = sum(Pos_sp_end_uniprot_1 == Length),
    percentage_identical =
      mean(Pos_sp_end_uniprot_1 == Length) * 100
  )

