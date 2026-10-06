# Analysis of signal peptide annotations across UniProt, signalpeptide.de,
# and the SignalP training set for viruses, bacteria, and mammals.
#
# The script:
# - compares UniProt and signalpeptide.de datasets;
# - analyses signal peptide length distributions;
# - generates Venn diagrams;
# - analyses ECO evidence codes;
# - compares the datasets with the SignalP training set.

# Necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggvenn)
library(ggpubr)
library(ggpointdensity)
library(viridis)


####################### Viruses ################################################

#### Files of input 

# Read the file from UniProt
uniprot_viruses_ft_sign_metadata <- read_csv(
  "data/processed/uniprot/uniprotkb_taxonomy_viruses_AND_ft_sign_before_update_metadata.csv"
)

# Read the file from signalpeptide.de
signalpeptide_viruses <- read_csv(
  "data/raw/signalpeptide_de/signalpeptide_viruses.csv"
)
signalpeptide_viruses$Length <- as.numeric(signalpeptide_viruses$Length)

# signalpeptide.de proteins with metadata retrieved from UniProt
signalpeptide_viruses_uniprot_meta <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_viruses_all_metadata.csv"
)

### To filter the lines with NA in Primary Accession number in UniProt (this protein doesn't exist)
signalpeptide_viruses_uniprot_meta_filtered <- signalpeptide_viruses_uniprot_meta[!is.na(signalpeptide_viruses_uniprot_meta$PrimaryAccession), ]


#### Intersection between signalpeptide.de et UniProt

# Merge signalpeptide_viruses and signalpeptide_viruses_uniprot_meta 

signalpeptide_meta_merged_viruses <- merge(signalpeptide_viruses_uniprot_meta_filtered, signalpeptide_viruses, by.x = "ProteinId_input",
                                                         by.y = "Accession Number")

signalpeptide_meta_merged_viruses_uniprot <- merge(signalpeptide_meta_merged_viruses, uniprot_viruses_ft_sign_metadata, by = "PrimaryAccession") 


# Venn diagram

signalpeptide_uniprot_meta_merged_list_viruses <- list(
  "signalpeptide.de with 
  UniProt metadata (n = 1678)" = signalpeptide_viruses_uniprot_meta_filtered$PrimaryAccession,
  "UniProt with annotated
  signal peptide (n = 1579)" = uniprot_viruses_ft_sign_metadata$PrimaryAccession)

venn_signalpeptide_uniprot_meta_merged_viruses <- ggvenn(signalpeptide_uniprot_meta_merged_list_viruses, 
                                              fill_color = c("#0073C2FF", "pink"),
                                              set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalpeptide_uniprot_meta_merged_viruses <- venn_signalpeptide_uniprot_meta_merged_viruses + 
  labs(title = "Viruses") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalpeptide_uniprot_meta_merged_viruses)

ggsave(
  filename = "results/figures/venn_signalpeptide_uniprot_meta_merged_viruses.png",
  plot = venn_signalpeptide_uniprot_meta_merged_viruses,
  width = 6.93,
  height = 5.51,
  dpi = 300
)

# Statistics

# signalpeptide.de length vs uniprot with annotated signal peptide
# Intersection (N=1085)

test_corr_viruses_int <- cor.test(signalpeptide_meta_merged_viruses_uniprot$Length, signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y, method = "spearman", exact = FALSE)
test_corr_viruses_int

plot_signal_len_compare_viruses <- ggplot(
  signalpeptide_meta_merged_viruses_uniprot,
  aes(x = Length, y = Pos_sp_end_uniprot_1.y, na.rm = TRUE)
) +
  geom_point() +
  geom_pointdensity(show.legend = TRUE) +
  theme_classic() +
  xlab("Length in signalpeptide.de") +
  ylab("Length in UniProt") +
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
  labs(title = "Signal peptide lengths in Viruses (N = 1085)") +
  scale_color_viridis() +
  scale_x_continuous(limits = c(0, 150), breaks = seq (0, 150, by = 10)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq (0, 150, by = 10)) +
  stat_cor(
    method = "spearman",
    label.x = 10,
    label.y = 130
  ) +
  coord_fixed(ratio = 1)


plot_signal_len_compare_viruses

ggsave(
  filename = "results/figures/plot_signal_len_compare_viruses.png",
  plot = plot_signal_len_compare_viruses,
  width = 6.91,
  height = 5.53,
  dpi = 300
)

# To compare the distributions of data from UniProt and signalpeptide.de (intersection)

###Shapiro-Wilk Test for Normality ####
### If the test is non-significant (p>. 05) it tells us that the distribution of 
### the sample is not significantly different from a normal distribution. 
### If, however, the test is significant (p < . 05) then the distribution in 
### question is significantly different from a normal distribution.


shapiro.test(na.omit(signalpeptide_meta_merged_viruses_uniprot$Length)) ### Length from signalpeptide.de (p-value < 2.2e-16 not normal distribuion) 
shapiro.test(na.omit(signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y)) ### Length from UniProt (p-value < 2.2e-16 not normal distribuion)
hist(na.omit(signalpeptide_meta_merged_viruses_uniprot$Length)) ### Distribution of the lengths in signalpaptide.de
hist(na.omit(signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y))### Distribution of the lengths in UniProt
mean(na.omit(signalpeptide_meta_merged_viruses_uniprot$Length)) ### 22.69111
mean(na.omit(signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y)) ### 22.71127
length(na.omit(signalpeptide_meta_merged_viruses_uniprot$Length)) ### 1091  
length(na.omit(signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y))  ### 1091   
wilcox.test(signalpeptide_meta_merged_viruses_uniprot$Length, signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y, paired = TRUE) ### No differences, p-value = 0.2844

### Boxplot (intersection between signalpeptide.de and UniProt (signalpeptide annotated))


png(
  filename = "results/figures/distribution_length_boxplot_viruses_intersection.png",
  width = 6.43,
  height = 6.47,
  units = "in",
  res = 300
)

boxplot(
  signalpeptide_meta_merged_viruses_uniprot$Length,
  signalpeptide_meta_merged_viruses_uniprot$Pos_sp_end_uniprot_1.y,
  names = c("signalpeptide.de", "UniProt"),
  ylab = "Length",
  main = "Distribution comparison in Viruses",
  col = c("#0073C2FF", "pink")
)

dev.off()

#### ECO evidences in Virus dataset from UniProt ###############################

eco_uni_virus <- uniprot_viruses_ft_sign_metadata$Evidence_signal_peptide_1

eco_uni_virus <- as.character(eco_uni_virus)
eco_uni_virus[is.na(eco_uni_virus)] <- "NA"

counts_eco_uni_virus <- sort(table(eco_uni_virus))

png(
  filename = "results/figures/counts_eco_uni_virus.png",
  width = 8.99,
  height = 5.55,
  units = "in",
  res = 300
)

par(mar = c(5, 10, 4, 2), pty = "s") 

bp_counts_eco_uni_virus <- barplot(
  counts_eco_uni_virus,
  las = 1,
  horiz = TRUE,
  xlab = "Count",
  main = "UniProt signal peptide annotated in Viruses (n = 1579), evidences",
  col = c("lightblue", "pink", "white", "darkgreen"),
  xlim = c(0, 2000)
)

text(
  x = counts_eco_uni_virus + max(counts_eco_uni_virus) * 0.01,
  y = bp_counts_eco_uni_virus,
  labels = counts_eco_uni_virus,
  cex = 0.8,
  pos = 4
)

dev.off()



####################### Bacteria ###############################################

#### Files of input 

# Read the file from UniProt
uniprot_bacteria_ft_sign_exp_metadata <- read_csv(
  "data/processed/uniprot/uniprot_bacteria_signal_confirmed_exp_before_update_metadata.csv"
)

# Read confirmed signal peptides from signalpeptide.de
signalpeptide_bacteria_confirmed <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_bacteria_confirmed.csv"
)
signalpeptide_bacteria_confirmed$Length <- as.numeric(
  signalpeptide_bacteria_confirmed$Length
)

# signalpeptide.de proteins with metadata retrieved from UniProt
signalpeptide_bacteria_confirmed_uniprot_meta <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_bacteria_uniprot_meta.csv"
)

# Filter proteins without a valid PrimaryAccession in UniProt
signalpeptide_bacteria_confirmed_uniprot_meta_filtered <-
  signalpeptide_bacteria_confirmed_uniprot_meta[
    !is.na(signalpeptide_bacteria_confirmed_uniprot_meta$PrimaryAccession),
  ]

#### Intersection between signalpeptide.de et UniProt

# Merge signalpeptide_bacteria and signalpeptide_bacteria_uniprot_meta 

signalpeptide_meta_merged_bacteria <- merge(signalpeptide_bacteria_confirmed_uniprot_meta_filtered, signalpeptide_bacteria_confirmed, by.x = "ProteinId_input",
                                           by.y = "Accession Number")

signalpeptide_meta_merged_bacteria_uniprot <- merge(signalpeptide_meta_merged_bacteria, uniprot_bacteria_ft_sign_exp_metadata, by = "PrimaryAccession") 

# Venn diagram

signalpeptide_uniprot_meta_merged_list_bacteria <- list(
  "signalpeptide.de with 
  UniProt metadata (n = 1139)" = signalpeptide_bacteria_confirmed_uniprot_meta_filtered$PrimaryAccession,
  "UniProt with experimentally verified
  signal peptide (n = 824)" = uniprot_bacteria_ft_sign_exp_metadata$PrimaryAccession)

venn_signalpeptide_uniprot_meta_merged_bacteria <- ggvenn(signalpeptide_uniprot_meta_merged_list_bacteria, 
                                                         fill_color = c("#0073C2FF", "pink"),
                                                         set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalpeptide_uniprot_meta_merged_bacteria <- venn_signalpeptide_uniprot_meta_merged_bacteria + 
  labs(title = "Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalpeptide_uniprot_meta_merged_bacteria)

ggsave(
  filename = "results/figures/venn_signalpeptide_uniprot_meta_merged_bacteria.png",
  plot = venn_signalpeptide_uniprot_meta_merged_bacteria,
  width = 7.79,
  height = 5.70,
  dpi = 300
)

# Statistics

# signalpeptide.de length vs uniprot with annotated signal peptide
# Intersection (N=646)

plot_signal_len_compare_bacteria <- ggplot(
  signalpeptide_meta_merged_bacteria_uniprot,
  aes(x = Length, y = Pos_sp_end_uniprot_1.y, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in signalpeptide.de") +
  ylab("Length in UniProt") +
  scale_color_viridis() +
  labs(title = "Signal peptide lengths in Bacteria (N = 646)") +
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
  scale_x_continuous(limits = c(0, 70), breaks = seq (0, 70, by = 5)) +
  scale_y_continuous(limits = c(0, 70), breaks = seq (0, 70, by = 5)) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  coord_fixed(ratio = 1)


plot_signal_len_compare_bacteria

ggsave(
  filename = "results/figures/plot_signal_len_compare_bacteria.png",
  plot = plot_signal_len_compare_bacteria,
  width = 6.91,
  height = 5.53,
  dpi = 300
)



# To compare the distributions of data from UniProt and signalpeptide.de (intersection)

###Shapiro-Wilk Test for Normality ####
### If the test is non-significant (p>. 05) it tells us that the distribution of 
### the sample is not significantly different from a normal distribution. 
### If, however, the test is significant (p < . 05) then the distribution in 
### question is significantly different from a normal distribution.


shapiro.test(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Length)) ### Length from signalpeptide.de (p-value < 2.2e-16 not normal distribuion) 
shapiro.test(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Pos_sp_end_uniprot_1.y)) ### Length from UniProt (p-value < 2.2e-16 not normal distribuion)
hist(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Length)) ### Distribution of the lengths in signalpaptide.de
hist(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Pos_sp_end_uniprot_1.y))### Distribution of the lengths in UniProt
mean(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Length)) ### 27.89474
mean(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Pos_sp_end_uniprot_1.y)) ### 27.89628
length(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Length)) ### 646  
length(na.omit(signalpeptide_meta_merged_bacteria_uniprot$Pos_sp_end_uniprot_1.y))  ### 646   
wilcox.test(signalpeptide_meta_merged_bacteria_uniprot$Length, signalpeptide_meta_merged_bacteria_uniprot$Pos_sp_end_uniprot_1.y, paired = TRUE) ### No differences, p-value = 1; most paired differences are zero

### Boxplot (intersection between signalpeptide.de and UniProt (signalpeptide annotated))


png(
  filename = "results/figures/distribution_length_boxplot_bacteria_intersection.png",
  width = 6.43,
  height = 6.47,
  units = "in",
  res = 300
)

boxplot(
  signalpeptide_meta_merged_bacteria_uniprot$Length,
  signalpeptide_meta_merged_bacteria_uniprot$Pos_sp_end_uniprot_1.y,
  names = c("signalpeptide.de", "UniProt"),
  ylab = "Length",
  main = "Distribution comparison in Bacteria",
  col = c("#0073C2FF", "pink")
)

dev.off()

#### ECO evidences in Bacteria dataset from UniProt ###############################

eco_uni_bacteria <- uniprot_bacteria_ft_sign_exp_metadata$Evidence_signal_peptide_1
  
eco_uni_bacteria <- as.character(eco_uni_bacteria)
eco_uni_bacteria[is.na(eco_uni_bacteria)] <- "NA"

counts_eco_uni_bacteria <- sort(table(eco_uni_bacteria))

png(
  filename = "results/figures/counts_eco_uni_bacteria.png",
  width = 8.99,
  height = 5.55,
  units = "in",
  res = 300
)

par(mar = c(5, 10, 4, 2), pty = "s") 

bp_counts_eco_uni_bacteria <- barplot(
  counts_eco_uni_bacteria,
  las = 1,
  horiz = TRUE,
  xlab = "Count",
  main = "UniProt signal peptide annotated in Bacteria (n = 824), evidences",
  col = c("pink", "lightblue","white"),
  xlim = c(0, 1000)
)

text(
  x = counts_eco_uni_bacteria + max(counts_eco_uni_bacteria) * 0.01,
  y = bp_counts_eco_uni_bacteria,
  labels = counts_eco_uni_bacteria,
  cex = 0.8,
  pos = 4
)

dev.off()



####################### Mammalia ###############################################


#### Files of input 

# Read the file from UniProt
uniprot_mammalia_ft_sign_exp_metadata <- read_csv(
  "data/processed/uniprot/uniprot_mammalia_signal_confirmed_exp_before_update_metadata.csv"
)

# Read confirmed signal peptides from signalpeptide.de
signalpeptide_mammalia_confirmed <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_mammalia_confirmed.csv"
)
signalpeptide_mammalia_confirmed$Length <- as.numeric(
  signalpeptide_mammalia_confirmed$Length
)

# signalpeptide.de proteins with metadata retrieved from UniProt
signalpeptide_mammalia_confirmed_uniprot_meta <- read_csv(
  "data/processed/signalpeptide_de/signalpeptide_mammalia_uniprot_meta.csv"
)

### To filter the lines with NA in Primary Accession number in UniProt (this protein doesn't exist)
signalpeptide_mammalia_confirmed_uniprot_meta_filtered <- signalpeptide_mammalia_confirmed_uniprot_meta[!is.na(signalpeptide_mammalia_confirmed_uniprot_meta$PrimaryAccession), ]

#### Intersection between signalpeptide.de et UniProt

# Merge signalpeptide_mammalia and signalpeptide_mammalia_uniprot_meta 

signalpeptide_meta_merged_mammalia <- merge(
  signalpeptide_mammalia_confirmed_uniprot_meta_filtered,
  signalpeptide_mammalia_confirmed,
  by.x = "ProteinId_input",
  by.y = "ProteinId"
)

signalpeptide_meta_merged_mammalia_uniprot <- merge(signalpeptide_meta_merged_mammalia, uniprot_mammalia_ft_sign_exp_metadata, by = "PrimaryAccession") 


# Venn diagram

signalpeptide_uniprot_meta_merged_list_mammalia <- list(
  "signalpeptide.de with 
  UniProt metadata (n = 2104)" = signalpeptide_mammalia_confirmed_uniprot_meta_filtered$PrimaryAccession,
  "UniProt with experimentally verified
  signal peptide (n = 1463)" = uniprot_mammalia_ft_sign_exp_metadata$PrimaryAccession)

venn_signalpeptide_uniprot_meta_merged_mammalia <- ggvenn(signalpeptide_uniprot_meta_merged_list_mammalia, 
                                                          fill_color = c("#0073C2FF", "pink"),
                                                          set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalpeptide_uniprot_meta_merged_mammalia <- venn_signalpeptide_uniprot_meta_merged_mammalia + 
  labs(title = "Mammalia") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalpeptide_uniprot_meta_merged_mammalia)

ggsave(
  filename = "results/figures/venn_signalpeptide_uniprot_meta_merged_mammalia.png",
  plot = venn_signalpeptide_uniprot_meta_merged_mammalia,
  width = 7.79,
  height = 5.70,
  dpi = 300
)

# Statistics

# signalpeptide.de length vs uniprot with annotated signal peptide
# Intersection (N=1350)

plot_signal_len_compare_mammalia <- ggplot(
  signalpeptide_meta_merged_mammalia_uniprot,
  aes(x = Length, y = Pos_sp_end_uniprot_1.y, na.rm = TRUE)
) +
  geom_point() +
  geom_pointdensity(show.legend = TRUE) +
  theme_classic() +
  scale_color_viridis() +
  xlab("Length in signalpeptide.de") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in Mammalia (N = 1350)") +
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
  scale_x_continuous(limits = c(0, 70), breaks = seq (0, 70, by = 5)) +
  scale_y_continuous(limits = c(0, 70), breaks = seq (0, 70, by = 5)) +
  stat_cor(
    method = "spearman",
    label.x = 5,
    label.y = 65
  ) +
  coord_fixed(ratio = 1)


plot_signal_len_compare_mammalia

ggsave(
  filename = "results/figures/plot_signal_len_compare_mammalia.png",
  plot = plot_signal_len_compare_mammalia,
  width = 6.91,
  height = 5.53,
  dpi = 300
)



# To compare the distributions of data from UniProt and signalpeptide.de (intersection)

###Shapiro-Wilk Test for Normality ####
### If the test is non-significant (p>. 05) it tells us that the distribution of 
### the sample is not significantly different from a normal distribution. 
### If, however, the test is significant (p < . 05) then the distribution in 
### question is significantly different from a normal distribution.


shapiro.test(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Length)) ### Length from signalpeptide.de (p-value < 2.2e-16 not normal distribuion) 
shapiro.test(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Pos_sp_end_uniprot_1.y)) ### Length from UniProt (p-value < 2.2e-16 not normal distribuion)
hist(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Length)) ### Distribution of the lengths in signalpaptide.de
hist(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Pos_sp_end_uniprot_1.y))### Distribution of the lengths in UniProt
mean(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Length)) ### 23.01627
mean(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Pos_sp_end_uniprot_1.y)) ### 23.09689
length(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Length)) ### 1352  
length(na.omit(signalpeptide_meta_merged_mammalia_uniprot$Pos_sp_end_uniprot_1.y))  ### 1352   
wilcox.test(signalpeptide_meta_merged_mammalia_uniprot$Length, signalpeptide_meta_merged_mammalia_uniprot$Pos_sp_end_uniprot_1.y, paired = TRUE) ### p-value = 0.04317; difference

### Boxplot (intersection between signalpeptide.de and UniProt (signalpeptide annotated))


png(
  filename = "results/figures/distribution_length_boxplot_mammalia_intersection.png",
  width = 6.43,
  height = 6.47,
  units = "in",
  res = 300
)

boxplot(
  signalpeptide_meta_merged_mammalia_uniprot$Length,
  signalpeptide_meta_merged_mammalia_uniprot$Pos_sp_end_uniprot_1.y,
  names = c("signalpeptide.de", "UniProt"),
  ylab = "Length",
  main = "Distribution comparison in Mammalia",
  col = c("#0073C2FF", "pink")
)

dev.off()

#### ECO evidences in Mammalia dataset from UniProt ###############################

eco_uni_mammalia <- uniprot_mammalia_ft_sign_exp_metadata$Evidence_signal_peptide_1

eco_uni_mammalia <- as.character(eco_uni_mammalia)
eco_uni_mammalia[is.na(eco_uni_mammalia)] <- "NA"

counts_eco_uni_mammalia <- sort(table(eco_uni_mammalia))

png(
  filename = "results/figures/counts_eco_uni_mammalia.png",
  width = 8.99,
  height = 5.55,
  units = "in",
  res = 300
)

par(mar = c(5, 10, 4, 2), pty = "s") 

bp_counts_eco_uni_mammalia <- barplot(
  counts_eco_uni_mammalia,
  las = 1,
  horiz = TRUE,
  xlab = "Count",
  main = "UniProt signal peptide annotated in Mammalia (n = 1463), evidences",
  col = c("pink", "lightblue","white"),
  xlim = c(0, 2000)
)

text(
  x = counts_eco_uni_mammalia + max(counts_eco_uni_mammalia) * 0.01,
  y = bp_counts_eco_uni_mammalia,
  labels = counts_eco_uni_mammalia,
  cex = 0.8,
  pos = 4
)

dev.off()

################################################################################
################ Signal peptide length by group — UniProt ######################
################################################################################

uniprot_sp_length_groups <- bind_rows(
  data.frame(
    Group = "Viruses",
    Signal_peptide_length = uniprot_viruses_ft_sign_metadata$Pos_sp_end_uniprot_1
  ),
  data.frame(
    Group = "Mammalia",
    Signal_peptide_length = uniprot_mammalia_ft_sign_exp_metadata$Pos_sp_end_uniprot_1
  ),
  data.frame(
    Group = "Bacteria",
    Signal_peptide_length = uniprot_bacteria_ft_sign_exp_metadata$Pos_sp_end_uniprot_1
  )
)

uniprot_sp_length_groups <- uniprot_sp_length_groups %>%
  filter(!is.na(Signal_peptide_length))

uniprot_sp_length_groups$Group <- factor(
  uniprot_sp_length_groups$Group,
  levels = c("Bacteria", "Mammalia", "Viruses")
)


plot_uniprot_all_group_sp_length <- ggplot(
  uniprot_sp_length_groups,
  aes(
    x = Group,
    y = Signal_peptide_length,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.75,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 140),
    breaks = seq(0, 140, by = 10)
  ) +
  scale_fill_manual(
    values = c(
      "Viruses" = "#619CFF",
      "Mammalia" = "#00BA38",
      "Bacteria" = "#F8766D"
    )
  ) +
  labs(
    title = "Signal peptide length by group (UniProt Database)",
    x = NULL,
    y = "Signal peptide length"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13)
  )

plot_uniprot_all_group_sp_length


ggsave(
  filename = "results/figures/uniprot_all_group_sp_length.png",
  plot = plot_uniprot_all_group_sp_length,
  width = 10.5,
  height = 7,
  dpi = 300
)


### French version
plot_uniprot_all_group_sp_length_fr <- ggplot(
  uniprot_sp_length_groups,
  aes(
    x = Group,
    y = Signal_peptide_length,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.75,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 140),
    breaks = seq(0, 140, by = 10)
  ) +
  scale_x_discrete(
    labels = c(
      "Viruses" = "Virus",
      "Mammalia" = "Mammifères",
      "Bacteria" = "Bactéries"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Viruses" = "#619CFF",
      "Mammalia" = "#00BA38",
      "Bacteria" = "#F8766D"
    )
  ) +
  labs(
    title = "Longueur des peptides signaux par groupe (base de données UniProt)",
    x = NULL,
    y = "Longueur du peptide signal"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13)
  )

plot_uniprot_all_group_sp_length_fr


ggsave(
  filename = "results/figures/uniprot_all_group_sp_length_fr.png",
  plot = plot_uniprot_all_group_sp_length_fr,
  width = 10.5,
  height = 7,
  dpi = 300
)

#####################################################################################
####################### Intersection with SignalP training set ######################
#####################################################################################

# Read the file from SignalP training set (metadata form UniProt)
protein_ids_dataset_signalP_metadata <- read_csv(
  "data/processed/signalp_training_set/protein_ids_dataset_signalP_training_set_metadata.csv"
)

# Venn diagram / SignalP training set vs UniProt viruses ######################

signalp_trainset_uniprot_list_viruses <- list(
  "SignalP training set (n = 20263)" = protein_ids_dataset_signalP_metadata$PrimaryAccession,
  "UniProt with experimentally verified
  signal peptide in Virus (n = 1579)" = uniprot_viruses_ft_sign_metadata$PrimaryAccession)

venn_signalp_trainset_uniprot_list_viruses <- ggvenn(signalp_trainset_uniprot_list_viruses, 
                                                          fill_color = c("#0073C2FF", "pink"),
                                                          set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalp_trainset_uniprot_list_viruses <- venn_signalp_trainset_uniprot_list_viruses + 
  labs(title = "SignalP training set vs Uniprot Viruses") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalp_trainset_uniprot_list_viruses)


ggsave(
  filename = "results/figures/venn_signalp_trainset_uniprot_list_viruses.png",
  plot = venn_signalp_trainset_uniprot_list_viruses,
  width = 10.13,
  height = 6.67,
  dpi = 300
)

# Venn diagram / SignalP training set vs UniProt mammalia######################


signalp_trainset_uniprot_list_mammalia <- list(
  "SignalP training set (n = 20263)" = protein_ids_dataset_signalP_metadata$PrimaryAccession,
  "UniProt with experimentally verified
  signal peptide in Mammalia (n = 1463)" = uniprot_mammalia_ft_sign_exp_metadata$PrimaryAccession)

venn_signalp_trainset_uniprot_list_mammalia <- ggvenn(signalp_trainset_uniprot_list_mammalia, 
                                                     fill_color = c("#0073C2FF", "pink"),
                                                     set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalp_trainset_uniprot_list_mammalia <- venn_signalp_trainset_uniprot_list_mammalia + 
  labs(title = "SignalP training set vs Uniprot Mammalia") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalp_trainset_uniprot_list_mammalia)

ggsave(
  filename = "results/figures/venn_signalp_trainset_uniprot_list_mammalia.png",
  plot = venn_signalp_trainset_uniprot_list_mammalia,
  width = 10.13,
  height = 6.67,
  dpi = 300
)

# Venn diagram / SignalP training set vs UniProt bacteria ######################


signalp_trainset_uniprot_list_bacteria <- list(
  "SignalP training set (n = 20263)" = protein_ids_dataset_signalP_metadata$PrimaryAccession,
  "UniProt with experimentally verified
  signal peptide in Bacteria (n = 824)" = uniprot_bacteria_ft_sign_exp_metadata$PrimaryAccession)

venn_signalp_trainset_uniprot_list_bacteria <- ggvenn(signalp_trainset_uniprot_list_bacteria, 
                                                      fill_color = c("#0073C2FF", "pink"),
                                                      set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalp_trainset_uniprot_list_bacteria <- venn_signalp_trainset_uniprot_list_bacteria + 
  labs(title = "SignalP training set vs Uniprot Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalp_trainset_uniprot_list_bacteria)

ggsave(
  filename = "results/figures/venn_signalp_trainset_uniprot_list_bacteria.png",
  plot = venn_signalp_trainset_uniprot_list_bacteria,
  width = 10.13,
  height = 6.67,
  dpi = 300
)


# Venn diagram / SignalP training set vs signalpeptide.de Virus ######################

signalp_trainset_spde_list_viruses <- list(
  "SignalP training set (n = 20263)" = protein_ids_dataset_signalP_metadata$PrimaryAccession,
  "signalpeptide.de in Virus (n = 6198)" = signalpeptide_viruses$'Accession Number')

venn_signalp_trainset_spde_list_viruses <- ggvenn(signalp_trainset_spde_list_viruses, 
                                                     fill_color = c("#0073C2FF", "pink"),
                                                     set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalp_trainset_spde_list_viruses <- venn_signalp_trainset_spde_list_viruses + 
  labs(title = "SignalP training set vs signalpeptide.de Viruses") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalp_trainset_spde_list_viruses)


ggsave(
  filename = "results/figures/venn_signalp_trainset_spde_list_viruses.png",
  plot = venn_signalp_trainset_spde_list_viruses,
  width = 10.13,
  height = 6.67,
  dpi = 300
)

# Venn diagram / SignalP training set vs signalpeptide.de Mammalia ######################

signalp_trainset_spde_list_mammalia <- list(
  "SignalP training set (n = 20263)" = protein_ids_dataset_signalP_metadata$PrimaryAccession,
  "signalpeptide.de in Mammalia (cfrm = 2109)" = signalpeptide_mammalia_confirmed$ProteinId)

venn_signalp_trainset_spde_list_mammalia <- ggvenn(signalp_trainset_spde_list_mammalia, 
                                                  fill_color = c("#0073C2FF", "pink"),
                                                  set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalp_trainset_spde_list_mammalia <- venn_signalp_trainset_spde_list_mammalia + 
  labs(title = "SignalP training set vs signalpeptide.de Mammalia") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalp_trainset_spde_list_mammalia)


ggsave(
  filename = "results/figures/venn_signalp_trainset_spde_list_mammalia.png",
  plot = venn_signalp_trainset_spde_list_mammalia,
  width = 10.13,
  height = 6.67,
  dpi = 300
)



# Venn diagram / SignalP training set vs signalpeptide.de Bacteria ######################

signalp_trainset_spde_list_bacteria <- list(
  "SignalP training set (n = 20263)" = protein_ids_dataset_signalP_metadata$PrimaryAccession,
  "signalpeptide.de in Bacteria (cfrm = 1161)" = signalpeptide_bacteria_confirmed$'Accession Number')

venn_signalp_trainset_spde_list_bacteria <- ggvenn(signalp_trainset_spde_list_bacteria, 
                                                   fill_color = c("#0073C2FF", "pink"),
                                                   set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalp_trainset_spde_list_bacteria <- venn_signalp_trainset_spde_list_bacteria + 
  labs(title = "SignalP training set vs signalpeptide.de Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

print(venn_signalp_trainset_spde_list_bacteria)


ggsave(
  filename = "results/figures/venn_signalp_trainset_spde_list_bacteria.png",
  plot = venn_signalp_trainset_spde_list_bacteria,
  width = 10.13,
  height = 6.67,
  dpi = 300
)

######## Evidences in SignalP training dataset ################################# 

x <- protein_ids_dataset_signalP_metadata$Evidence_signal_peptide_1

x <- as.character(x)
x[is.na(x)] <- "NA"

counts_ECO <- sort(table(x))

png(
  filename = "results/figures/protein_ids_dataset_signalP_metadata.png",
  width = 8.99,
  height = 5.55,
  units = "in",
  res = 300
)

par(mar = c(5, 10, 4, 2), pty = "s") 

bp_trainsst_signalp <- barplot(
  counts_ECO,
  las = 1,
  horiz = TRUE,
  xlab = "Count",
  main = "SignalP training dataset (n = 20290), evidences",
  col = c("lightblue", "pink", "white", "darkgreen"),
  xlim = c(0, 20000)
)

text(
  x = counts_ECO + max(counts_ECO) * 0.01,
  y = bp_trainsst_signalp,
  labels = counts_ECO,
  cex = 0.8,
  pos = 4
)

dev.off()

png(
  filename = "results/figures/protein_ids_dataset_signalP_metadata_fr.png",
  width = 8.99,
  height = 5.55,
  units = "in",
  res = 300
)

par(mar = c(5, 10, 4, 2), pty = "s") 

bp_trainsst_signalp_fr <- barplot(
  counts_ECO,
  las = 1,
  horiz = TRUE,
  xlab = "Nombre",
  main = "Jeu de données d’entraînement de SignalP (n = 20 290) : types de preuves",
  col = c("lightblue", "pink", "white", "darkgreen"),
  xlim = c(0, 20000)
)

text(
  x = counts_ECO + max(counts_ECO) * 0.01,
  y = bp_trainsst_signalp,
  labels = counts_ECO,
  cex = 0.8,
  pos = 4
)

dev.off()


