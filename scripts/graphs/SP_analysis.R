getwd()
setwd("/Users/tatianakrylova/Documents/STAGE_LIRMM/Travail/data/R_analysis")
getwd()

# Necessary libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggvenn)
library(ggpubr)

################################################################################
############################## Viruses #########################################
################################################################################


# Read the file from SignalP 6.0 results ("output.gff3")
signalP6_table_viruses <- read.table(
  "viruses/output_signalp6_all_sp.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

# Select columns from signalP6_table
signalP6_table_viruses_an <- signalP6_table_viruses[0:6]


# Add the headers 
colnames(signalP6_table_viruses_an) <- c("ProteinId", "App", "Type", "Pos_start", "Pos_end", "Score")

# Download the tables for megre

# Table from signalpeptide.de
signalpeptide_viruses <- read.table(
  "viruses/signalpeptide_viruses.csv",
  sep = ",",
  header = TRUE,
  quote = "\"",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

signalpeptide_viruses$Length <- as.numeric(signalpeptide_viruses$Length)

mean_by_group <- signalpeptide_viruses %>%
  group_by(SP.Status) %>%
  summarise(
    Length = mean(Length, na.rm = TRUE)
  )



# List of all ProteinId which present in Uniprot
fasta_full_peptides_viruses <- readLines("viruses/peptides_full.fasta")
uniprot_ids_viruses <- fasta_full_peptides_viruses[grepl("^>", fasta_full_peptides_viruses)]
uniprot_ids_viruses <- sub("^>", "", uniprot_ids_viruses)
length(uniprot_ids_viruses) ### intersection between signalpeptide.de proteins and UniProt (actual information)

str(uniprot_viruses_sppr$SPuniprot)
summary(uniprot_viruses_sppr$SPuniprot)


names(uniprot_viruses_sppr)

par(mar = c(2, 2, 2, 1))
boxplot(uniprot_viruses_sppr[["SPuniprot"]],
        ylab = "Signal peptide length")


##### To calculate ECO in UniProt data viruses #################################
str(uniprot_viruses_sppr$Signal.peptide)
library(tidyverse)

tibble(raw = uniprot_viruses_sppr$Signal.peptide) %>%
  mutate(ECO = str_extract_all(raw, "ECO:\\d+")) %>%
  mutate(ECO = map(ECO, ~ if (length(.x) == 0) NA_character_ else .x)) %>%
  unnest(ECO) %>%
  count(ECO, sort = TRUE)


library(stringr)

uniprot_viruses_sppr %>%
  mutate(ECO = str_extract_all(`Signal.peptide`, "ECO:\\d+")) %>%
  mutate(ECO = map(ECO, ~ if (length(.x) == 0) NA_character_ else .x)) %>%
  unnest(ECO) %>%
  distinct(`Accession.Number`, ECO) %>%
  count(ECO, sort = TRUE)


################################################################################



## Merge output from signalpeptide.de and output from SignalP 6.0
spP6_spSite_merged <- merge(signalP6_table_viruses_an, signalpeptide_viruses, by.x = "ProteinId",
                     by.y = "Accession.Number", all = TRUE)

## Delete the lignes where ProteinId = '<<' # 124 proteins dropped
spP6_spSite_merged <- spP6_spSite_merged %>%
  filter(ProteinId != "<<")


## New column yes/no about the presence of UniProtID in modern UniProt Database (compare to table from site)
spP6_spSite_merged$ActualUniProt <- ifelse(spP6_spSite_merged$ProteinId %in% uniprot_ids_viruses, "yes", "no" )

## Number of SP which presents in Uniprot but was refused by SignalP6
sum(is.na(spP6_spSite_merged$App) & spP6_spSite_merged$ActualUniProt == "yes") #### 621

## Number of SP which presents in Uniprot and has a status confirmed (in the site signalprptide.de) but was refused by SignalP6 ### 0
sum(is.na(spP6_spSite_merged$App) & spP6_spSite_merged$ActualUniProt == "yes" & spP6_spSite_merged$`SP Status` == "confirmed")

## Number of SP which presents in Uniprot and with the status confirmed in signalpeptide database (N=0)
sum(spP6_spSite_merged$'SP Status' == 'confirmed' & spP6_spSite_merged$ActualUniProt == "yes")

spP6_spSite_merged$Length <- as.numeric(spP6_spSite_merged$Length)


## Graph: correlation between the length from signalpeptide.de and SigalP 6 prediction
ggplot(aes(x=spP6_spSite_merged$Pos_end, y=spP6_spSite_merged$Length), data=spP6_spSite_merged) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("Position predicted by SignalP6") +
  ylab("Position in signalpeptide map") +
  theme_classic() +
  ylim(0,65) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)


## Table downloaded from UniProt with confirmed presence of signal peptide (viruses) - bash command via query Uniprot
uniprot_viruses_signal_confirmed <- read.table("viruses/uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11.tsv", sep = "\t",
                                               header = TRUE,
                                               quote = "",
                                               fill = TRUE,
                                               comment.char = "",
                                               stringsAsFactors = FALSE)
uniprot_viruses_signal_confirmed$SPuniprot<- sub(".*\\.\\.(\\d+).*",
                                                 "\\1", uniprot_viruses_signal_confirmed$`Signal.peptide`) ## extraction of end position from UniProt
uniprot_viruses_signal_confirmed$SPuniprot <- as.numeric(uniprot_viruses_signal_confirmed$SPuniprot)

mean(uniprot_viruses_signal_confirmed$SPuniprot, na.rm = TRUE)



### Venn graph for all datasets (viruses)


signalpeptide_viruses_ids_clean <- signalpeptide_viruses$Accession.Number[
  !is.na(signalpeptide_viruses$Accession.Number) &
    !grepl("<<", signalpeptide_viruses$Accession.Number)
]

SignalPeptides_viruses_db_sp_all <- list(
  "SignalP (n = 5553)" = signalP6_table_viruses_an$ProteinId,
  "Signal Peptide DB (n = 6197)" = signalpeptide_viruses_ids_clean,
  "UniProt (n = 1579)" = uniprot_viruses_signal_confirmed$Entry) ### presence of siignal peptide confirmed
len_db_sp_venn_vir_all <- ggvenn(SignalPeptides_viruses_db_sp_all, 
                             fill_color = c("#0073C2FF", "pink", "white"),
                             set_name_size = 5,text_size = 5,show_percentage = FALSE)

len_db_sp_venn_vir_all

ggsave(
  filename = "../../../Images/venn_signalpeptides_viruses_all.png",
  plot = len_db_sp_venn_vir_all,
  width = 7.55,
  height = 6.79,
  dpi = 300
)

library(ggvenn)
SignalPeptides_viruses_db_sp_venn <- list(
  "Signal Peptide DB " = signalpeptide_viruses_ids_clean,
  "UniProt " = uniprot_viruses_sppr$Accession.Number)
len_db_sp_venn_vir_db_sp2 <- ggvenn(SignalPeptides_viruses_db_sp_venn, 
                                 fill_color = c("#0073C2FF", "pink"),
                                 set_name_size = 5,text_size = 5,show_percentage = FALSE)

len_db_sp_venn_vir_db_sp2
ggsave(
  filename = "../../../Images/len_db_sp_venn_vir_db_sp2.png",
  plot = len_db_sp_venn_vir_db_sp2,
  width = 6.63,
  height = 4.46,
  dpi = 300
)

## Merge between table from UniProt and spP6_spSite_merged 
spP6_spSite_merged_uniprot_vir <- merge(spP6_spSite_merged, uniprot_viruses_signal_confirmed, by.x = "ProteinId",
                            by.y = "Entry", all = TRUE)
spP6_spSite_merged_uniprot_vir

spP6_spSite_merged_uniprot_vir$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", spP6_spSite_merged_uniprot_vir$`Signal.peptide`) ## extraction of end position from UniProt
spP6_spSite_merged_uniprot_vir$SPuniprot <- as.numeric(spP6_spSite_merged_uniprot_vir$SPuniprot)


length(spP6_spSite_merged_uniprot_vir$Pos_end) ### 6698
length(spP6_spSite_merged_uniprot_vir$SPuniprot) ### 6698

### Verification of the presence of NA
length(na.omit(spP6_spSite_merged_uniprot_vir$Pos_end)) ### 5553
length(na.omit(spP6_spSite_merged_uniprot_vir$SPuniprot)) #### 1524

length( uniprot_viruses_signal_confirmed$Entry) ### 1579
length(spP6_spSite_merged_uniprot_vir$ProteinId) ### 6698



#### Graph Length in SignalP vs Length in UniProt


signal_peptide_lengths_viruses_all <- ggplot(
  spP6_spSite_merged_uniprot_vir,
  aes(x = Pos_end, y = SPuniprot)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in Viruses") +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 68
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )


signal_peptide_lengths_viruses_all


ggsave(
  filename = "../../../Images/signal_peptide_lengths_viruses_all.png",
  plot = signal_peptide_lengths_viruses_all,
  width = 7.48,
  height = 4.85,
  dpi = 300
)


### intersection between SignalP results and UniProt
merged_spP6_uniprot_vir_intersection <- merge(signalP6_table_viruses_an, uniprot_viruses_signal_confirmed, by.x = "ProteinId",
                                        by.y = "Entry") 



#### Graph Length in SignalP vs Length in UniProt - intersection
signal_peptide_lengths_viruses_intersection <- ggplot(merged_spP6_uniprot_vir_intersection, aes ( x = merged_spP6_uniprot_vir_intersection$Pos_end, y = merged_spP6_uniprot_vir_intersection$SPuniprot)) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in Viruses (intersection, n=875)") +
  geom_smooth(method = 'lm',
              se = TRUE, na.rm = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 68
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )


signal_peptide_lengths_viruses_intersection

ggsave(
  filename = "../../../Images/signal_peptide_lengths_viruses_intersection.png",
  plot = signal_peptide_lengths_viruses_intersection,
  width = 7.48,
  height = 4.85,
  dpi = 300
)



### Boxplots for viruses - intersection


png(
  filename = "../../../Images/boxplot_viruses_signalp_uniprot_intersection.png",
  width = 6.40,
  height = 6.53,
  units = "in",
  res = 300
)

ymax <- max(c(merged_spP6_uniprot_vir_intersection$Pos_end, merged_spP6_uniprot_vir_intersection$SPuniprot), na.rm = TRUE) * 1.25

boxplot(list(Pos_end = (na.omit(merged_spP6_uniprot_vir_intersection$Pos_end)),
             SPuniprot = (na.omit(merged_spP6_uniprot_vir_intersection$SPuniprot))),
        ylab = "Length of SP",
        main = "Signal peptide lengths in Viruses (SignalP vs UniProt, n = 875)", ###intersection
        names = c("SignalP prediction", "UniProt"), col = c("skyblue", "pink"),
        boxwex = 0.4,
        ylim = c(0, ymax)
)


y <- max(c(merged_spP6_uniprot_vir_intersection$Pos_end, merged_spP6_uniprot_vir_intersection$SPuniprot), na.rm = TRUE) * 1.08

segments(1, y, 2, y)
segments(1, y, 1, y - 3)
segments(2, y, 2, y - 3)

text(
  x = 1.5,
  y = y + 3.5,
  labels = "***",
  cex = 1.5
)

dev.off()



## Statistics

###Shapiro-Wilk Test for Normality ####
### If the test is non-significant (p>. 05) it tells us that the distribution of 
### the sample is not significantly different from a normal distribution. 
### If, however, the test is significant (p < . 05) then the distribution in 
### question is significantly different from a normal distribution.



shapiro.test(na.omit(merged_spP6_uniprot_vir_intersection$Pos_end)) ### p-value < 2.2e-16 not normal distribuion
shapiro.test(na.omit(merged_spP6_uniprot_vir_intersection$SPuniprot)) ### p-value < 2.2e-16 not normal distribuion
hist(na.omit(merged_spP6_uniprot_vir_intersection$Pos_end))
hist(na.omit(merged_spP6_uniprot_vir_intersection$SPuniprot))
length(na.omit(merged_spP6_uniprot_vir_intersection$Pos_end)) #### 875 
length(na.omit(merged_spP6_uniprot_vir_intersection$SPuniprot)) ### 875

merged_spP6_uniprot_vir_intersection

wilcox.test(merged_spP6_uniprot_vir_intersection$Pos_end, merged_spP6_uniprot_vir_intersection$SPuniprot, paired = TRUE) # statistically different p-value = 8.458e-10


summary(merged_spP6_uniprot_vir_intersection$Pos_end) ### Mean 21.34
summary(merged_spP6_uniprot_vir_intersection$SPuniprot) ### Mean 21.07


### FALSE and TRUE prediction of the SP length

comparison_viruses <- merged_spP6_uniprot_vir_intersection %>%
  filter(!is.na(Pos_end), !is.na(SPuniprot)) %>%
  mutate(
    match = Pos_end == SPuniprot
  )

table(comparison_viruses$match) ### FALSE = 268 (30.6%), TRUE = 607 (69.3%)


###### The difference between the predicted length by SignalP 6 and actual information in UniProt

spP6_spSite_merged_uniprot_vir$diff = spP6_spSite_merged_uniprot_vir$Pos_end/spP6_spSite_merged_uniprot_vir$SPuniprot



png(
  filename = "../../../Images/hist_diff_viruses.png",
  width = 604,
  height = 533
)

hist(spP6_spSite_merged_uniprot_vir$diff, main = "The difference between SignalP and UniProt
     (ratio of length, viruses)", xlab = "Ratio")

dev.off()


### Histograms for precision between 0.9 and 1.0
hist(
  spP6_spSite_merged_uniprot_vir$diff[
    spP6_spSite_merged_uniprot_vir$diff >= 0.9 &
      spP6_spSite_merged_uniprot_vir$diff <= 1.0
  ],
  main = "The difference between SignalP and UniProt\n(ratio of length, viruses)",
  xlab = "Ratio",
  breaks = seq(0.9, 1.0, by = 0.01)
)


png(
  filename = "../../../Images/hist_diff_viruses_small.png",
  width = 604,
  height = 533
)


hist(
  spP6_spSite_merged_uniprot_vir$diff[
    spP6_spSite_merged_uniprot_vir$diff >= 0.99 &
      spP6_spSite_merged_uniprot_vir$diff <= 1.0
  ],
  main = "The difference between SignalP and UniProt\n(ratio of length, viruses)",
  xlab = "Ratio",
  breaks = seq(0.99, 1.0, by = 0.001)
)

dev.off()


# Table, which contains ProteinId, length from the site signalpeptide.de, name and !!!taxonId of the host!!!
hosts_taxID <- read.table("viruses/hosts_taxID.csv", 
                          sep = ",",
                          header = TRUE)

### Table for make a boxplot of length per Organism Host
hosts_long_viruses <- hosts_taxID %>%
  separate_rows(Host, TaxonID, sep = ",\\s*")



################################################################################
############################## Mammalia ########################################
################################################################################

# Read the file from SignalP 6.0 results ("output.gff3")
signalP6_table_mammalia <- read.table("mammalia/output_mammalia_confirmed.gff3", sep = "\t", header = FALSE, 
                                      comment.char = "#", stringsAsFactors = FALSE )
colnames(signalP6_table_mammalia) <- c("ProteinId", "App", "Type", "Pos_start", "Pos_end", "Score")
signalP6_table_mammalia_an <- signalP6_table_mammalia[0:6]
hist(signalP6_table_mammalia_an$Pos_end)
signalP6_table_mammalia_an$ProteinID <- sub("\\s*\\|.*$", "", signalP6_table_mammalia_an$ProteinId)
signalP6_table_mammalia_an$ProteinId <- NULL

# Read the file from signalpeptide.de ("signalpeptide_mammalia_confirmed.csv"), confirmed
signalpeptide_mammalia_confirmed <- read.csv(
  "mammalia/signalpeptide_mammalia_confirmed.csv",
  quote = "\"",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)


mean_by_group_mammalia <- signalpeptide_mammalia_confirmed %>%
  group_by(SP.Status) %>%
  summarise(
    Length = mean(Length, na.rm = TRUE)
  )

mean_by_group_mammalia


## Merge output from signalpeptide.de (confirmed) and output from SignalP 6.0
merged_sp_mammalia <- merge(signalP6_table_mammalia_an, signalpeptide_mammalia_confirmed, by.x = "ProteinID", by.y = 'Accession.Number')
merged_sp_mammalia$Length <- as.numeric(merged_sp_mammalia$Length)

## Graph Signal peptide lengths in Mammalia (signalpeptide.de vs SignalP6)
ggplot(merged_sp_mammalia, aes ( x = merged_sp_mammalia$Pos_end, y = merged_sp_mammalia$Length)) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP6") +
  ylab("Length in DB") +
  labs(title = "Signal peptide lengths in Mammalia") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 68
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )



## Table downloaded from UniProt with confirmed presence of signal peptide (mammalia) - bash command via query Uniprot / only with evidence ECO:0000269 !!!!!
uniprot_mammalia_signal_confirmed_exp <- read_tsv(
  "mammalia/uniprot_mammalia_signal_confirmed_exp.tsv",
  show_col_types = FALSE)


uniprot_mammalia_signal_confirmed_exp$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", uniprot_mammalia_signal_confirmed_exp$`Signal peptide`) ## extraction of end position from UniProt
uniprot_mammalia_signal_confirmed_exp$SPuniprot <- as.numeric(uniprot_mammalia_signal_confirmed_exp$SPuniprot)

signalpeptide_mammalia_ids_clean <- signalpeptide_mammalia_confirmed$Accession.Number[
  !is.na(signalpeptide_mammalia_confirmed$Accession.Number) &
    !grepl("<<", signalpeptide_mammalia_confirmed$Accession.Number)
]

SignalPeptides_mammalia_db_sp_all <- list(
  "Signal Peptide DB " = signalpeptide_mammalia_ids_clean,
  "UniProt " = uniprot_mammalia_signal_confirmed_exp$Entry) ### presence of siignal peptide confirmed
len_db_sp_venn_mamm_uniprot <- ggvenn(SignalPeptides_mammalia_db_sp_all, 
                                 fill_color = c("#0073C2FF", "pink", "white"),
                                 set_name_size = 5,text_size = 5,show_percentage = FALSE)

len_db_sp_venn_mamm_uniprot

ggsave(
  filename = "../../../Images/len_db_sp_venn_mamm_uniprot.png",
  plot = len_db_sp_venn_mamm_uniprot,
  width = 7.55,
  height = 6.79,
  dpi = 300
)















##### To calculate ECO in UniProt data mammalia #################################
str(uniprot_viruses_sppr$Signal.peptide)
library(tidyverse)

tibble(raw = uniprot_mammalia_signal_confirmed_exp$`Signal peptide`) %>%
  mutate(ECO = str_extract_all(raw, "ECO:\\d+")) %>%
  mutate(ECO = map(ECO, ~ if (length(.x) == 0) NA_character_ else .x)) %>%
  unnest(ECO) %>%
  count(ECO, sort = TRUE)


library(stringr)

uniprot_mammalia_signal_confirmed_exp %>%
  mutate(ECO = str_extract_all(`Signal peptide`, "ECO:\\d+")) %>%
  mutate(ECO = map(ECO, ~ if (length(.x) == 0) NA_character_ else .x)) %>%
  unnest(ECO) %>%
  distinct(Entry, ECO) %>%
  count(ECO, sort = TRUE)

################################################################################


uniprot_mammalia_signal_confirmed_exp$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", uniprot_mammalia_signal_confirmed_exp$`Signal peptide`) 
uniprot_mammalia_signal_confirmed_exp$SPuniprot <- as.numeric(uniprot_mammalia_signal_confirmed_exp$SPuniprot)

mean(uniprot_mammalia_signal_confirmed_exp$SPuniprot, na.rm = TRUE)



nrow(uniprot_mammalia_signal_confirmed_exp)

### List for venn diagram of dataset for mammalia
SignalPeptides_mammalia_all <- list(
  "UniProt (n = 1463)" = uniprot_mammalia_signal_confirmed_exp$Entry,
  "Signal Peptide DB (n = 2109)" = signalpeptide_mammalia_confirmed$Accession.Number,
  "SignalP (n = 2061)" = signalP6_table_mammalia_an$ProteinID
)

##Vienn graph for intersection between signalpeptide.de DB, Uniprot and SignalP

mamm_sp_venn_all <- ggvenn(SignalPeptides_mammalia_all, 
       fill_color = c("#0073C2FF", "pink", "white"),
       set_name_size = 5,text_size = 5,show_percentage = FALSE)
mamm_sp_venn_all

ggsave(
  filename = "../../../Images/venn_signalpeptides_mammalia_all.png",
  plot = mamm_sp_venn_all,
  width = 7.32,
  height = 6.34,
  dpi = 300
)


length(uniprot_mammalia_signal_confirmed_exp$Entry) ## nb of proteins (mammalia, confirmed in UniProt)
length(signalpeptide_mammalia_confirmed$Accession.Number) ## nb of proteins (mammalia, confirmed in signalpeptide DB)
length(signalP6_table_mammalia_an$ProteinID) ## nb of proteins (mammalia, annotated by SignalP)

## Merge merged_sp_mammalia with confirmed presence of signal peptide exp (mammalia) in Uniprot
merged_sp_mammalia_uniprot <- merge(merged_sp_mammalia, uniprot_mammalia_signal_confirmed_exp, by.x = "ProteinID", by.y = "Entry") # intersection of 3 datasets
merged_sp_mammalia_uniprot$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", merged_sp_mammalia_uniprot$`Signal peptide`) ## extraction of end position from UniProt
merged_sp_mammalia_uniprot$SPuniprot <- as.numeric(merged_sp_mammalia_uniprot$SPuniprot)


length(merged_sp_mammalia_uniprot$Pos_end)
length(merged_sp_mammalia_uniprot$SPuniprot)

### Verification of the presence of NA
length(na.omit(merged_sp_mammalia_uniprot$Pos_end))
length(na.omit(merged_sp_mammalia_uniprot$SPuniprot))

length(uniprot_mammalia_signal_confirmed_exp$Entry)
length(merged_sp_mammalia$ProteinID)



## Graph Signal peptide lengths in Mammalia (signalpeptide.de vs UniProt)
range_xy <- range(
  merged_sp_mammalia_uniprot$Length.x,
  merged_sp_mammalia_uniprot$SPuniprot,
  na.rm = TRUE
)

ggplot(merged_sp_mammalia_uniprot, aes ( x = merged_sp_mammalia_uniprot$Length.x, y = merged_sp_mammalia_uniprot$SPuniprot)) +
  geom_point() +
  theme_classic() +
  xlab("Length in DB") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in Mammalia") +
  coord_equal(xlim = range_xy, ylim = range_xy) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 50
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )


## Graph Signal peptide lengths in Mammalia (SignalP vs UniProt)
signal_peptide_lengths_mammalia <- ggplot(merged_sp_mammalia_uniprot, aes ( x = merged_sp_mammalia_uniprot$Pos_end, y = merged_sp_mammalia_uniprot$SPuniprot)) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in Mammalia (n=1242, confirmed exp)") +
  geom_smooth(method = 'lm',
              se = TRUE, na.rm = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 65
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )

signal_peptide_lengths_mammalia

ggsave(
  filename = "../../../Images/signal_peptide_lengths_mammalia.png",
  plot = signal_peptide_lengths_mammalia,
  width = 7.48,
  height = 4.85,
  dpi = 300
)


length(merged_sp_mammalia_uniprot$Pos_end)
length(merged_sp_mammalia_uniprot$SPuniprot)


### Boxplot for Mammalia SignalP vs UniProt

png(
  filename = "../../../Images/boxplot_mammalia_signalp_uniprot.png",
  width = 6.8,
  height = 5.00,
  units = "in",
  res = 300
)

boxplot(list(Pos_end = merged_sp_mammalia_uniprot$Pos_end,
             SPuniprot = merged_sp_mammalia_uniprot$SPuniprot),
        ylab = "Length of SP",
        main = "Signal peptide lengths in Mammalia (SignalP vs UniProt, n = 1242)", ###intersection
        names = c("SignalP prediction", "UniProt"), col = c("skyblue", "pink"),
        boxwex = 0.4
)

dev.off()

length(merged_sp_mammalia_uniprot$SPuniprot)


## Statistics

###Shapiro-Wilk Test for Normality ####
### If the test is non-significant (p>. 05) it tells us that the distribution of 
### the sample is not significantly different from a normal distribution. 
### If, however, the test is significant (p < . 05) then the distribution in 
### question is significantly different from a normal distribution.



shapiro.test(na.omit(merged_sp_mammalia_uniprot$Pos_end)) ### p-value < 2.2e-16 not normal distribuion
shapiro.test(na.omit(merged_sp_mammalia_uniprot$SPuniprot)) ### p-value < 2.2e-16 not normal distribuion
hist(na.omit(merged_sp_mammalia_uniprot$Pos_end))
hist(na.omit(merged_sp_mammalia_uniprot$SPuniprot))
length(na.omit(merged_sp_mammalia_uniprot$Pos_end))  
length(na.omit(merged_sp_mammalia_uniprot$SPuniprot)) 

wilcox.test(merged_sp_mammalia_uniprot$Pos_end, merged_sp_mammalia_uniprot$SPuniprot, paired = TRUE) # statistically not different p-value = 0.3622


summary(merged_sp_mammalia_uniprot$Pos_end)
summary(merged_sp_mammalia_uniprot$SPuniprot)


# Graph length SP vs Organism (boxplot) Length - from signalpeptide.de mammalia
ggplot(data = merged_sp_mammalia, aes(x = Organism, y = Length, fill = Organism)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(title = "Mammalia SP length vs Organism") +
  xlab("Organism") +
  ylab("Length of peptide (db signalpeptide.de)")


### For do labels like a number of observations
nbOrganism <- merged_sp_mammalia_uniprot %>%
  add_count(Organism.x)

length(nbOrganism)

nbOrganism$n <- as.character(nbOrganism$n)
table(nbOrganism$n)
length(nbOrganism$n)
table(merged_sp_mammalia_uniprot$Organism.x) ### have a number of observations
organism_counts <- as.data.frame(table(merged_sp_mammalia_uniprot$Organism.x))

### to add a new label (column) which contains information about the organism and the number of frequency
organism_counts$label <- paste0(
  organism_counts$Var1,
  " (n=", organism_counts$Freq, ")" 
)

merged_sp_mammalia_uniprot2 <- merged_sp_mammalia_uniprot %>%
  left_join(organism_counts, by = c("Organism.x" = "Var1"))



# Graph length SP vs Organism (boxplot) Length - from SignalP mammalia
ggplot(data = merged_sp_mammalia_uniprot2, aes(x = label, y = Pos_end, fill = Organism.x)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(title = "Mammalia SP length vs Organism (n = 1242) ") +
  xlab("Organism") +
  ylab("Length of signal peptide (SignalP)") +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )

ggsave(
  filename = "../../../Images/sp_length_mammalia_uniprot_signalp_organism.png",
  width = 7.48,
  height = 4.85,
  dpi = 300)

length(merged_sp_mammalia_uniprot2$SPuniprot)


# Graph length SP vs Organism (boxplot) Length - from Uniprot mammalia
ggplot(data = merged_sp_mammalia_uniprot2, aes(x = label, y = SPuniprot, fill = Organism.x)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(title = "Mammalia SP length vs Organism (n = 1242) ") +
  xlab("Organism") +
  ylab("Length of signal peptide (UniProt)") +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )

ggsave(
  filename = "../../../Images/sp_length_mammalia_uniprot_organism.png",
  width = 7.48,
  height = 4.85,
  dpi = 300)

### FALSE and TRUE prediction of the SP length

comparison <- merged_sp_mammalia_uniprot2 %>%
  filter(!is.na(Pos_end), !is.na(SPuniprot)) %>%
  mutate(
    match = Pos_end == SPuniprot
  )

table(comparison$match) ### FALSE = 162 (13%), TRUE = 1080 (87%)


###### The difference between the predicted length by SignalP 6 and actual information in UniProt (mammalia)

merged_sp_mammalia_uniprot$diff = merged_sp_mammalia_uniprot$Pos_end/merged_sp_mammalia_uniprot$SPuniprot



png(
  filename = "../../../Images/hist_diff_mammalia.png",
  width = 604,
  height = 533
)

hist(merged_sp_mammalia_uniprot$diff, main = "The difference between SignalP and UniProt
     (ratio of length, mammalia)", xlab = "Ratio")

dev.off()


### Histograms for precision between 0.9 and 1.0
hist(
  merged_sp_mammalia_uniprot$diff[
    merged_sp_mammalia_uniprot$diff >= 0.9 &
      merged_sp_mammalia_uniprot$diff <= 1.0
  ],
  main = "The difference between SignalP and UniProt\n(ratio of length, mammalia)",
  xlab = "Ratio",
  breaks = seq(0.9, 1.0, by = 0.01)
)


png(
  filename = "../../../Images/hist_diff_mammalia_small.png",
  width = 604,
  height = 533
)

hist(
  merged_sp_mammalia_uniprot$diff[
    merged_sp_mammalia_uniprot$diff >= 0.99 &
      merged_sp_mammalia_uniprot$diff <= 1.0
  ],
  main = "The difference between SignalP and UniProt\n(ratio of length, mammalia)",
  xlab = "Ratio",
  breaks = seq(0.99, 1.0, by = 0.001)
)

dev.off()


################################################################################
############################## Bacteria ########################################
################################################################################

# Read the file from SignalP 6.0 results ("output.gff3")
signalP6_table_bacteria <- read.table("bacteria/output_bacteria_confirmed.gff3", sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE )
colnames(signalP6_table_bacteria) <- c("ProteinId", "App", "Type", "Pos_start", "Pos_end", "Score")
signalP6_table_bacteria_an <- signalP6_table_bacteria[0:6] ### 1118 confirmed
hist(signalP6_table_bacteria_an$Pos_end)


# Read the file from signalpeptide.de ("signalpeptide_bacteria.csv")
signalpeptide_bacteria <- read.table(
  "bacteria/signalpeptide_bacteria.csv",
  sep = ",",
  header = TRUE,
  quote = "\"",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

# Filter the table 'signalpeptide_bacteria' and choose the lines where 'SP.Status' == 'confirmed' (1161 proteins) 
signalpeptide_bacteria_confirmed <- signalpeptide_bacteria %>%
  filter(`SP.Status` == "confirmed")


## Merge output from signalpeptide.de (confirmed) and output from SignalP 6.0
merged_sp_bacteria <- merge(signalP6_table_bacteria_an, signalpeptide_bacteria_confirmed,
                            by.x = "ProteinId",
                            by.y = "Accession.Number",
                            all = TRUE )

merged_sp_bacteria$Length <- as.numeric(merged_sp_bacteria$Length)


## Table downloaded from UniProt with confirmed presence of signal peptide (bacteria) - bash command via query Uniprot / only with evidence ECO:0000269 !!!!!
uniprot_bacteria_signal_confirmed_exp <- read_tsv(
  "bacteria/uniprotkb_taxonomy_id_2_AND_ft_sign_exp_2026_05_11.tsv",
  show_col_types = FALSE)



uniprot_bacteria_signal_confirmed_exp$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", uniprot_bacteria_signal_confirmed_exp$`Signal peptide`) ## extraction of end position from UniProt
uniprot_bacteria_signal_confirmed_exp$SPuniprot <- as.numeric(uniprot_bacteria_signal_confirmed_exp$SPuniprot)


signalpeptide_bacteria_ids_clean <- signalpeptide_bacteria_confirmed$Accession.Number[
  !is.na(signalpeptide_bacteria_confirmed$Accession.Number) &
    !grepl("<<", signalpeptide_bacteria_confirmed$Accession.Number)
]

SignalPeptides_bacteria_db_sp_all <- list(
  "Signal Peptide DB " = signalpeptide_bacteria_ids_clean,
  "UniProt " = uniprot_bacteria_signal_confirmed_exp$Entry) ### presence of siignal peptide confirmed
len_db_sp_venn_bact_uniprot2 <- ggvenn(SignalPeptides_bacteria_db_sp_all, 
                                      fill_color = c("#0073C2FF", "pink"),
                                      set_name_size = 5,text_size = 5,show_percentage = FALSE)

len_db_sp_venn_bact_uniprot2

ggsave(
  filename = "../../../Images/len_db_sp_venn_bact_uniprot2.png",
  plot = len_db_sp_venn_bact_uniprot2,
  width = 7.55,
  height = 6.79,
  dpi = 300
)








##### To calculate ECO in UniProt data bacteria ################################
uniprot_bacteria_signal_confirmed_exp %>%
  mutate(ECO = str_extract_all(`Signal peptide`, "ECO:\\d+")) %>%
  mutate(ECO = map(ECO, ~ if (length(.x) == 0) NA_character_ else .x)) %>%
  unnest(ECO) %>%
  distinct(Entry, ECO) %>%
  count(ECO, sort = TRUE)
################################################################################



uniprot_bacteria_signal_confirmed_exp$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", uniprot_bacteria_signal_confirmed_exp$`Signal peptide`) 
uniprot_bacteria_signal_confirmed_exp$SPuniprot <- as.numeric(uniprot_bacteria_signal_confirmed_exp$SPuniprot)
mean(uniprot_bacteria_signal_confirmed_exp$SPuniprot, na.rm = TRUE)


nrow(uniprot_bacteria_signal_confirmed_exp)
length(na.omit(uniprot_bacteria_signal_confirmed_exp$Entry))
length(na.omit(signalpeptide_bacteria_confirmed$Accession.Number))


### Venn graph for all datasets (bacteria)
SignalPeptides_bacteria_db_sp_all <- list(
  "SignalP (n = 1118)" = signalP6_table_bacteria_an$ProteinId,
  "Signal Peptide DB (n = 1161)" = signalpeptide_bacteria_confirmed$Accession.Number,
  "UniProt (n = 824)" = uniprot_bacteria_signal_confirmed_exp$Entry)
len_db_sp_venn_all <- ggvenn(SignalPeptides_bacteria_db_sp_all, 
                             fill_color = c("#0073C2FF", "pink", "white"),
                             set_name_size = 5,text_size = 5,show_percentage = FALSE)

len_db_sp_venn_all

ggsave(
  filename = "../../../Images/venn_signalpeptides_bac_all.png",
  plot = len_db_sp_venn_all,
  width = 7.55,
  height = 6.79,
  dpi = 300
)


## Graph Signal peptide lengths in Bacteria (signalpeptide.de vs SignalP6)
ggplot(merged_sp_bacteria, aes ( x = merged_sp_bacteria$Pos_end, y = merged_sp_bacteria$Length)) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP6") +
  ylab("Length in DB") +
  labs(title = "Signal peptide lengths in Bacteria") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 65
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )

ggsave(
  filename = "../../../Images/sp_length_bacteria_signalp_signal_db.png",
  plot = last_plot(),
  width = 8,
  height = 5,
  dpi = 300
)




## Merge merged_sp_bacteria with confirmed presence of signal peptide (bacteria)
merged_sp_bacteria_uniprot <- merge(merged_sp_bacteria, uniprot_bacteria_signal_confirmed_exp, by.x = "ProteinId", by.y = "Entry") # intersection 
merged_sp_bacteria_uniprot$SPuniprot<- sub(".*\\.\\.(\\d+);.*", "\\1", merged_sp_bacteria_uniprot$`Signal peptide`) ## extraction of end position from UniProt
merged_sp_bacteria_uniprot$SPuniprot <- as.numeric(merged_sp_bacteria_uniprot$SPuniprot)


length(merged_sp_bacteria_uniprot$Pos_end)
length(merged_sp_bacteria_uniprot$SPuniprot)
length(na.omit(merged_sp_bacteria_uniprot$Pos_end))
length(na.omit(merged_sp_bacteria_uniprot$SPuniprot))


length(na.omit(uniprot_bacteria_signal_confirmed_exp$Entry))
length(na.omit(merged_sp_bacteria$ProteinId))

length(uniprot_mammalia_signal_confirmed_exp$Entry)
length(merged_sp_mammalia$ProteinID)


## Graph Signal peptide lengths in bacteria (signalpeptide.de vs UniProt)
range_xy <- range(
  merged_sp_bacteria_uniprot$Length.x,
  merged_sp_bacteria_uniprot$SPuniprot,
  na.rm = TRUE
)

ggplot(merged_sp_bacteria_uniprot, aes ( x = merged_sp_bacteria_uniprot$Length.x, y = merged_sp_bacteria_uniprot$SPuniprot)) +
  geom_point() +
  theme_classic() +
  xlab("Length in DB") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in bacteria") +
  coord_equal(xlim = range_xy, ylim = range_xy) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE)


## Graph Signal peptide lengths in bacteria (SignalP vs UniProt)
ggplot(merged_sp_bacteria_uniprot, aes ( x = merged_sp_bacteria_uniprot$Pos_end, y = merged_sp_bacteria_uniprot$SPuniprot)) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in bacteria") +
  geom_smooth(method = 'lm',
              se = TRUE, na.rm = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  stat_cor(
    method = "pearson",
    label.x = 5,
    label.y = 65
  ) +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  ) +
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )


ggsave(
  filename = "../../../Images/sp_length_bacteria_signalp_uniprot.png",
  plot = last_plot(),
  width = 8,
  height = 5,
  dpi = 300
)


length(na.omit(merged_sp_bacteria_uniprot$Pos_end)) ### 629
length(na.omit(merged_sp_bacteria_uniprot$SPuniprot)) ### 643


### FALSE and TRUE prediction of the SP length

comparison_bact <- merged_sp_bacteria_uniprot %>%
  filter(!is.na(Pos_end), !is.na(SPuniprot)) %>%
  mutate(
    match = Pos_end == SPuniprot
  )

table(comparison_bact$match) ### FALSE = 56 (9%), TRUE = 573 (91%)



### Boxplot for bacteria SignalP vs UniProt

png(
  filename = "../../../Images/boxplot_bacteria_signalp_uniprot.png",
  width = 6.40,
  height = 6.53,
  units = "in",
  res = 300
)


ymax <- max(c(merged_sp_bacteria_uniprot$Pos_end, merged_sp_bacteria_uniprot$SPuniprot), na.rm = TRUE) * 1.25

boxplot(list(Pos_end = (na.omit(merged_sp_bacteria_uniprot$Pos_end)),
             SPuniprot = (na.omit(merged_sp_bacteria_uniprot$SPuniprot))),
        ylab = "Length of SP",
        main = "Signal peptide lengths in Bacteria (SignalP vs UniProt, n = 629)", ###intersection
        names = c("SignalP prediction", "UniProt"), col = c("skyblue", "pink"),
        boxwex = 0.4,
        ylim = c(0, ymax)
)


y <- max(c(merged_sp_bacteria_uniprot$Pos_end, merged_sp_bacteria_uniprot$SPuniprot), na.rm = TRUE) * 1.08

segments(1, y, 2, y)
segments(1, y, 1, y - 3)
segments(2, y, 2, y - 3)

text(
  x = 1.5,
  y = y + 3.5,
  labels = "***",
  cex = 1.5
)

dev.off()


### Statistics ################################################################


shapiro.test(na.omit(merged_sp_bacteria_uniprot$Pos_end)) ### p-value < 2.2e-16 not normal distribuion
shapiro.test(na.omit(merged_sp_bacteria_uniprot$SPuniprot)) ### p-value < 2.2e-16 not normal distribuion
hist(na.omit(merged_sp_bacteria_uniprot$Pos_end))
hist(na.omit(merged_sp_bacteria_uniprot$SPuniprot))
length(na.omit(merged_sp_bacteria_uniprot$Pos_end))  
length(na.omit(merged_sp_bacteria_uniprot$SPuniprot)) 

mean(na.omit(merged_sp_bacteria_uniprot$Pos_end)) 
median(na.omit(merged_sp_bacteria_uniprot$Pos_end)) 
mean(na.omit(merged_sp_bacteria_uniprot$SPuniprot)) 
median(na.omit(merged_sp_bacteria_uniprot$SPuniprot)) 




wilcox.test(merged_sp_bacteria_uniprot$Pos_end, merged_sp_bacteria_uniprot$SPuniprot, paired = TRUE) # statistically different p-value = 3.333e-07


summary(merged_sp_mammalia_uniprot$Pos_end)
summary(merged_sp_mammalia_uniprot$SPuniprot)

###### The difference between the predicted length by SignalP 6 and actual information in UniProt

merged_sp_bacteria_uniprot$diff = merged_sp_bacteria_uniprot$Pos_end/merged_sp_bacteria_uniprot$SPuniprot



png(
  filename = "../../../Images/hist_diff_bacteria.png",
  width = 604,
  height = 533
)

hist(merged_sp_bacteria_uniprot$diff, main = "The difference between SignalP and UniProt
     (ratio of length, bacteria)", xlab = "Ratio")

dev.off()

### Histograms for precision between 0.9 and 1.0
hist(
  merged_sp_bacteria_uniprot$diff[
    merged_sp_bacteria_uniprot$diff >= 0.9 &
      merged_sp_bacteria_uniprot$diff <= 1.0
  ],
  main = "The difference between SignalP and UniProt\n(ratio of length, bacteria)",
  xlab = "Ratio",
  breaks = seq(0.9, 1.0, by = 0.01)
)



png(
  filename = "../../../Images/hist_diff_bacteria_small.png",
  width = 604,
  height = 533
)

hist(
  merged_sp_bacteria_uniprot$diff[
    merged_sp_bacteria_uniprot$diff >= 0.99 &
      merged_sp_bacteria_uniprot$diff <= 1.0
  ],
  main = "The difference between SignalP and UniProt\n(ratio of length, bacteria)",
  xlab = "Ratio",
  breaks = seq(0.99, 1.0, by = 0.001)
)

dev.off()

################# BOXPLOTS UNIPROT #############################################

library(tidyverse)

df_box <- bind_rows(
  uniprot_viruses_sppr  %>% transmute(Group = "Viruses",  SPuniprot = SPuniprot),
  uniprot_mammalia_signal_confirmed_exp %>% transmute(Group = "Mammalia", SPuniprot = SPuniprot),
  uniprot_bacteria_signal_confirmed_exp %>% transmute(Group = "Bacteria", SPuniprot = SPuniprot)
)

ggplot(df_box, aes(x = Group, y = SPuniprot, fill = Group)) +
  geom_boxplot() +
  labs(
    title = "Signal peptide length by group (UniProt Database)",
    x = "",
    y = "Signal peptide length"
  ) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 150, by = 10)) +
  coord_flip() +
  theme(legend.position = "none")
  

ggsave(
  filename = "../../../Images/uniprot_all_group_sp_length.png",
  plot = last_plot(),
  width = 6.63,
  height = 4.66,
  dpi = 300
)

































