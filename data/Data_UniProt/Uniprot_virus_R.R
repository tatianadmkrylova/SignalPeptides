getwd()
setwd("/home/tatiana/Documents/Data_UniProt")
getwd()

# Necessary libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggvenn)
install.packages("ggpubr", dependenc)
library(ggpubr)
library(readr)


############## Viruses UniProt ################################################

# Read the file from SignalP 6.0 results ("output.gff3")
signalP6_table_viruses_uniprot_sppr <- read.table(
  "Res_SignalP_uniProt_viruses/output_signalp_viruses_uniprot.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

# Select columns from signalP6_table
signalP6_table_viruses_uniprot_sppr_an <- signalP6_table_viruses_uniprot_sppr[0:6]


# Add the headers 
colnames(signalP6_table_viruses_uniprot_sppr_an) <- c("ProteinId", "App", "Type", "Pos_start", "Pos_end", "Score")


# Table downloaded from UniProt with confirmed presence of signal peptide (viruses) - bash command via query Uniprot
uniprot_viruses_sppr <- read.table("uniprotkb_taxonomy_viruses_AND_ft_sign.tsv", sep = "\t",
                                               header = TRUE,
                                               quote = "",
                                               fill = TRUE,
                                               comment.char = "",
                                               stringsAsFactors = FALSE)
uniprot_viruses_sppr$SPuniprot<- sub(".*\\.\\.(\\d+).*",
                                                 "\\1", uniprot_viruses_sppr$`Signal.peptide`) ## extraction of end position from UniProt
uniprot_viruses_sppr$SPuniprot <- as.numeric(uniprot_viruses_sppr$SPuniprot)


### Venn graph for all datasets (viruses) CHANGER

signalPeptides_viruses_uniprot <- list(
  "SignalP (n = 1202)" = signalP6_table_viruses_uniprot_sppr_an$ProteinId,
  "UniProt (n = 1579)" = uniprot_viruses_sppr$Accession.Number) ### presence of siignal peptide confirmed
venn_signalPeptides_viruses_uniprot <- ggvenn(signalPeptides_viruses_uniprot, 
                                 fill_color = c("#0073C2FF", "pink"),
                                 set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalPeptides_viruses_uniprot

## Merge output from UniProt and output from SignalP 6.0 
spP6_uniprot_merged_viruses <- merge(signalP6_table_viruses_uniprot_sppr_an, uniprot_viruses_sppr, by.x = "ProteinId",
                            by.y = "Accession.Number")


no_signalp_detection_uniprot <-  spP6_uniprot_merged_viruses %>% ### signalp6 didn't recognize
  filter(is.na(spP6_uniprot_merged_viruses$App))


#### Graph Length in SignalP vs Length in UniProt

ggplot(
  spP6_uniprot_merged_viruses,
  aes(x = Pos_end, y = SPuniprot, na.rm = TRUE)
) +
  geom_point() +
  theme_classic() +
  xlab("Length in SignalP") +
  ylab("Length in UniProt") +
  labs(title = "Signal peptide lengths in Viruses / UniProt data") +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  scale_x_continuous(limits = c(0, 70), breaks = seq (0, 70, by = 5)) +
  scale_y_continuous(limits = c(0, 70), breaks = seq (0, 70, by = 5)) 


#### Boxplots of the distribution for non detected proteins
boxplot( data = no_signalp_detection_uniprot, x = no_signalp_detection_uniprot$SPuniprot, xlab = "Not detected by SignalP", ylab = "Length of signal peptide", main = "UniProt Viruses")
boxplot( data = no_signalp_detection_uniprot, x = no_signalp_detection_uniprot$Length, xlab = "Not detected by SignalP", ylab = "Length of protein", main = "UniProt Viruses")

#### Boxplots of the distribution for detected proteins
boxplot( data = spP6_uniprot_merged_viruses_intersection, x = spP6_uniprot_merged_viruses_intersection$SPuniprot, xlab = "Detected by SignalP", ylab = "Length of signal peptide", main = "UniProt Viruses")
boxplot( data = spP6_uniprot_merged_viruses_intersection, x = spP6_uniprot_merged_viruses_intersection$Length, xlab = "Detected by SignalP", ylab = "Length of protein", main = "UniProt Viruses")




## Statistics

###Shapiro-Wilk Test for Normality ####
### If the test is non-significant (p>. 05) it tells us that the distribution of 
### the sample is not significantly different from a normal distribution. 
### If, however, the test is significant (p < . 05) then the distribution in 
### question is significantly different from a normal distribution.



shapiro.test(na.omit(no_signalp_detection_uniprot$SPuniprot)) ### p-value < 2.2e-16 not normal distribuion
shapiro.test(na.omit(no_signalp_detection_uniprot$Length)) ### p-value < 2.2e-16 not normal distribuion
hist(na.omit(no_signalp_detection_uniprot$SPuniprot))
hist(na.omit(no_signalp_detection_uniprot$Length))
length(na.omit(no_signalp_detection_uniprot$SPuniprot)) #### 372
length(na.omit(no_signalp_detection_uniprot$Length)) ### 377

wilcox.test(no_signalp_detection_uniprot$Length, spP6_uniprot_merged_viruses_intersection$Length, paired = FALSE) # statistically different p-value = 1.671e-11
wilcox.test(no_signalp_detection_uniprot$SPuniprot, spP6_uniprot_merged_viruses_intersection$SPuniprot, paired = FALSE) # statistically different p-value < 2.2e-16

summary(merged_spP6_uniprot_vir_intersection$Pos_end) ### Mean 21.34
summary(merged_spP6_uniprot_vir_intersection$SPuniprot) ### Mean 21.07



min(no_signalp_detection_uniprot$SPuniprot, na.rm = TRUE)
max(no_signalp_detection_uniprot$SPuniprot, na.rm = TRUE)
mean(no_signalp_detection_uniprot$SPuniprot, na.rm = TRUE)
median(no_signalp_detection_uniprot$SPuniprot, na.rm = TRUE)

min(no_signalp_detection_uniprot$Length, na.rm = TRUE)
max(no_signalp_detection_uniprot$Length, na.rm = TRUE)
mean(no_signalp_detection_uniprot$Length, na.rm = TRUE)
median(no_signalp_detection_uniprot$Length, na.rm = TRUE)


min( spP6_uniprot_merged_viruses_intersection$SPuniprot, na.rm = TRUE)
max( spP6_uniprot_merged_viruses_intersection$SPuniprot, na.rm = TRUE)
mean( spP6_uniprot_merged_viruses_intersection$SPuniprot, na.rm = TRUE)
median( spP6_uniprot_merged_viruses_intersection$SPuniprot, na.rm = TRUE)

min( spP6_uniprot_merged_viruses_intersection$Length, na.rm = TRUE)
max( spP6_uniprot_merged_viruses_intersection$Length, na.rm = TRUE)
mean( spP6_uniprot_merged_viruses_intersection$Length, na.rm = TRUE)
median( spP6_uniprot_merged_viruses_intersection$Length, na.rm = TRUE)

summary(spP6_uniprot_merged_viruses_intersection$SPuniprot)
summary(spP6_uniprot_merged_viruses_intersection$Length)

# Graph length SP vs Organism (boxplot) Length - from SignalP mammalia
ggplot(data = no_signalp_detection_uniprot, aes(x = no_signalp_detection_uniprot$Gene.Names, y = no_signalp_detection_uniprot$SPuniprot, fill = Gene.Names)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(title = "Mammalia SP length vs Organism / UniProt data ") +
  xlab("Organism") +
  ylab("Length of signal peptide (SignalP)") +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )


library(stringr)
sum(str_count(no_signalp_detection_uniprot$Gene.Names), na.rm = TRUE)
table(no_signalp_detection_uniprot$Gene.Names, na.rm = TRUE)


res_nb_genes <- no_signalp_detection_uniprot %>%
  group_by(Gene.Names) %>%
  summarise(n=n())

res_nb_proteins <- no_signalp_detection_uniprot %>%
  group_by(Protein.names) %>%
  summarise(n=n()) 


res_nb_genes_filt_1 <- res_nb_genes %>%
  filter(n > 1)

res_nb_proteins_filt_1 <- res_nb_proteins %>%
  filter(n > 1)

no_signalp_detection_uniprot_names <- no_signalp_detection_uniprot %>%
  select(Gene.Names, Protein.names, ProteinId, Organism, SPuniprot)

no_signalp_detection_uniprot_names$Organism_short <- trimws(sub("\\s*\\(.*$", "", no_signalp_detection_uniprot_names$Organism))


# Graph length SP vs Organism (boxplot) Length - from Uniprot viruses
ggplot(data = no_signalp_detection_uniprot_names, aes(x = Organism_short, y = SPuniprot, fill = Organism_short)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(title = "Virus SP length vs Organism (UniProt), no SignalP ") +
  xlab("Organism") +
  ylab("Length of signal peptide (UniProt)") +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )


res_nb_organisms_short <- no_signalp_detection_uniprot_names %>%
  group_by(Organism_short) %>%
  summarise(n=n())

res_nb_organisms_short
sum(res_nb_organisms_short$n)

res_nb_organisms_filt_1 <- res_nb_organisms_short %>%
  filter(n > 1)



ggplot(data = no_signalp_detection_uniprot_names, aes(x = Organism_short, y = SPuniprot, fill = Organism_short)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(title = "Virus SP length vs Organism (UniProt), no SignalP ") +
  xlab("Organism") +
  ylab("Length of signal peptide (UniProt)") +
  scale_y_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, by = 5)
  )




## Merge output from UniProt and output from SignalP 6.0 #intersection
spP6_uniprot_merged_viruses_intersection <- merge(signalP6_table_viruses_uniprot_sppr_an, uniprot_viruses_sppr, by.x = "ProteinId",
                                     by.y = "Accession.Number")



spP6_uniprot_merged_viruses_intersection$Organism_short <- trimws(sub("\\s*\\(.*$", "",spP6_uniprot_merged_viruses_intersection$Organism))


### Venn graph for all datasets (viruses) for no and yes SignalP prediction

signalPeptides_viruses_uniprot_yes_no <- list(
  "SignalP (n = 1202)" = spP6_uniprot_merged_viruses_intersection$Protein.names,
  "Without SignalP (n = 377)" = no_signalp_detection_uniprot_names$Protein.names) ### presence of siignal peptide confirmed
venn_signalPeptides_viruses_uniprot_yes_no <- ggvenn(signalPeptides_viruses_uniprot_yes_no, 
                                              fill_color = c("#0073C2FF", "pink"),
                                              set_name_size = 5,text_size = 5,show_percentage = FALSE)

venn_signalPeptides_viruses_uniprot_yes_no


## Merge with and without SignalP 6.0 #intersection
intersection_with_without_signalp_pr_names <- merge(spP6_uniprot_merged_viruses_intersection, no_signalp_detection_uniprot_names, by = "Protein.names")


common_names <- intersect(
  unique(spP6_uniprot_merged_viruses_intersection$Protein.names), unique(no_signalp_detection_uniprot_names$Protein.names))
length(common_names)

sp6_intersection_28 <- spP6_uniprot_merged_viruses_intersection[spP6_uniprot_merged_viruses_intersection$Protein.names %in% common_names, ]



################################################################################
################################################################################
################################################################################
################################################################################

uniprot_full_seq_prot_viruses_no_doublons <- read_csv("../scripts/uniprot_full_seq_prot_viruses_no_doublons.csv")
View(uniprot_full_seq_prot_viruses_no_doublons)



## Merge output from UniProt and output from SignalP 6.0 and UniProt without doublons - intersection
 merge_uniprot_without_doublons_sp_uni<- merge(spP6_uniprot_merged_viruses, uniprot_full_seq_prot_viruses_no_doublons, by = "ProteinId")




 ### Venn graph for all datasets (viruses) CHANGER
 
 signalPeptides_viruses_uniprot_2 <- list(
   "SignalP (n = 1202)" = signalP6_table_viruses_uniprot_sppr_an$ProteinId,
   "UniProt (n = 1579)" = uniprot_viruses_sppr$Accession.Number, ### presence of siignal peptide confirmed
   "UniProt without doublons (n = 1505)" = uniprot_full_seq_prot_viruses_no_doublons$ProteinId)
 venn_signalPeptides_viruses_uniprot2 <- ggvenn(signalPeptides_viruses_uniprot_2, 
                                               fill_color = c("#0073C2FF", "pink", "white"),
                                               set_name_size = 5,text_size = 5,show_percentage = FALSE)
 
 venn_signalPeptides_viruses_uniprot2



 ## Merge output from UniProt and output from SignalP 6.0 - data from uniprot
 spP6_uniprot_merged_viruses2 <- merge(signalP6_table_viruses_uniprot_sppr_an, uniprot_viruses_sppr, by.x = "ProteinId",
                                      by.y = "Accession.Number", all.y = TRUE)

 
 
 
 
 ## Merge output from UniProt and output from SignalP 6.0 and UniProt without doublons
 merge_uniprot_without_doublons_sp_uni2<- merge(spP6_uniprot_merged_viruses2, uniprot_full_seq_prot_viruses_no_doublons, by = "ProteinId", all.y = TRUE)
 
no_signalp_annotation_without_doublons <- merge_uniprot_without_doublons_sp_uni2 %>%
   filter(is.na(merge_uniprot_without_doublons_sp_uni2$App))
 
no_signalp_annotation_without_doublons_table <- no_signalp_annotation_without_doublons %>%
  select(ProteinId, Sequence, SPuniprot)

write.csv(no_signalp_annotation_without_doublons_table, "no_signalp_annotation_without_doublons_table.csv", row.names = FALSE) 
 
 
 
 
table(no_signalp_annotation_without_doublons$Virus.hosts)


proteinId_hosts_virus <- no_signalp_annotation_without_doublons %>%
  select(ProteinId, Virus.hosts)

library(stringr)
library(tidyr)
proteinId_hosts_virus_taxid <- proteinId_hosts_virus %>%
  mutate(TaxID_host = str_extract_all(Virus.hosts, "(?<=TaxID: )\\d+")) %>%
  unnest(TaxID_host)
 
proteinId_hosts_virus_taxid$TaxID_host
 
 
proteinId_hosts_virus_taxid_for_search <- proteinId_hosts_virus_taxid %>%
  select(TaxID_host)
 
proteinId_hosts_virus_taxid_for_search$TaxID_host <- as.numeric(proteinId_hosts_virus_taxid_for_search$TaxID_host) 
write.table(proteinId_hosts_virus_taxid_for_search, "proteinId_hosts_virus_taxid_for_search", row.names = FALSE, col.names = FALSE) 


taxid_host_viruses_merge <- merge(virus_nodublons_class_taxid_host, proteinId_hosts_virus_taxid , by.x = "TaxID", by.y = "TaxID_host")

table(taxid_host_viruses_merge$Class)

