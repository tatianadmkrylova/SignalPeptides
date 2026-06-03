getwd()
setwd("/Users/tatianakrylova/Documents/STAGE_LIRMM/Travail")
getwd()
gff <- read.table(
  "output_signalp6_all_sp.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)
colnames(gff) <- c(
  "ProteinId", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

merged <- merge(
  gff,
  signalpeptide_viruses,
  by.x = "ProteinId",
  by.y = "Accession Number",
  all.x = TRUE
)

merged2 <- merge(
  table_merged_for_analysis,
  hosts[, c("ProteinId","Host")],
  by.x = "ProteinId",
  by.y = "ProteinId",
  all.x = TRUE
)

merged2$humanHost <- ifelse(
  !is.na(merged2$Host) & grepl("Homo sapiens", merged2$Host, ignore.case = TRUE),
  "yes",
  "no"
) ### for verification if human present like a host

table(merged2$humanHost)

library(ggplot2)

model <- lm(merged2$Pos_end ~ merged2$Length, data = merged_copy)
merged_copy <- merged2
merged_copy$predicted <- predict(model)
merged_copy$residuals <- merged_copy$y - merged_copy$predicted

rmse <- sqrt(mean(merged_copy$residuals^2, na.rm = TRUE))

ggplot(merged_copy, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = paste("Linear regression, RMSE =", round(rmse, 3)),
    x = "x",
    y = "y"
  )


library("ggpubr")
df3 <- ggscatter(merged2, x = "Pos_end", y = "Length", 
                 add = "reg.line",add.params = list(color = "blue", fill = "lightgray"), 
                 conf.int = TRUE) +
    ylim(0,60)
df3
df3 + stat_cor(method = "pearson", label.x = 0, label.y = 50)

hist(merged2$Pos_end)
hist(merged2$Length)
x <- na.omit(merged2$Pos_end)
shapiro.test(x)
length(na.omit(merged2$Pos_end))

colours <- c("pink","pink")
ggplot(data = merged2) +
  geom_bar(aes(x = humanHost, fill = App)) +
  scale_fill_manual(values = colours) +
  theme_classic()
 

merged2 %>%
  filter(humanHost == "yes") %>%
  ggplot(aes(x = humanHost, fill = `SP Status`)) +
  geom_bar() +
  theme_classic()

human_yes <- merged2[merged2$humanHost == "yes", ]

ggplot(data = human_yes, aes(x = humanHost, fill = `SP Status`)) +
  geom_bar() +
  theme_classic() +
  geom_text(
    aes(label = after_stat(count), group = `SP Status`),
    stat = "count",
    position = position_stack(vjust = 0.5)
  )

ggplot(data = merged2) +
  geom_bar(aes(x = Host)) +
  coord_flip()

hasHost <- merged2[!is.na(merged2$Host), ]
str(hasHost)

library(tidyverse)

hasHostTable <- hasHost %>%
  select(ProteinId, App, Length, 'SP Status', Host)


hasHostTable_long <- hasHostTable %>%
  separate_rows(Host, sep = ",") %>%
  mutate(Host = trimws(Host))

ggplot(data = hasHostTable_long) +
  geom_bar(aes(x = Host, fill = App)) +
  coord_flip()

nbhasHost_filtered <- hasHostTable_long %>%
  add_count(Host) %>%
  filter(n >= 10)

ggplot(data = nbhasHost_filtered, aes(x = fct_infreq(Host), fill = App)) +
  geom_bar() +
  geom_text(
    aes(label = after_stat(count), group = App),
    stat = "count",
    size = 3,
    position = position_stack(vjust = 0.5)
  ) +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme_classic()

write.csv(hasHostTable, "table_for_looking_taxons.csv", row.names = FALSE)

hosts_taxID_copy_correct <- separate_rows(hosts_taxID_copy, Host, TaxonID, sep = ",\\s*", convert = TRUE)
write.csv(hosts_taxID_copy_correct, "table_hosts_taxons.csv", row.names = FALSE)

colnames(taxonomy) <- c("TaxonId","TaxonId", "type1", "type2", "taxonomy1", "taxonomy2")

taxonomy_split <- strsplit(taxonomy$taxonomy1, ",")
taxonomy_split[[1]]
library(dplyr)

library(tidyverse)
tax_long <- taxonomy %>%
  separate_longer_delim(taxonomy1, delim = ",") %>%
  separate(taxonomy1, into = c("rank", "name"), sep = "_", extra = "merge")


tax_wide <- tax_long %>%
  pivot_wider(
    id_cols = c(TaxonId, TaxonId.1),
    names_from = rank,
    values_from = name,
    values_fn = dplyr::first
  )


merged_taxID_class <- merge(
  table_hosts_taxons,
  tax_wide,
  by.x = "TaxonID",
  by.y = "TaxonId",
  all.x = TRUE
)



ggplot(data = merged_taxID_class, aes(x = fct_infreq(class))) +
  geom_bar(fill = "steelblue") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    size = 3,
    position = position_stack(vjust = 0.5)
  ) +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme_classic()

ggplot(data = merged_taxID_class, aes(x = fct_infreq(domain))) +
  geom_bar(fill = "steelblue") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    size = 3,
    position = position_stack(vjust = 0.5)
  ) +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme_classic()

library(ggplot2)

# graph length vs class (boxplot)
ggplot(merged_taxID_class, aes(x = fct_infreq(class), y = Length, fill = class)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme(legend.position = "none")

# graph length vs domain (boxplot)
ggplot(merged_taxID_class, aes(x = fct_infreq(domain), y = Length, fill = domain)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme(legend.position = "none")

# graph length vs App (boxplot)
ggplot(table_merged_for_analysis, aes(x = fct_infreq(App), y = Length, fill = App)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme(legend.position = "none")

# graph length vs domain (boxplot) with classes 
ggplot(merged_taxID_class, aes(x = fct_infreq(domain), y = Length, fill = class)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27)



library(dplyr)



merged_comparaison <- table_merged_for_analysis %>%
  left_join(
    uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11,
    by = c("ProteinId" = "Entry")
  ) %>%
  mutate(
    status = ifelse(
      ProteinId %in% uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$Entry,
      "in both",
      "not in Uniprot"
    )
  )

table(merged_comparaison$status)


merged_comparaison %>%
  filter(status == "in both") %>%
  count(`SP Status`)


merged_comparaison %>%
  filter(status == "in both") %>%
  count(`SP Status`) %>%
  ggplot(aes(x = `SP Status`, y = n, fill = n)) +
  geom_col(fill = c("indianred3", "lightskyblue", "lightpink", "lightsalmon"),
           color = "black", width = 0.5 ) + #color = border of bar
  labs(
    x = "SP Status",
    y = "Number of proteins"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text(
    aes(label = n),
    vjust = -0.3
  )


library(ggvenn)
library(ggplot2)

SignalPeptides <- list(
  "UniProt (N = 1579)" = uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$Entry,
  "Signal Peptide DB (N = 6197)" = table_merged_for_analysis$ProteinId
)

##Vienn graph for intersection between signalpeptide.de DB and Uniprot

ggvenn(SignalPeptides, 
       fill_color = c("#0073C2FF", "pink"),
       set_name_size = 5)

merged_comparaison2 <- merged_comparaison %>%
  filter( `SP Status`== "confirmed")

table(merged_comparaison2$status)

merged_comparaison2 %>%
  count(`status`) %>%
  ggplot(aes(x = `status`, y = n, fill = n)) +
  geom_col(fill = c("lightskyblue", "pink"), color = "black", width = 0.3) + #color = border of bar
  labs(
    x = "Status",
    y = "Number of proteins (confirmed)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text(
    aes(label = n),
    vjust = -0.3
  ) 

sum(is.na(merged_comparaison2$App)) # nb of proteins which SignalP didn't recognize
sum(!is.na(merged_comparaison2$App)) # nb of proteins which SignalP recognized


merged_comparaison2 %>% # count nb of proteins which presents in Uniprot and recognized by SignalP
  filter(status == "in both") %>%
  count(App)

merged_comparaison2 %>%  # count nb of proteins which presents in Uniprot
  filter(status == "in both") %>%
  count(status)

merged_comparaison2 %>%
  filter(status == "not in Uniprot") %>%
  count()
## graph for 215 signal peptides (confirmed in signalpeptide.de) + data from UniProt
ggplot(aes( x = `Length.x`,  y = `Length.y`), data = merged_comparaison2) + 
  geom_point(colour = "black",na.rm = TRUE) +
  theme_classic() +
  xlab("Length of signal peptide") +
  ylab("Length of peptide") +
  scale_x_continuous(n.breaks = 25) +
  scale_y_continuous(n.breaks = 20)


test <- cor.test(merged_comparaison2$Length.x, merged_comparaison2$Length.y)
test

sum(is.na(merged_comparaison2$Length.y))
typeof(merged_comparaison2$Length.y)


merged_taxID_organism <- merge(
  organism_taxID,
  taxonomy_unique_output_clean,
  by.x = "TaxonID_org",
  by.y = "TaxID",
  all.x = TRUE
)

class_labels <- merged_taxID_organism %>%
  count(Class) %>%
  mutate(label = paste0(Class, " (n=", n, ")"))


merged_taxID_organism2 <- merged_taxID_organism %>%
  left_join(class_labels, by = "Class")

# graph length SP vs class_organism (boxplot) Length - from signalpeptide.de
ggplot(data = merged_taxID_organism2, aes(x = label, y = Length, fill = Class)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme(legend.position = "none") +
  labs(title = "Virus SP length vs class of organism") +
  xlab("Class") +
  ylab("Length of peptide (db signalpeptide.de)")

table(merged_taxID_organism$Class)

## Extract the SP length from UniProt

uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$SPuniprot <- sub(".*\\.\\.(\\d+).*", "\\1", uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$`Signal peptide`)
table(uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$SPuniprot)

uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$SPuniprot <- as.numeric(uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$SPuniprot)

merged_taxID_organism_uniprot <- merge(
  merged_taxID_organism,
  uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11,
  by.x = "ProteinId",
  by.y = "Entry",
  all.x = TRUE
)

head(merged_taxID_organism$ProteinId)
head(uniprotkb_taxonomy_id_10239_AND_ft_sign_2026_05_11$Entry)

merged_taxID_organism_uniprot2 <- merged_taxID_organism_uniprot %>%
  left_join(class_labels, by = "Class")

# graph length SP vs class_organism (boxplot) Length - from uniprot
ggplot(data = merged_taxID_organism_uniprot2, aes(x = label, y = SPuniprot, fill = Class)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme(legend.position = "none") +
  labs(title = "Virus SP length uniprot vs class of organism (n = 6174)") +
  xlab("Class") +
  ylab("Length of SP from UniProt")

length(merged_taxID_organism_uniprot2$SPuniprot)





merged_taxID_organism_uniprot_signalP <- merge(
  merged_taxID_organism_uniprot,
  table_merged_for_analysis,
  by.x = "ProteinId",
  by.y = "ProteinId",
  all.x = TRUE
)

merged_taxID_organism_uniprot_signalP2 <- merged_taxID_organism_uniprot_signalP %>%
  left_join(class_labels, by = "Class")


# graph length SP vs class_organism (boxplot) Length - from SignalP6
ggplot(merged_taxID_organism_uniprot_signalP2, aes(x = label, y = Pos_end, fill = Class)) +
  geom_boxplot() +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(n.breaks = 27) +
  theme(legend.position = "none") +
  labs(title = "Virus SP length SignalP vs class of organism") +
  xlab("Class")

library(ggplot2)

hist(merged_taxID_organism_uniprot_signalP$SPuniprot)
hist(merged_taxID_organism_uniprot_signalP$Pos_end)

SPuniprot_virus <- merged_taxID_organism_uniprot_signalP$SPuniprot
SPuniprot_virus
table(SPuniprot_virus)
boxplot(SPuniprot_virus, na.rm = TRUE, ylab = "Length of SP", xlab = "UniProt data")
length(na.omit(SPuniprot_virus))

Pos_end_virus <- merged_taxID_organism_uniprot_signalP$Pos_end
Pos_end_virus
boxplot(Pos_end_virus, na.rm = TRUE, ylab = "Length of SP", xlab = "SignalP prediction")
length(na.omit(Pos_end_virus))

#### Ditribution in SP length in all dataset (6174 entries), SignalP prediction and UniProt data ###


png(
  filename = "Images/sp_length_viruses_uniprot_signalp_all_dist.png",
  width = 6.19,
  height = 5.88,
  units = "in",
  res = 300
)

ymax <- max(c(Pos_end_virus, SPuniprot_virus), na.rm = TRUE) * 1.25

boxplot(
  Pos_end_virus,
  SPuniprot_virus,
  na.rm = TRUE,
  names = c("SignalP prediction (n=5553)", "UniProt (n=1078)"),
  ylab = "Length of SP",
  main = "Signal peptide lengths in Virus",
  col = c("skyblue", "pink"),
  ylim = c(0, ymax),
  boxwex = 0.4
)

y <- max(c(Pos_end_virus, SPuniprot_virus), na.rm = TRUE) * 1.08

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


### statistics ###
pval <- test_all$p.value
shapiro.test(na.omit(Pos_end_virus)) # not distributed normally, distribution of signalP
shapiro.test(na.omit(SPuniprot_virus)) # not distributed normally, distribution of UniProt
hist(na.omit(Pos_end_virus))
hist(na.omit(SPuniprot_virus))

test_all <- wilcox.test(
  na.omit(Pos_end_virus),
  na.omit(SPuniprot_virus) , paired = FALSE
) ### p-value < 2.2e-16 statistically significant



library(dplyr)

#### Ditribution in SP length for confirmed (signal peptides), SignalP prediction and UniProt data ###

confirmed_merged_signalP_uniprot <- merged_taxID_organism_uniprot_signalP %>%
  filter(`SP Status` == "confirmed")


SPuniprot_virus_confirmed <- confirmed_merged_signalP_uniprot$SPuniprot
Pos_end_virus_confirmed <- confirmed_merged_signalP_uniprot$Pos_end
length(na.omit(SPuniprot_virus_confirmed))
length(na.omit(Pos_end_virus_confirmed))

png(
  filename = "Images/sp_length_viruses_uniprot_signalp_confirmed_dist.png",
  width = 5.34,
  height = 5.65,
  units = "in",
  res = 300
)


boxplot(Pos_end_virus_confirmed, SPuniprot_virus_confirmed, na.rm = TRUE, names = c("SignalP prediction (n = 180)", "UniProt (n = 195)"), ylab = "Length of SP",
        main = "Signal peptide lengths in Virus (confirmed = 214)", col = c("skyblue", "pink"),
        boxwex = 0.4)


dev.off()

text(
  x = 1.5,
  y = max(c(Pos_end_virus_confirmed, SPuniprot_virus_confirmed), na.rm = TRUE) * 0.95,
  labels = paste0("Mann-Whitney p = ", signif(pval_confirmed, 3))
)


library(dplyr)

long_peptides <- confirmed_merged_signalP_uniprot %>%
         filter(SPuniprot > 70) # only 2 proteins (confirmed) which above 70 in length


long_peptides


table(merged_taxID_organism_uniprot_signalP$`SP Status`)

#### Ditribution in SP length for potential (signal peptides), SignalP prediction and UniProt data ###


potential_merged_signalP_uniprot <- merged_taxID_organism_uniprot_signalP %>%
  filter(`SP Status` == "potential")

SPuniprot_virus_potential <- potential_merged_signalP_uniprot$SPuniprot
Pos_end_virus_potential <- potential_merged_signalP_uniprot$Pos_end
length(na.omit(SPuniprot_virus_potential))
length(na.omit(Pos_end_virus_potential))

png(
  filename = "Images/sp_length_viruses_uniprot_signalp_potential_dist.png",
  width = 6.19,
  height = 5.88,
  units = "in",
  res = 300
)

ymax <- max(c(Pos_end_virus_potential, SPuniprot_virus_potential), na.rm = TRUE) * 1.25



boxplot(Pos_end_virus_potential, SPuniprot_virus_potential, na.rm = TRUE, names = c("SignalP prediction (n = 5284)", "UniProt (n = 781)"), ylab = "Length of SP",
        main = "Signal peptide lengths in Virus (potential = 5852)", col = c("skyblue", "pink"),  boxwex = 0.4, ylim = c(0, ymax))



y <- max(c(Pos_end_virus, SPuniprot_virus), na.rm = TRUE) * 1.08

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



shapiro.test(na.omit(Pos_end_virus_potential))
shapiro.test(na.omit(SPuniprot_virus_potential)) ### the distribution is not normal
hist(na.omit(Pos_end_virus_potential))
hist(na.omit(SPuniprot_virus_potential))
summary(Pos_end_virus_potential)


range_xy <- range(
  confirmed_merged_signalP_uniprot$Pos_end,
  confirmed_merged_signalP_uniprot$SPuniprot,
  na.rm = TRUE
)

ggplot(
  data = confirmed_merged_signalP_uniprot,
  aes(x = Pos_end, y = SPuniprot)
) +
  geom_point(na.rm = TRUE) +
  labs(
    title = "Signal peptide lengths in Virus (confirmed = 214)",
    x = "SignalP prediction",
    y = "UniProt"
  ) +
  theme_classic() +
  geom_smooth(
    method = "lm",
    se = TRUE,
    na.rm = TRUE
  ) +
  coord_equal(
    xlim = range_xy,
    ylim = range_xy
  )
 
ggsave(
  filename = "Images/sp_length_viruses_uniorit_signalp.png",
  plot = last_plot(),
  width = 8,
  height = 5,
  dpi = 300
)

table(table_merged_for_analysis$ActualUniProt)

filtered_data_uniprot_yes <- table_merged_for_analysis %>%
  filter(`ActualUniProt` == "yes")
length(filtered_data_uniprot_yes)
filtered_data_uniprot_yes



filtered_data_uniprot_yes2 <- filtered_data_uniprot_yes %>%
  filter(`SP Status` == "potential")

length(filtered_data_uniprot_yes2)



### Mann -whitney test for independent samples (and not distributed normally) ####

##for potential
mntest_potential <- wilcox.test(
  na.omit(Pos_end_virus_potential),
  na.omit(SPuniprot_virus_potential), paired = FALSE,
) #### distribution is different statistically, p-value < 2.2e-16

##for confirmed
mntest_confirmed <- wilcox.test(
  na.omit(Pos_end_virus_confirmed),
  na.omit(SPuniprot_virus_confirmed) , paired = FALSE
) #### distribution is not different statistically, p-value = 0.2877

pval_potential <- mntest_potential$p.value
pval_confirmed <- mntest_confirmed$p.value


install.packages("ggpubr")
library(ggplot2)
library(ggpubr)


table(filtered_data_uniprot_yes$`SP Status`)

















