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
) ### for verification if human prresent like a host

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
