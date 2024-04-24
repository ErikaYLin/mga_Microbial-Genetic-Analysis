
# Load `tidyverse` and requisite packages from Bioconductor
library(tidyverse)
library(BiocStyle)
library(BiocGenerics)
library(phyloseq)

# Load packages for data visualization
library(ggplot2)
library(cowplot)
library(rphylopic)

# Data use follows the workflow by Callahan et al. (2017):
# https://doi.org/10.12688/f1000research.8986.2
# The highly-overlapping Illumina Miseq 2x250 16S amplicon sequences were made available by Schloss et al. (2012):
# https://doi.org/10.4161/gmic.21008
# The original study investigated the stability of murine gut microbiome after weaning.

# 16S Silva v132 reference library ----

# Load in the data
load(file = "vignettes/RDS/silva132_mga_results.rda")
results.df <- readRDS(file = "vignettes/RDS/results.df_silva132_mga.rds")

## Restructuring data for easier subsequent handling ----

# Extract all phyloseq objects from SILVA.OLD into a list
ps.list <- list()
for (i in 1:length(SILVA.OLD)) {
  ps.list[[i]] <- SILVA.OLD[[i]]$ps
}

# Merge phyloseq OTU tables to return single combined OTU table
ps_1 <- ps.list[[1]]
merge.otu <- merge_phyloseq(otu_table(ps_1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.list)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list[[i]]
  merge.otu <- merge_phyloseq(otu_table(ps_otu), merge.otu)
}

# Merge taxonomy
merge.tax <- merge_phyloseq(tax_table(ps_1))

for (i in 1:length(ps.list)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list[[i]]
  merge.tax <- merge_phyloseq(tax_table(ps_otu), merge.tax)
}

rownames(results.df) <- results.df$sample.ID  # Needed for to create the phyloseq-class object

# Create a new phyloseq-class object with all combined elements
# Note that the `merge_phyloseq()` does not work if phylogenetic trees have different numbers of tips

merge.ps <- phyloseq::phyloseq(otu_table(merge.otu),
                               tax_table(merge.tax),
                               sample_data(results.df))

# Colour palettes for figures
palette <- c("#c44601", "#FCC9B5", "#E1B239", "#FCF2C7", "#A3D8C6", "#329973", "#7D99E6", "#E0D2EB", "#98669F", "#353A70", "#814B08", "gray60", "black")
palette2 <- c("#005F73", "#0A9396", "#94D2BD", "#7A9EC6", "#E9D8A6", "#fbb13c", "#CA6702", "#9B2226", "#9F72AA", "#6d597a", "#355070")


## Figures of diversity measures ----

### preliminary visualization of correlation between metrics ----

# Simpson vs Shannon diversity
metric1 <-
ggplot(data = results.df, aes(x = Simpson, y = Shannon)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "A", x = "Simpson Index", y = "Shannon Index", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.background = element_blank(),
        legend.position = c(0.27,0.84))

## Shannon vs PD diversity
metric2 <-
ggplot(data = results.df, aes(x = PD, y = Shannon)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "B", x = "Phylogenetic Diversity", y = "Shannon Index", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.position = "none")

## Simpson vs PD diversity
metric3 <-
ggplot(data = results.df, aes(x = Simpson, y = PD)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "C", x = "Simpson Index", y = "Phylogenetic Diversity", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.position = "none")

## Network Connectivity vs Connectance
metric4 <-
ggplot(data = results.df, aes(x = connectance, y = connectivity)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "D", x = "Network Connectance", y = "Network Connectivity", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.position = "none")

metrics <- plot_grid(metric1, metric2, metric3, metric4, nrow = 2)
ggsave(file = "vignettes/figures/preliminary_metrics.png", plot = metrics, width = 6.86, height = 6.5)

### Shannon diversity ----

# Linear regression of Shannon diversity and age
FIT1 <- lm(Shannon ~ age, data = results.df)
summary(FIT1)

FIG1 <-
ggplot(data = results.df, aes(x = age, y = Shannon)) +
  geom_point() +  # y = c(PD, rich, Shannon, Simpson)
  geom_smooth(aes(col = "black"), method = "lm", se = F) +  # add regression lines
  labs(x = "Host Age (days)", y = "Shannon Index", fill = "Family Relationship") +
  scale_y_log10() +
  # scale_fill_manual(values = palette +
  # facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        # axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/Shannon_age.png", width = 6.86, height = 5)

## Shannon boxplot

plot1 <-
  ggplot(data = results.df, aes(x = dpwgroup, y = Shannon)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_point(shape = 1) +
  # scale_y_log10() +
  # scale_fill_manual(values = palette) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Shannon Index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/Shannon_boxplot_silva132.png", width = 6.86, height = 5)

### Simpson diversity ----

FIG2 <-
ggplot(data = results.df) +
  geom_point(aes(x = age, y = Simpson)) +  # y = c(PD, rich, Shannon, Simpson)
  labs(x = "Age (days)", y = "Shannon Index", fill = "Family Relationship") +
  scale_y_log10() +
  # scale_fill_manual(values = palette) +
  # facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        # axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("figures/Simpson_age.png", width = 6.86, height = 5)

plot2 <-
  ggplot(data = results.df, aes(x = dpwgroup, y = Simpson)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_point(shape = 1) +
  # scale_y_log10() +
  # scale_fill_manual(values = palette) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Simpson Index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/Simpson_boxplot_silva132.png", width = 6.86, height = 5)

### Phylogenetic diversity ----
FIG3 <-
  ggplot(data = results.df) +
  geom_point(aes(x = age, y = Shannon)) +  # y = c(PD, rich, Shannon, Simpson)
  labs(x = "Age (days)", y = "Shannon Index", fill = "Family Relationship") +
  scale_y_log10() +
  # scale_fill_manual(values = palette +
  # facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        # axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

plot3 <-
  ggplot(data = results.df, aes(x = dpwgroup, y = PD)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_jitter(shape = 1, width = 0.2) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Phylogenetic Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/PD_boxplot_silva132.png", width = 6.86, height = 5)

### Species richness ----
plot4 <-
  ggplot(data = results.df, aes(x = dpwgroup, y = rich)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5, outlier.shape = 1) +  # y = c(PD, rich, Shannon, Simpson)
  geom_point(shape = 1) +
  # scale_y_log10() +
  # scale_fill_manual(values = palette) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Species Richness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/rich_boxplot_silva132.png", width = 6.86, height = 5)


## Abundance bar plots ----

# Abundance plots were built following the workflow of Hui (2021):
# https://www.yanh.org/2021/01/01/microbiome-r/#abundance-bar-plot

### Class ----
# Calculate relative taxa abundance
ps.rel <- transform_sample_counts(merge.ps, function(x) x/sum(x)*100)
# Agglomerate samples by taxon of choice
agglomerated <- tax_glom(ps.rel, taxrank = 'Class', NArm = FALSE)
# Melt into data frame
ps.melt <- psmelt(agglomerated)

ps.melt <- ps.melt %>%
  group_by(Class) %>%
  mutate(median = median(Abundance))
rare <- unique(ps.melt$Class[ps.melt$median > 1])
ps.melt$Class[!(ps.melt$Class %in% rare)] <- "< 1%"
ps.sum <- ps.melt %>%
  group_by(sample.ID, dpwgroup, Class) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial classes are <1%
other <- ps.melt$Class[!(ps.melt$Class %in% rare)]
length(other)

# Abundance bar plot (Class)
FIG <-
  ggplot(ps.sum, aes(x = sample.ID, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", aes(fill = Class)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
    axis.text.y = element_text(size = 8, family = "sans"),
    axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
    axis.line.y = element_line(linewidth = 0.2),
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks.x = element_line(linewidth = 0.2),
    legend.text = element_text(size = 7, family = "sans", face = "bold"),
    legend.title = element_text(size = 8, family = "sans", face = "bold"),
    legend.key.height = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.15, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("vignettes/figures/abundance_silva132_class.png", plot = FIG, width = 6.86, height = 4.5)


### Family ----
# Agglomerate samples by taxon of choice
agglomerated2 <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)
# Melt into data frame
ps.melt2 <- psmelt(agglomerated2)

ps.melt2 <- ps.melt2 %>%
  group_by(Family) %>%
  mutate(median = median(Abundance))
rare2 <- unique(ps.melt2$Family[ps.melt2$median > 1])
ps.melt2$Family[!(ps.melt2$Family %in% rare2)] <- "< 1%"
ps.sum2 <- ps.melt2 %>%
  group_by(sample.ID, dpwgroup, Family) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial families are <1%
other2 <- ps.melt2$Family[!(ps.melt2$Family %in% rare2)]
length(other2)

# Abundance bar plot
FIG.fam <-
  ggplot(ps.sum2, aes(x = sample.ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", aes(fill = Family)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("vignettes/figures/abundance_silva132_Family.png", plot = FIG.fam, width = 6.86, height = 4.5)

# Retrieve mouse phylopic
mouse <- get_phylopic("92989e35-4e68-4a2d-b3a2-191ba9da671a")

icons <-
  ggplot() +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  add_phylopic(mouse, alpha = 1, x = 24.2, y = 96, ysize = 4.5, color = "#F8766D") +
  add_phylopic(mouse, alpha = 1, x = 62, y = 96.7, ysize = 7, color = "#00BFC4") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())

# Layer mouse icons on top of abundance bar plot
# ABUNDANCE <-
#   ggdraw(FIG) +
#   draw_plot(icons,
#             x = 0,
#             y = 0,
#             width = 1,
#             height = 1)
# ggsave("vignettes/figures/abundance_barplot_silva132.png", plot = ABUNDANCE, width = 6.86, height = 4.5)

ABUNDANCE2 <-
  ggdraw(FIG.fam) +
  draw_plot(icons,
            x = 0,
            y = 0,
            width = 0.96,
            height = 1)
ggsave("vignettes/figures/abundance_barplot2_silva132.png", plot = ABUNDANCE2, width = 6.86, height = 4.5)


# 16S Silva NR99 v138.1 reference library with species ----

# Load in the data
load(file = "vignettes/RDS/silva138.1S_mga_results.rda")
results.df2 <- readRDS(file = "vignettes/RDS/results.df_silva138.1S_mga.rds")

## Restructuring data for easier subsequent handling ----

# Extract all phyloseq objects from SILVA.NEWS into a list
ps.list2 <- list()
for (i in 1:length(SILVA.NEWS)) {
  ps.list2[[i]] <- SILVA.NEWS[[i]]$ps
}

# Merge phyloseq OTU tables to return single combined OTU table
ps_1 <- ps.list2[[1]]
merge.otu <- merge_phyloseq(otu_table(ps_1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.list2)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list2[[i]]
  merge.otu <- merge_phyloseq(otu_table(ps_otu), merge.otu)
}

# Merge taxonomy
merge.tax <- merge_phyloseq(tax_table(ps_1))

for (i in 1:length(ps.list2)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list2[[i]]
  merge.tax <- merge_phyloseq(tax_table(ps_otu), merge.tax)
}

rownames(results.df2) <- results.df2$sample.ID  # Needed for to create the phyloseq-class object

# Create a new phyloseq-class object with all combined elements
# Note that the `merge_phyloseq()` does not work if phylogenetic trees have different numbers of tips

merge.ps2 <- phyloseq::phyloseq(otu_table(merge.otu),
                               tax_table(merge.tax),
                               sample_data(results.df2))


## Figures of diversity measures ----

### preliminary visualization of correlation between metrics ----

# Simpson vs Shannon diversity
metric1 <-
  ggplot(data = results.df2, aes(x = Simpson, y = Shannon)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "A", x = "Simpson Index", y = "Shannon Index", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.background = element_blank(),
        legend.position = c(0.27,0.84))

## Shannon vs PD diversity
metric2 <-
  ggplot(data = results.df2, aes(x = PD, y = Shannon)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "B", x = "Phylogenetic Diversity", y = "Shannon Index", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.position = "none")

## Simpson vs PD diversity
metric3 <-
  ggplot(data = results.df2, aes(x = Simpson, y = PD)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "C", x = "Simpson Index", y = "Phylogenetic Diversity", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.position = "none")

## Network Connectivity vs Connectance
metric4 <-
  ggplot(data = results.df2, aes(x = connectance, y = connectivity)) +
  geom_point(aes(col = dpwgroup), size = 2) +
  geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
  labs(tag = "D", x = "Network Connectance", y = "Network Connectivity", col = "Days Post Weaning") +
  theme_classic() +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.position = "none")

metrics <- plot_grid(metric1, metric2, metric3, metric4, nrow = 2)
ggsave(file = "vignettes/figures/preliminary_metrics2.png", plot = metrics, width = 6.86, height = 6.5)

### Shannon diversity ----

# Linear regression of Shannon diversity and age
FIT1 <- lm(Shannon ~ age, data = results.df2)
summary(FIT1)

FIG1 <-
  ggplot(data = results.df2, aes(x = age, y = Shannon)) +
  geom_point() +  # y = c(PD, rich, Shannon, Simpson)
  geom_smooth(aes(col = "black"), method = "lm", se = F) +  # add regression lines
  labs(x = "Host Age (days)", y = "Shannon Index", fill = "Family Relationship") +
  scale_y_log10() +
  # scale_fill_manual(values = palette +
  # facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        # axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

# ggsave("vignettes/figures/Shannon_age.png", width = 6.86, height = 5)

## Shannon boxplot

plot1 <-
  ggplot(data = results.df2, aes(x = dpwgroup, y = Shannon)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_point(shape = 1) +
  # scale_y_log10() +
  # scale_fill_manual(values = palette) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Shannon Index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/Shannon_boxplot_silva138.1S.png", width = 6.86, height = 5)

### Simpson diversity ----

FIG2 <-
  ggplot(data = results.df2) +
  geom_point(aes(x = age, y = Simpson)) +  # y = c(PD, rich, Shannon, Simpson)
  labs(x = "Age (days)", y = "Shannon Index", fill = "Family Relationship") +
  scale_y_log10() +
  # scale_fill_manual(values = palette) +
  # facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        # axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("figures/Simpson_age.png", width = 6.86, height = 5)

plot2 <-
  ggplot(data = results.df2, aes(x = dpwgroup, y = Simpson)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_point(shape = 1) +
  # scale_y_log10() +
  # scale_fill_manual(values = palette) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Simpson Index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/Simpson_boxplot_silva138.1S.png", width = 6.86, height = 5)

### Phylogenetic diversity ----
FIG3 <-
  ggplot(data = results.df2) +
  geom_point(aes(x = age, y = Shannon)) +  # y = c(PD, rich, Shannon, Simpson)
  labs(x = "Age (days)", y = "Shannon Index", fill = "Family Relationship") +
  scale_y_log10() +
  # scale_fill_manual(values = palette +
  # facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        # axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

plot3 <-
  ggplot(data = results.df2, aes(x = dpwgroup, y = PD)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_jitter(shape = 1, width = 0.2) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Phylogenetic Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/PD_boxplot_silva138.1S.png", width = 6.86, height = 5)

### Species richness ----
plot4 <-
  ggplot(data = results.df2, aes(x = dpwgroup, y = rich)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5, outlier.shape = 1) +  # y = c(PD, rich, Shannon, Simpson)
  geom_point(shape = 1) +
  # scale_y_log10() +
  # scale_fill_manual(values = palette) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Species Richness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/rich_boxplot_silva138.1S.png", width = 6.86, height = 5)


## Abundance bar plots ----

# Abundance plots were built following the workflow of Hui (2021):
# https://www.yanh.org/2021/01/01/microbiome-r/#abundance-bar-plot

### Class ----
# Calculate relative taxa abundance
ps.rel <- transform_sample_counts(merge.ps2, function(x) x/sum(x)*100)
# Agglomerate samples by taxon of choice
agglomerated <- tax_glom(ps.rel, taxrank = 'Class', NArm = FALSE)
# Melt into data frame
ps.melt <- psmelt(agglomerated)

ps.melt <- ps.melt %>%
  group_by(Class) %>%
  mutate(median = median(Abundance))
rare <- unique(ps.melt$Class[ps.melt$median > 1])
ps.melt$Class[!(ps.melt$Class %in% rare)] <- "< 1%"
ps.sum <- ps.melt %>%
  group_by(sample.ID, dpwgroup, Class) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial classes are <1%
other <- ps.melt$Class[!(ps.melt$Class %in% rare)]
length(other)

# Abundance bar plot (Class)
FIG <-
  ggplot(ps.sum, aes(x = sample.ID, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", aes(fill = Class)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("vignettes/figures/abundance_silva138.1S_class.png", plot = FIG, width = 6.86, height = 4.5)


### Family ----
# Agglomerate samples by taxon of choice
agglomerated2 <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)
# Melt into data frame
ps.melt2 <- psmelt(agglomerated2)

ps.melt2 <- ps.melt2 %>%
  group_by(Family) %>%
  mutate(median = median(Abundance))
rare2 <- unique(ps.melt2$Family[ps.melt2$median > 1])
ps.melt2$Family[!(ps.melt2$Family %in% rare2)] <- "< 1%"
ps.sum2 <- ps.melt2 %>%
  group_by(sample.ID, dpwgroup, Family) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial families are <1%
other2 <- ps.melt2$Family[!(ps.melt2$Family %in% rare2)]
length(other2)

# Abundance bar plot
FIG.fam <-
  ggplot(ps.sum2, aes(x = sample.ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", aes(fill = Family)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("vignettes/figures/abundance_silva138.1S_Family.png", plot = FIG.fam, width = 6.86, height = 4.5)

# Layer mouse icons on top of abundance bar plot
# ABUNDANCE <-
#   ggdraw(FIG) +
#   draw_plot(icons,
#             x = 0,
#             y = 0,
#             width = 1,
#             height = 1)
# ggsave("vignettes/figures/abundance_barplot_silva138.1S.png", plot = ABUNDANCE, width = 6.86, height = 4.5)

ABUNDANCE2 <-
  ggdraw(FIG.fam) +
  draw_plot(icons,
            x = 0,
            y = 0,
            width = 0.96,
            height = 1)
ggsave("vignettes/figures/abundance_barplot2_silva138.1S.png", plot = ABUNDANCE2, width = 6.86, height = 4.5)


# 16S Silva NR99 v138.1 reference library (no species, no network) ----

# Load in the data
load(file = "vignettes/RDS/silva138.1_mga_noNet.rda")
results.df3 <- readRDS(file = "vignettes/RDS/results.df_silva138.1_mga_noNet.rds")

## Restructuring data for easier subsequent handling ----

# Extract all phyloseq objects from SILVA.NEW into a list
ps.list3 <- list()
for (i in 1:length(SILVA.NEW)) {
  ps.list3[[i]] <- SILVA.NEW[[i]]$ps
}

# Merge phyloseq OTU tables to return single combined OTU table
ps_1 <- ps.list3[[1]]
merge.otu <- merge_phyloseq(otu_table(ps_1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.list3)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list3[[i]]
  merge.otu <- merge_phyloseq(otu_table(ps_otu), merge.otu)
}

# Merge taxonomy
merge.tax <- merge_phyloseq(tax_table(ps_1))

for (i in 1:length(ps.list3)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list3[[i]]
  merge.tax <- merge_phyloseq(tax_table(ps_otu), merge.tax)
}

rownames(results.df3) <- results.df3$sample.ID  # Needed for to create the phyloseq-class object

# Create a new phyloseq-class object with all combined elements
# Note that the `merge_phyloseq()` does not work if phylogenetic trees have different numbers of tips

merge.ps3 <- phyloseq::phyloseq(otu_table(merge.otu),
                               tax_table(merge.tax),
                               sample_data(results.df3))


# ## Figures of diversity measures
#
# ### preliminary visualization of correlation between metrics
#
# # Simpson vs Shannon diversity
# metric1 <-
#   ggplot(data = results.df3, aes(x = Simpson, y = Shannon)) +
#   geom_point(aes(col = dpwgroup), size = 2) +
#   geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
#   labs(tag = "A", x = "Simpson Index", y = "Shannon Index", col = "Days Post Weaning") +
#   theme_classic() +
#   theme(plot.tag = element_text(size = 14, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.title = element_text(size = 8, face = "bold"),
#         legend.background = element_blank(),
#         legend.position = c(0.27,0.84))
#
# ## Shannon vs PD diversity
# metric2 <-
#   ggplot(data = results.df3, aes(x = PD, y = Shannon)) +
#   geom_point(aes(col = dpwgroup), size = 2) +
#   geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
#   labs(tag = "B", x = "Phylogenetic Diversity", y = "Shannon Index", col = "Days Post Weaning") +
#   theme_classic() +
#   theme(plot.tag = element_text(size = 14, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.title = element_text(size = 8, face = "bold"),
#         legend.position = "none")
#
# ## Simpson vs PD diversity
# metric3 <-
#   ggplot(data = results.df3, aes(x = Simpson, y = PD)) +
#   geom_point(aes(col = dpwgroup), size = 2) +
#   geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
#   labs(tag = "C", x = "Simpson Index", y = "Phylogenetic Diversity", col = "Days Post Weaning") +
#   theme_classic() +
#   theme(plot.tag = element_text(size = 14, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.title = element_text(size = 8, face = "bold"),
#         legend.position = "none")
#
# ## Network Connectivity vs Connectance
# metric4 <-
#   ggplot(data = results.df3, aes(x = connectance, y = connectivity)) +
#   geom_point(aes(col = dpwgroup), size = 2) +
#   geom_smooth(aes(col = dpwgroup), method = "lm", se = F, linewidth = 1) +
#   labs(tag = "D", x = "Network Connectance", y = "Network Connectivity", col = "Days Post Weaning") +
#   theme_classic() +
#   theme(plot.tag = element_text(size = 14, face = "bold"),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.title = element_text(size = 8, face = "bold"),
#         legend.position = "none")
#
# metrics <- plot_grid(metric1, metric2, metric3, metric4, nrow = 2)
# ggsave(file = "vignettes/figures/preliminary_metrics.png", plot = metrics, width = 6.86, height = 6.5)
#
# ### Shannon diversity
#
# # Linear regression of Shannon diversity and age
# FIT1 <- lm(Shannon ~ age, data = results.df3)
# summary(FIT1)
#
# FIG1 <-
#   ggplot(data = results.df3, aes(x = age, y = Shannon)) +
#   geom_point() +  # y = c(PD, rich, Shannon, Simpson)
#   geom_smooth(aes(col = "black"), method = "lm", se = F) +  # add regression lines
#   labs(x = "Host Age (days)", y = "Shannon Index", fill = "Family Relationship") +
#   scale_y_log10() +
#   # scale_fill_manual(values = palette +
#   # facet_wrap(vars(site)) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 8, family = "sans"),
#         axis.text.x = element_text(size = 8, family = "sans"),
#         # axis.ticks.x = element_line(colour = "transparent"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         panel.spacing = unit(0.2, "lines"),
#         strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# ggsave("vignettes/figures/Shannon_age.png", width = 6.86, height = 5)
#
# ## Shannon boxplot
#
# plot1 <-
#   ggplot(data = results.df3, aes(x = dpwgroup, y = Shannon)) +
#   geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
#   geom_point(shape = 1) +
#   # scale_y_log10() +
#   # scale_fill_manual(values = palette) +
#   labs(fill = NULL, x = "Days Post Weaning", y = "Shannon Index") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.tag = element_text(size = 12, family = "sans", face = "bold"),
#         # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 10, family = "sans"),
#         axis.text.x = element_text(size = 10, family = "sans"),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# ggsave("vignettes/figures/Shannon_boxplot.png", width = 6.86, height = 5)
#
# ### Simpson diversity
#
# FIG2 <-
#   ggplot(data = results.df3) +
#   geom_point(aes(x = age, y = Simpson)) +  # y = c(PD, rich, Shannon, Simpson)
#   labs(x = "Age (days)", y = "Shannon Index", fill = "Family Relationship") +
#   scale_y_log10() +
#   # scale_fill_manual(values = palette) +
#   # facet_wrap(vars(site)) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 8, family = "sans"),
#         axis.text.x = element_text(size = 8, family = "sans"),
#         # axis.ticks.x = element_line(colour = "transparent"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         panel.spacing = unit(0.2, "lines"),
#         strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# ggsave("figures/Simpson_age.png", width = 6.86, height = 5)
#
# plot2 <-
#   ggplot(data = results.df3, aes(x = dpwgroup, y = Simpson)) +
#   geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
#   geom_point(shape = 1) +
#   # scale_y_log10() +
#   # scale_fill_manual(values = palette) +
#   labs(fill = NULL, x = "Days Post Weaning", y = "Simpson Index") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.tag = element_text(size = 12, family = "sans", face = "bold"),
#         # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 10, family = "sans"),
#         axis.text.x = element_text(size = 10, family = "sans"),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
#         legend.position = "none",
#         panel.background = element_blank(),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         plot.background = element_blank(),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# ggsave("vignettes/figures/Simpson_boxplot.png", width = 6.86, height = 5)
#
# ### Phylogenetic diversity
# FIG3 <-
#   ggplot(data = results.df3) +
#   geom_point(aes(x = age, y = Shannon)) +  # y = c(PD, rich, Shannon, Simpson)
#   labs(x = "Age (days)", y = "Shannon Index", fill = "Family Relationship") +
#   scale_y_log10() +
#   # scale_fill_manual(values = palette +
#   # facet_wrap(vars(site)) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 8, family = "sans"),
#         axis.text.x = element_text(size = 8, family = "sans"),
#         # axis.ticks.x = element_line(colour = "transparent"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         panel.spacing = unit(0.2, "lines"),
#         strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# plot3 <-
#   ggplot(data = results.df3, aes(x = dpwgroup, y = PD)) +
#   geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
#   geom_point(shape = 1) +
#   # scale_y_log10() +
#   # scale_fill_manual(values = palette) +
#   labs(fill = NULL, x = "Days Post Weaning", y = "Phylogenetic Diversity") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.tag = element_text(size = 12, family = "sans", face = "bold"),
#         # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 10, family = "sans"),
#         axis.text.x = element_text(size = 10, family = "sans"),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
#         legend.position = "none",
#         panel.background = element_blank(),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         plot.background = element_blank(),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# ggsave("vignettes/figures/PD_boxplot.png", width = 6.86, height = 5)
#
# ### Species richness
# plot3 <-
#   ggplot(data = results.df3, aes(x = dpwgroup, y = rich)) +
#   geom_boxplot(aes(fill = dpwgroup), size = 0.5, outlier.shape = 1) +  # y = c(PD, rich, Shannon, Simpson)
#   geom_point(shape = 1) +
#   # scale_y_log10() +
#   # scale_fill_manual(values = palette) +
#   labs(fill = NULL, x = "Days Post Weaning", y = "Species Richness") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.tag = element_text(size = 12, family = "sans", face = "bold"),
#         # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
#         axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
#         axis.text.y = element_text(size = 10, family = "sans"),
#         axis.text.x = element_text(size = 10, family = "sans"),
#         axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
#         legend.position = "none",
#         panel.background = element_blank(),
#         panel.border = element_rect(linewidth = 0.2, fill = NA),
#         plot.background = element_blank(),
#         plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left
#
# ggsave("vignettes/figures/rich_boxplot.png", width = 6.86, height = 5)


## Abundance bar plots ----

# Abundance plots were built following the workflow of Hui (2021):
# https://www.yanh.org/2021/01/01/microbiome-r/#abundance-bar-plot

### Class ----
# Calculate relative taxa abundance
ps.rel <- transform_sample_counts(merge.ps3, function(x) x/sum(x)*100)
# Agglomerate samples by taxon of choice
agglomerated <- tax_glom(ps.rel, taxrank = 'Class', NArm = FALSE)
# Melt into data frame
ps.melt <- psmelt(agglomerated)

ps.melt <- ps.melt %>%
  group_by(Class) %>%
  mutate(median = median(Abundance))
rare <- unique(ps.melt$Class[ps.melt$median > 1])
ps.melt$Class[!(ps.melt$Class %in% rare)] <- "< 1%"
ps.sum <- ps.melt %>%
  group_by(sample.ID, dpwgroup, Class) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial classes are <1%
other <- ps.melt$Class[!(ps.melt$Class %in% rare)]
length(other)

# Abundance bar plot (Class)
FIG <-
  ggplot(ps.sum, aes(x = sample.ID, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", aes(fill = Class)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("vignettes/figures/abundance_bacterial_class.png", plot = FIG, width = 6.86, height = 4.5)


### Family ----
# Agglomerate samples by taxon of choice
agglomerated2 <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)
# Melt into data frame
ps.melt2 <- psmelt(agglomerated2)

ps.melt2 <- ps.melt2 %>%
  group_by(Family) %>%
  mutate(median = median(Abundance))
rare2 <- unique(ps.melt2$Family[ps.melt2$median > 1])
ps.melt2$Family[!(ps.melt2$Family %in% rare2)] <- "< 1%"
ps.sum2 <- ps.melt2 %>%
  group_by(sample.ID, dpwgroup, Family) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial families are <1%
other2 <- ps.melt2$Family[!(ps.melt2$Family %in% rare2)]
length(other2)

# Abundance bar plot
FIG.fam3 <-
  ggplot(ps.sum2, aes(x = sample.ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", aes(fill = Family)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("vignettes/figures/abundance_silva138.1_Family.png", plot = FIG.fam3, width = 6.86, height = 4.5)

# Layer mouse icons on top of abundance bar plot
# ABUNDANCE <-
#   ggdraw(FIG) +
#   draw_plot(icons,
#             x = 0,
#             y = 0,
#             width = 1,
#             height = 1)
# ggsave("vignettes/figures/abundance_barplot_bacteria.png", plot = ABUNDANCE, width = 6.86, height = 4.5)

ABUNDANCE2 <-
  ggdraw(FIG.fam3) +
  draw_plot(icons,
            x = 0,
            y = 0,
            width = 0.96,
            height = 1)
ggsave("vignettes/figures/abundance_barplot_silva138.1.png", plot = ABUNDANCE2, width = 6.86, height = 4.5)


# 16S Greengenes v13.8 reference library ----

# Load in the data
load(file = "vignettes/RDS/gg13.8_mga_results.rda")
results.df4 <- readRDS(file = "vignettes/RDS/results.df_gg13.8_mga.rds")

## Restructuring data for easier subsequent handling ----

# Extract all phyloseq objects from GREEN into a list
ps.list4 <- list()
for (i in 1:length(GREEN)) {
  ps.list4[[i]] <- GREEN[[i]]$ps
}

# Merge phyloseq OTU tables to return single combined OTU table
ps_1 <- ps.list4[[1]]
merge.otu <- merge_phyloseq(otu_table(ps_1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.list4)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list4[[i]]
  merge.otu <- merge_phyloseq(otu_table(ps_otu), merge.otu)
}

# Merge taxonomy
merge.tax <- merge_phyloseq(tax_table(ps_1))

for (i in 1:length(ps.list4)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list4[[i]]
  merge.tax <- merge_phyloseq(tax_table(ps_otu), merge.tax)
}

rownames(results.df4) <- results.df4$sample.ID  # Needed for to create the phyloseq-class object

# Create a new phyloseq-class object with all combined elements
# Note that the `merge_phyloseq()` does not work if phylogenetic trees have different numbers of tips

merge.ps4 <- phyloseq::phyloseq(otu_table(merge.otu),
                                tax_table(merge.tax),
                                sample_data(results.df4))


## Abundance bar plots ----

# Abundance plots were built following the workflow of Hui (2021):
# https://www.yanh.org/2021/01/01/microbiome-r/#abundance-bar-plot

### Family ----
# Calculate relative taxa abundance
ps.rel <- transform_sample_counts(merge.ps3, function(x) x/sum(x)*100)
# Agglomerate samples by taxon of choice
agglomerated2 <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)
# Melt into data frame
ps.melt2 <- psmelt(agglomerated2)

ps.melt2 <- ps.melt2 %>%
  group_by(Family) %>%
  mutate(median = median(Abundance))
rare2 <- unique(ps.melt2$Family[ps.melt2$median > 1])
ps.melt2$Family[!(ps.melt2$Family %in% rare2)] <- "< 1%"
ps.sum2 <- ps.melt2 %>%
  group_by(sample.ID, dpwgroup, Family) %>%
  summarise(Abundance = sum(Abundance))

# How many bacterial families are <1%
other2 <- ps.melt2$Family[!(ps.melt2$Family %in% rare2)]
length(other2)

# Abundance bar plot
FIG.fam4 <-
  ggplot(ps.sum2, aes(x = sample.ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", aes(fill = Family)) +
  labs(x = "Days Post Weaning", y= "Abundance (%)") +
  facet_wrap(~dpwgroup, scales = "free_x", strip.position = "bottom") + # nrow = 1
  theme_classic() +
  scale_fill_manual(values = palette2) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 45, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 8, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

# ggsave("vignettes/figures/abundance_gg13.8_Family.png", plot = FIG.fam3, width = 6.86, height = 4.5)

ABUNDANCE2 <-
  ggdraw(FIG.fam4) +
  draw_plot(icons,
            x = 0,
            y = 0,
            width = 0.96,
            height = 1)
ggsave("vignettes/figures/abundance_barplot_gg13.8.png", plot = ABUNDANCE2, width = 6.86, height = 4.5)


## Figures of diversity measures ----

### Phylogenetic diversity ----
plot3 <-
  ggplot(data = results.df4, aes(x = dpwgroup, y = PD)) +
  geom_boxplot(aes(fill = dpwgroup), size = 0.5) +  # y = c(PD, rich, Shannon, Simpson)
  geom_jitter(shape = 1, width = 0.2) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Phylogenetic Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/PD_boxplot_gg13.8.png", width = 6.86, height = 5)


# Sensitivity Analysis: Reference Libraries ----

# Remove row names from results data frame
div.metrics <- results.df[,c(1:5,27)]
rownames(div.metrics) <- NULL
div.metrics2 <- results.df2[,c(1:5,27)]
rownames(div.metrics2) <- NULL
div.metrics4 <- results.df4[,c(1:5,27)]
rownames(div.metrics4) <- NULL

# Extract diversity indices and dpw group for each ref library
PD.df <- list(Silva_old = div.metrics,
              Silva_new = div.metrics2,
              Greengenes = div.metrics4)

# Add ref libraries
PD.df$Silva_old$ref <- rep("Silva_132", 19)
PD.df$Silva_new$ref <- rep("Silva_138.1", 19)
PD.df$Greengenes$ref <- rep("Greengenes_13.8", 19)

# Collapse into data frame
PD.df <- do.call(rbind, PD.df)
rownames(PD.df) <- NULL # remove row names

# Boxplot of PD across different ref libraries
PD_box <- ggplot(data = PD.df, aes(x = dpwgroup, y = PD)) +
  geom_boxplot(size = 0.5, aes(fill = dpwgroup)) +
  # scale_fill_manual(values = palette3) +
  geom_jitter(shape = 1, width = 0.2) +
  facet_wrap(~ref) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Phylogenetic Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1),
        strip.text = element_text(size = 12, family = "sans"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/PD_boxplot_sensitivity.png", plot = PD_box, width = 6.86, height = 5.2)


# Boxplot of Shannon diversity across different ref libraries
Shannon_box <- ggplot(data = PD.df, aes(x = dpwgroup, y = Shannon)) +
  geom_boxplot(size = 0.5, aes(fill = dpwgroup)) +
  # scale_fill_manual(values = palette3) +
  geom_jitter(shape = 1, width = 0.2) +
  facet_wrap(~ref) +
  labs(fill = NULL, x = "Days Post Weaning", y = "Shannon Index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12, family = "sans", face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.text.x = element_text(size = 10, family = "sans"),
        axis.title.x = element_text(size = 12, family = "sans", face = "bold", vjust = 1),
        strip.text = element_text(size = 12, family = "sans"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/Shannon_boxplot_sensitivity.png", plot = Shannon_box, width = 6.86, height = 5.2)


# Single sample test for species agglomeration ----

# Load in the data
load(file = "vignettes/RDS/TEST_mga_results.rda")
results.df5 <- readRDS(file = "vignettes/RDS/results.df_TEST_mga.rds")

# Load in the data
load(file = "vignettes/RDS/TEST2_mga_results.rda")
results.df6 <- readRDS(file = "vignettes/RDS/results.df_TEST2_mga.rds")

phy <- TEST[[1]]$ps

glom_test <- as.data.frame(tax_table(phy))
glom_test$Species <- paste0(glom_test$Genus, glom_test$Species, sep = " ")

phy2 <- TEST[[1]]$ps
tax_table(phy2) <- tax_table(as.matrix(glom_test))

phy_glom2 <- tax_glom(phy2, taxrank = "Species", NArm = FALSE)


# Google Scholar search ----

google <- data.frame(Year = seq(from = 1990, to = 2019, by = 1),
                     Articles = c(47.2, 41.1, 42.6, 45.9, 48.5, 56.5,
                                  60.1, 67.2, 76.2, 76.7, 111, 108,
                                  122, 129, 152, 167, 186, 220, 239, 299,
                                  363, 380, 454, 506, 529, 541, 585, 562,
                                  627, 614))
# google$Year <- as.factor(google$Year)

# line graph of articles per year
line.gs <- ggplot(data = google, aes(x = Year, y = Articles)) +
  geom_line(col = "#D19651", linewidth = 2) +
  scale_x_continuous(limits = c(1990, 2020.5)) +
  labs(y = "Articles (thousands)",
       title = "Microbial Research Articles Published Per Year") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 16, family = "sans", face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 14, family = "sans"),
        axis.text.x = element_text(size = 14, family = "sans"),
        axis.title.x = element_text(size = 16, family = "sans", face = "bold", vjust = 1),
        strip.text = element_text(size = 16, family = "sans"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("vignettes/figures/articles_lineplot.png", plot = line.gs, width = 6.86, height = 5.8)
