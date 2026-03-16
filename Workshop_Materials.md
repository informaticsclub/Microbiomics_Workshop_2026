# 2026 Introduction to Microbiomics
2026-03-18

## Learning objectives:
- Understanding the 16S Pipeline
- Utilizing pipeline output in R
- Creating alpha and beta diversity statistics
- Creating graphs
## R Setup (Please complete before the Workshop)
We will be using R v4.5.0 and Phyloseq 5.2.1. In RStudio create a new project named “Microbiomics_Workshop_2026”. We will use this project file for the workshop.
``` {r setup}
# Install required packages:

install.packages(c("Phyloseq", "vegan", "ape", "dplyr", "ggplot2", "ggpubr", "corncob", "patchwork"))

```
After installing all packages, restart your RStudio. This allows for all libraries to load cleanly and ensures that any cached or corrupted session data doesn’t interfere with the newly installed packages.

## Download Data
On the day of the workshop, please load your packages and download the data. The microbiome data used will come from Phyloseq in R. Additional Ways to upload data, subset, filter, and manipulate phyloseq objects can be found here: https://bioconductor.posit.co/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-basics.html

``` {r download}
# Load packages
library(c("Phyloseq", "vegan", "ape", "dplyr", "ggplot2", "ggpubr", "corncob", "patchwork"))

# Download data from the corncob package
data("ibd_phylo_sample")
data("ibd_phylo_otu")
data("ibd_phylo_taxa")

# Create our variables
metadata <- sample_data(ibd_phylo_sample)
OTU <- otu_table(ibd_phylo_otu, taxa_are_rows = TRUE)
TAX <- tax_table(ibd_phylo_taxa)

# Create our phyloseq object
ibd <- phyloseq(OTU, TAX, metadata) # optional variable TREE
ibd

# What do these variables look like?
otu_table(ibd)
tax_table(ibd)
sam_data(ibd)
```
## Exploring Data
The data we will be using is published here: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0039242
``` {r phyloseq}
# Let's look at our sample metadata
colnames(ibd@sam_data)
unique(ibd@sam_data$DiseaseState)
unique(ibd@sam_data$activity)

View(ibd@sam_data)

# Create a smaller object
set.seed(1234)
myData <- subset_samples(ibd, DiseaseState %in% c("UC", "nonIBD"))
random_samples <- sample(sample_names(myData), 10)
myData <- prune_samples(random_samples, myData)
plot_bar(myData, fill = "Phylum")

# Relative Abundance for Bar Plot
myData.rel <- transform_sample_counts(myData, function(x) x / sum(x))

# Visualize data based on Disease State
p <- plot_bar(myData.rel, fill = "Phylum")
p
p2 <- p + facet_wrap(~DiseaseState, scales= "free_x")
p2

# Visualize data on other metadata
p3 <- p + facet_wrap(gender~DiseaseState, scales= "free_x")
p3

# Visualize specific Phyla
Firm <- subset_taxa(myData.rel, Phylum == "Firmicutes")

plot_bar(Firm, fill = "Genus") +
  facet_wrap(~DiseaseState, scales= "free_x")

plot_bar(Firm, fill = "Genus") +
  facet_wrap(~DiseaseState, scales= "free_x") +
  theme(legend.text  = element_text(size = 4),
        legend.key.size = unit(0.3, "cm"))

```
## Alpha Diversity
Calculating alpha diversity requires counts data, not relative abundance. Phyloseq plot_richness() will calulate and plot Observed, Chao1, ACE, Shannon, Simpson, Inverse-Simpson, and Fisher metrics for alpha diversity of your data.
``` {r phyloseq}
plot_richness(myData)
plot_richness(myData, x = "DiseaseState")
plot_richness(myData, measures = "Shannon", x = "DiseaseState")
ad <- plot_richness(myData, measures = "Shannon", x = "DiseaseState", color = "age", shape = "gender")
ad

# What happens when we use our relative abundance data?
plot_richness(myData.rel, x = "DiseaseState")
plot_richness(myData.rel, methods = "Shannon", x = "DiseaseState")
ad.rel <- plot_richness(myData.rel, measures = "Shannon", x = "DiseaseState", color = "age", shape = "gender")
ad / ad.rel
```
## Beta Diversity and PERMONOVA
There are multiple PCoA metrics that we can use to caluclate differences betweeen groups. 
``` {r phyloseq}
# Which distance methods are there?
distanceMethodList

# Let's use Bray-Curtis Beta Diversity
Dist <- distance(myData, method="Bray")
MDS <- ordinate(myData, "PCoA", distance=Dist)
p <- plot_ordination(myData, MDS, color="DiseaseState")
p
p + stat_ellipse()

# Let's use our ordination plot to run PERMANOVA test
adonis_res <- adonis2(Dist ~ myData@sam_data$DiseaseState)
adonis_res

p2 <- p + stat_ellipse() +
  labs(title = "Beta-Diversity by Disease State",
  subtitle = paste("Bray-Curtis PERMANOVA p =", format(adonis_res$`Pr(>F)`[1], digits = 3)))

ggsave()
```
