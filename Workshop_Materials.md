# 2026 Introduction to Microbiomics
2026-03-18

## Learning objectives:
- Understand the **16S pipeline** for microbiome analysis
- Work with pipeline outputs in R using phyloseq
- Calculate **alpha** and **beta diversity** metrics
- Create clear, publication‑quality visualizations

## R Setup (Please complete before the Workshop)
We will be using phyloseq version 5.2.1 during this workshop. Before attending, please open RStudio, create a new project titled “Microbiomics_Workshop_2026”, and install the required packages listed below. We will use this project as our working environment throughout the session.
``` {r setup}
# Install required packages:

install.packages(c("Phyloseq", "vegan", "ape", "dplyr", "ggplot2", "ggpubr", "corncob", "patchwork"))

```
After installing all packages, restart your RStudio. This allows for all libraries to load cleanly and ensures that any cached or corrupted session data doesn’t interfere with the newly installed packages.

## Download Data
On the day of the workshop, please ensure that you have loaded the required R packages. The microbiome data we will be using are provided in the corncob package, specifically fecal microbiome profiles from patients with and without inflammatory bowel disease (IBD), along with other relevant patient‑level metadata. These data originate from the following publication, which examined microbial community structure in relation to IBD status: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0039242
``` {r download}
# Load packages
library(c("Phyloseq",
          "vegan",
          "ape",
          "dplyr",
          "ggplot2",
          "ggpubr",
          "corncob",
          "patchwork"))

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
To familiarize ourselves with the structure of a phyloseq object, we will examine its key components: the OTU table, the taxonomic table, and the sample metadata associated with our dataset. For the purposes of this workshop, we will work with a subset of the full dataset to create a smaller, more manageable phyloseq object. Additional ways to upload data, subset, filter, and manipulate phyloseq objects can be found here: https://bioconductor.posit.co/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-basics.html
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
In microbiome studies, **alpha diversity** refers to how many different taxa are present (richness) and how evenly those taxa are distributed (evenness) within a single ecological community or sample. Alpha diversity metrics help characterize the complexity of a microbial community by quantifying not only the total number of species or operational taxonomic units (OTUs/ASVs), but also how dominant or rare those species are relative to one another.
Calculating alpha diversity requires counts data, not relative abundance, because most diversity estimators rely on the raw frequency of observations to accurately capture richness and detect rare taxa. When relative abundances are used, essential information about sampling depth is lost, which can lead to biased or uninterpretable diversity estimates.
In phyloseq, the function plot_richness() provides a convenient way to compute and visualize multiple alpha diversity measures simultaneously. This function will calculate and plot a range of commonly used metrics, including:
- **Observed richness** – the total number of unique taxa detected.
- **Chao1** – an estimator that accounts for undetected rare taxa to predict total richness.
- **ACE** (Abundance-based Coverage Estimator) – a richness estimator emphasizing low‑abundance taxa.
- **Shannon index** – a diversity metric incorporating both richness and evenness, weighting taxa by their proportional abundance.
- **Simpson index** – a diversity metric weighted towards dominant taxa.
- **Inverse Simpson** – a reciprocal of the Simpson metric to increase interpretability (higher values reflect higher diversity).
- **Fisher’s alpha** – a richness estimator based on the log-series distribution.

Together, these metrics offer complementary perspectives on community structure, enabling a more comprehensive characterization of within-sample microbial diversity.
``` {r phyloseq}
plot_richness(myData)
plot_richness(myData, x = "DiseaseState")
plot_richness(myData, measures = "Shannon", x = "DiseaseState")
ad <- plot_richness(myData, measures = "Shannon", x = "DiseaseState", color = "age", shape = "gender")
ad

# What happens when we use our relative abundance data?
plot_richness(myData.rel, x = "DiseaseState") # this should return an error
plot_richness(myData.rel, methods = "Shannon", x = "DiseaseState")
ad.rel <- plot_richness(myData.rel, measures = "Shannon", x = "DiseaseState", color = "age", shape = "gender")
ad / ad.rel
```
## Beta Diversity and PERMONOVA
Beta diversity measures differences in microbial community composition between samples, helping identify how treatments, environments, or conditions influence overall community structure. It complements alpha diversity by focusing on between‑sample variation rather than within‑sample richness. Principal Coordinates Analysis (**PCoA**) is a common ordination method used to visualize beta diversity. It takes a distance matrix and projects samples into a reduced number of dimensions so that the distances between points approximate their ecological dissimilarity.
A widely used distance metric for microbiome studies is **Bray–Curtis** dissimilarity, which quantifies differences based on shared taxa and their abundances. Values range from 0 (identical) to 1 (completely different), making it intuitive for evaluating compositional changes.
To formally test whether groups differ in beta diversity, researchers typically use **PERMANOVA**. This non‑parametric method assesses whether the multivariate centroids of groups differ significantly using permutation-based inference. In R, the vegan package implements PERMANOVA through adonis2().
``` {r phyloseq}
# Which distance methods are available
distanceMethodList

# Let's use Bray-Curtis dissimilarity
Dist <- distance(myData, method="Bray")
MDS <- ordinate(myData, "PCoA", distance=Dist)
p <- plot_ordination(myData, MDS, color="DiseaseState")
p
p + stat_ellipse()

# Let's use our ordination plot to run PERMANOVA
adonis_res <- adonis2(Dist ~ myData@sam_data$DiseaseState)
adonis_res

p2 <- p + stat_ellipse() +
  labs(title = "Beta-Diversity by Disease State",
  subtitle = paste("Bray-Curtis PERMANOVA p =", format(adonis_res$`Pr(>F)`[1], digits = 3)))

ggsave()
```
