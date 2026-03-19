####Setup####
# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(phyloseq)
library(ggplot2)


# Set data directory to bracken output location
data_dir <- "data/processed/bracken"

####Read bracken files####
files <- list.files(data_dir, pattern = "_bracken.txt", full.names = TRUE)

# Read and store in a list
bracken_list <- lapply(files, function(file) {
  df <- read_tsv(file)
  
  # Extract sample name from filename
  sample_name <- gsub("_bracken.txt", "", basename(file))
  
  df$sample <- sample_name
  return(df)
})

# Combine all into one dataframe
bracken_df <- bind_rows(bracken_list)

head(bracken_df)

####Create abundance matrix####
# rows - species, columns = samples, values = counts
abundance_table <- bracken_df %>%
  select(name, sample, new_est_reads) %>%
  pivot_wider(names_from = sample, values_from = new_est_reads, values_fill = 0)

# Check structure
dim(abundance_table)
head(abundance_table)

# Set rownames
abundance_mat <- as.data.frame(abundance_table)
rownames(abundance_mat) <- abundance_mat$name
abundance_mat$name <- NULL

write.csv(abundance_mat, "data/processed/species_abundance_matrix.csv")

####Create metadata####
metadata <- data.frame(
  sample_id = c(
    "SRR8146935", "SRR8146936", "SRR8146938",
    "SRR8146951", "SRR8146952", "SRR8146954"
  ),
  diet = c(
    "omnivore", "omnivore", "omnivore",
    "vegan", "vegan", "vegan"
  )
)

rownames(metadata) <- metadata$sample_id
metadata$sample_id <- NULL

write.csv(metadata, "data/processed/sample_metadata.csv")

# Check metadata to abundance matrix
all(colnames(abundance_mat) == rownames(metadata))

####Build phyloseq object####
otu <- otu_table(as.matrix(abundance_mat), taxa_are_rows = TRUE)
samples <- sample_data(metadata)

ps <- phyloseq(otu, samples)

# Verify
ps

####Transform to relative abundance####
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

####Melt phyloseq object to long format####
plot_df <- psmelt(ps_rel)

# Check columns
colnames(plot_df)
head(plot_df)

####Plot stacked bar chart####
# Keep only top 20 most abundant species
top_taxa <- plot_df %>%
  group_by(OTU) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:20) %>%
  pull(OTU)

plot_top <- plot_df %>%
  mutate(Taxon = ifelse(OTU %in% top_taxa, OTU, "Other"))

# Create bar chart
ggplot(plot_top, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ diet, scales = "free_x", space = "free_x") +
  labs(
    title = "Top 20 species by relative abundance",
    x = "Sample",
    y = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save figure
ggsave(
  filename = "results/species_relative_abundance_barplot.png",
  width = 12,
  height = 6,
  dpi = 300
)
