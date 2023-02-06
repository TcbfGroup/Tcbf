library(tidyverse)
profiles <- read_tsv("TAD_groups_count.tsv") |>
  select(-1)

species <- read_csv("class.csv",col_names = F)

species_order <- setNames(species$X1,
                          species$X2)
species_annotation <- data.frame(
  Species = species_order,
  Family =  species$X3
)
library(syntenet)
p1 <- plot_profiles(
  data,
  species_annotation,
  cluster_species = species_order,
  # dist_function = labdsv::dsvdis,
  # dist_params = list(index = "ruzicka")
  )