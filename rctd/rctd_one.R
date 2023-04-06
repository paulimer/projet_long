library(spacexr)
library(Matrix)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)

ref_counts <- read_csv("~/data/sc_mousebrain_allen/sc_allen/subsampled_cells_indiv_w_ensembl.csv")
metadata <- read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/metadata.csv") %>%
  select(sample_name, neighborhood_label) %>%
  mutate(neighborhood_label = str_replace_all(neighborhood_label, "/", "_"))
ref_counts <- left_join(ref_counts, metadata, by = "sample_name")
metadata_list <- ref_counts$neighborhood_label
names(metadata_list) <- ref_counts$sample_name
metadata_list <- as.factor(metadata_list)
ref_counts <- ref_counts %>% column_to_rownames("sample_name")

reference <- Reference(ref_counts, metadata_list)
