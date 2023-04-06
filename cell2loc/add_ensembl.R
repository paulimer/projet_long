library(BiocManager)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(biomaRt)

correspondance_ensembl <- function(mgi, ref_df) {
  matched <- left_join(tibble("mgi_symbol" = mgi), ref_df, by = "mgi_symbol", multiple = "first") %>%
    replace_na(list(mgi_symbol = "uninterested", ensembl_gene_id = "problem"))
  matched$ensembl_gene_id <- make.unique(matched$ensembl_gene_id, sep = "_")
  matched %>% pull("ensembl_gene_id")
}

subs_cells <- read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/subsampled_cells_indiv.csv")

ensembl.database <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
ensembl.attributes <- c('ensembl_gene_id', 'mgi_symbol')
ensembl_cor <- getBM(attributes = ensembl.attributes, mart = ensembl.database)
subs_cell_to_write <- subs_cells %>% rename_with(correspondance_ensembl, .cols = -sample_name, ensembl_cor) %>%
  dplyr::select(dplyr::starts_with("ENS"), sample_name) %>%
  dplyr::relocate(sample_name)

write_csv(subs_cell_to_write, "/home/paulet/data/sc_mousebrain_allen/sc_allen/subsampled_cells_indiv_w_ensembl.csv")
