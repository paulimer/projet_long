library(spacexr)
library(Matrix)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(janitor)

visium <- Load10X_Spatial("/home/paulet/data/sc_mousebrain_allen/st_mousebrain",
                          filename = "V1_Mouse_Brain_Sagittal_Anterior_Section_2_filtered_feature_bc_matrix.h5",
                          assay = "spatial",
                          slice = "V1_Mouse_Brain_Sagittal_Anterior_Section_2_image.tif")
plot1 <- VlnPlot(visium, features = "nCount_spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(visium, features = "nCount_spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

ref_counts <- read_csv("~/data/sc_mousebrain_allen/sc_allen/subsampled_cells_indiv.csv")

metadata <- read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/metadata.csv") %>%
  select(sample_name, neighborhood_label) %>%
  mutate(neighborhood_label = str_replace_all(neighborhood_label, "/", "_"))
metadata <- metadata[metadata$sample_name %in% ref_counts$sample_name, ]
metadata <- metadata %>% column_to_rownames("sample_name")
ref_counts <- ref_counts %>%
  column_to_rownames("sample_name") %>%
  t() %>%
  as.matrix()

ref <- CreateSeuratObject(
  counts = ref_counts,
  project = "ref",
  meta.data = metadata
)

ref <- UpdateSeuratObject(ref)
Idents(ref) <- "neighborhood_label"

# extract information to pass to the RCTD Reference function
counts <- ref@assays$RNA@counts
cluster <- as.factor(ref$neighborhood_label)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

counts_vis <- visium@assays$spatial@counts
coords <- GetTissueCoordinates(visium)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts_vis, colSums(counts_vis))

RCTD <- create.RCTD(query, reference, max_cores = 16)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
saveRDS(RCTD, file = "rctd_obj_full.rds")
visium <- AddMetaData(visium, metadata = RCTD@results$weights)
saveRDS(visium, file = "visium_w_rctd_full.rds")


# Load exploit results ------------------
rctd_res <- readRDS("./rctd/rctd_obj_full.rds")
weights_rctd <- rctd_res@results$weights %>% as.data.frame() %>% clean_names
weights_rctd_per <- t(apply(weights_rctd, 1, function(x){x/sum(x)}))



visium <- AddMetaData(visium, metadata = weights_rctd_per %>% as.data.frame())
sp_plot <- SpatialFeaturePlot(visium, features = colnames(weights_rctd_per))

# change to plot all types? how?
ggsave("rctd_full.png")
