library(STdeconvolve)
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

cd <- visium@assays$spatial@counts
pos <- GetTissueCoordinates(visium)
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 15, by = 1), ncores = 16)

# les deux méthodes de choix font des choix chelous, d'un côté choississant 4 topics (premier coude)
# De l'autre choississant la perplexité minimum (donc ici 14)
# donc nous on fait le choix entre les deux, 8, le second coude. En plus c'est proche de 7
# 9 minutes pour faire 7 modèles
# mais bon j'ai fait retourner avec by = 1, et 16 coeurs, et ça prend 11 minutes.
# Résultat je reste sur 8

my_best_model <- ldas$models$`8`
results <- getBetaTheta(my_best_model, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

# representation interne
plt <- vizAllTopics(theta = deconProp,
                   pos = pos,
                   r = 3,
                   lwd = 0,
                   showLegend = TRUE,
                   plotTitle = NA) +
  guides(fill=ggplot2::guide_legend(ncol=2)) +
  ## outer border
  geom_rect(data = data.frame(pos),
            ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                         ymin = min(y)-90, ymax = max(y)+90),
            fill = NA, color = "black", linetype = "solid", size = 0.5) +

  theme(
    plot.background = ggplot2::element_blank()
  ) +
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  guides(colour = "none")


# représentation via seurat
visium <- AddMetaData(visium, metadata = deconProp %>% as.data.frame())
sp_plot <- SpatialFeaturePlot(visium, features = colnames(deconProp))
