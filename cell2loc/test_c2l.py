#!/usr/bin/env python3

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
# if using the CPU uncomment this:
#os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,openmp=True,force_device=True'

import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

st_data_folder = "/home/paulet/data/sc_mousebrain_allen/st_mousebrain/"
sc_data_folder = "/home/paulet/data/sc_mousebrain_allen/sc_allen/"
results_data_folder = "/home/paulet/Documents/TAAAAF/Projet_long/cell2loc/results"
ref_run_name = f'{results_data_folder}/reference_signatures'
run_name = f'{results_data_folder}/cell2location_map'

adata_vis = sc.read_visium(path="/home/paulet/data/sc_mousebrain_allen/st_mousebrain/",
                           count_file="V1_Mouse_Brain_Sagittal_Anterior_Section_2_filtered_feature_bc_matrix.h5",
                           source_image_path="V1_Mouse_Brain_Sagittal_Anterior_Section_2_image.tif")

adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# find mt genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


# import ref data - évidemment pas en ENSEMBL putain comment on fait en python -> on fait en R
adata_ref = sc.read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/subsampled_cells_wo_index_w_ensembl.csv", first_column_names=True)
ref_metadata = pd.read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/metadata.csv")
neigh_data = ref_metadata.filter(["sample_name", "neighborhood_label"])
subs_cells = pd.read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/subsampled_cells_wo_index_w_ensembl.csv").iloc[:, 0]
subs_neigh_data = neigh_data[neigh_data["sample_name"].isin(subs_cells)]
adata_ref.obs["cell_neigh"] = pd.Categorical(subs_neigh_data["neighborhood_label"])


from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()


# prepare ann data for the regression model need to integrate metadata.csv at some point
# TODO : need to integrate metadata.csv at some point
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        # batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        # labels_key='Subset',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=['Method']
                       )


from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.view_anndata_setup()

# You are using a CUDA device ('NVIDIA GeForce RTX 3070') that has Tensor Cores. To properly utilize them, you should set `torch.set_float32_matmul_precision('medium' | 'high')` which will trade-off precision for performance. For more details, read https://pytorch.org/docs/stable/generated/torch.set_float32_matmul_precision.html#torch.set_float32_matmul_precision
from torch import set_float32_matmul_precision
set_float32_matmul_precision('high')
mod.train(max_epochs=250, use_gpu=True)


mod.plot_history(20)
plt.show()
# Ça a l'air d'aller

# Export the summary of the posterior distrib
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

mod.plot_QC()

# extracxt reference cell types as a pd.DataFrame
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]
