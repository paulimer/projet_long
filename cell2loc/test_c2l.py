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
import squidpy as sq

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

st_data_folder = "/home/paulet/data/sc_mousebrain_allen/st_mousebrain/"
sc_data_folder = "/home/paulet/data/sc_mousebrain_allen/sc_allen/"
results_data_folder = "/home/paulet/Documents/TAAAAF/Projet_long/cell2loc/results"
ref_run_name = f'{results_data_folder}/reference_signatures'
run_name = f'{results_data_folder}/cell2location_map'

adata_vis = sq.read.visium(path="/home/paulet/data/sc_mousebrain_allen/st_mousebrain/",
                           counts_file="V1_Mouse_Brain_Sagittal_Anterior_Section_2_filtered_feature_bc_matrix.h5",
                           source_image_path="/home/paulet/data/sc_mousebrain_allen/st_mousebrain/V1_Mouse_Brain_Sagittal_Anterior_Section_2_image.tif")

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
ref_metadata = ref_metadata.set_index("sample_name")
subs_ref_data = ref_metadata.loc[adata_ref.obs_names, ]
adata_ref.obs["neighborhood_label"] = pd.Categorical(subs_ref_data["neighborhood_label"])
adata_ref.obs["donor_name"] = pd.Categorical(subs_ref_data["external_donor_name_label"])


from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()


# prepare ann data for the regression model need to integrate metadata.csv at some point
# TODO : need to integrate metadata.csv at some point
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        # batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key="neighborhood_label",
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['donor_name']
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

#maybe add more cells
# mod.plot_QC()

# extracxt reference cell types as a pd.DataFrame
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, :]




# SPATIAL MAPPING
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )
# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.show()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

# mod.plot_QC()

fig = mod.plot_spatial_QC_across_batches()

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
# from cell2location.utils import select_slide
# slide = select_slide(adata_vis)

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sq.pl.spatial_scatter(adata_vis, img_cmap='magma',
                  # show first 8 cell types
                  color=['CGE', 'DG/SUB/CA', 'L2/3 IT', 'L4/5/6 IT Car3', 'MGE', 'NP/CT/L6b', 'Other', 'PT'],
                  ncols=4, size=1.3)#,
                  # img_key='hires')#,
                  # limit color scale at 99.2% quantile of cell abundance
                  # vmin=0, vmax='p99.2'
                 # )
