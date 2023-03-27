#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import pandas as pd

cells_metadata = pd.read_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/metadata.csv")

plt.bar(cells_metadata.groupby("neighborhood_label").size().index, cells_metadata.groupby("neighborhood_label").size())
plt.show()

gb_neigh = cells_metadata.groupby("neighborhood_label")
neigh_cells = [gb_neigh.get_group(x) for x in gb_neigh.groups]

subsm_df = pd.DataFrame()
for df in neigh_cells:
    subsm_df = pd.concat([subsm_df, df.sample(1000, axis=0)], axis=0)

subsm_sample_names = list(subsm_df["sample_name"])

part_csvs = ["/home/paulet/data/sc_mousebrain_allen/sc_allen/" + d for d in os.listdir("/home/paulet/data/sc_mousebrain_allen/sc_allen/", ) if d.startswith("gene_expression_matrix.00")]


cells_df = pd.DataFrame()
for csv in part_csvs:
    tmp_df = pd.read_csv(csv)
    cells_df = pd.concat([cells_df, tmp_df[tmp_df["sample_name"].isin(subsm_sample_names)]], axis=0)
    print("done concat")
    del tmp_df
    print("done del")

cells_df.to_csv("/home/paulet/data/sc_mousebrain_allen/sc_allen/subsampled_cells.csv", index=False)
