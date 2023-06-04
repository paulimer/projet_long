#!/usr/bin/env sh

#########################################################################################################
################## Stereoscope, tool to deconvolute signal from spatial transcriptomic ##################
########################### github repo : https://github.com/almaan/stereoscope #########################
#########################################################################################################


########################################################################################
### Prediction of NAd4 spots with SC Normal Adrenal integrated (dataset.combined) ######
########################################################################################


stereoscope run --sc_cnt ~/data/sc_mousebrain_allen/sc_allen/subsampled_cells_indiv.tsv \
--sc_labels ~/data/sc_mousebrain_allen/sc_allen/metadata_stsc.tsv \
-sce 75000 \
-o ~/Documents/TAAAAF/Projet_long/stereoscope/results \
-n 5000 \
--st_cnt ~/data/sc_mousebrain_allen/st_mousebrain/st_stereo_rn.tsv \
-ste 75000 \
-stb 5000 \
-scb 5000 \
--gpu
