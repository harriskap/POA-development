#!/usr/bin/env python

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mp
import math
import seaborn as sns
import os

# less common imports
from sklearn import preprocessing
from scipy.stats import f
from sklearn.cluster import AgglomerativeClustering
import matplotlib.cm as cm
sc.settings.verbosity = 0
sc.logging.print_header()

sns.set_context("paper")
import new_misc_code as nmc

n_splines = 12
n_grid_pts = 100
good_fits = nmc.load_obj("data/dev_deg_v7/AGE_gam_fits_12_grid100.filtered1000.pkl")
gam_fits = good_fits
logTMMs = nmc.load_obj(f"data/dev_deg_v7/logTMMs.pkl")
all_goods_df = pd.concat(list( good_fits.values()), keys=list(good_fits.keys()))
all_goods_df.index.names = ['cell_type', 'gene']
cmap = sns.color_palette( "ch:start=.2,rot=-.3", n_colors=1_000)
sns.set_context( context='paper')
sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'lightgrey'})

key_itr = "i-B1"
#key_itr = "e-A1"
sig_df = logTMMs[key_itr]
ax_age = [s.split("..")[1].split("_")[0] for s in sig_df.columns]
#ax_age
ai_age = []
# only for samples which are part of pseudobulk
sample_age = {}
for age in ax_age:
    i_age = 0
    if age[0] == "e":
        i_age = int(age[1:3])
    elif age[0:2] == "p0":
        i_age = 20
    elif age[0:2] == "p4":
        i_age = 24
    else:
        i_age = int(age[1:3]) + 20
    ai_age.append(i_age)
    sample_age[age] = i_age

stage_order = list(np.unique(ai_age))
print(stage_order)
#[16, 18, 20, 24, 30, 38, 48, 85]

ax = sns.clustermap(all_goods_df,
                    metric='euclidean',
                    method='ward',
                    col_cluster=False,
                    row_cluster=True,
                    cmap=cmap,
                    figsize=(15,15),
                    standard_scale=0)

nmc.save_obj(ax, "data/dev_deg_v7/ax.pkl")
