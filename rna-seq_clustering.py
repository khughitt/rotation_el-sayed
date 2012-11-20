#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
DESeq clustering
Keith Hughitt <khughitt@umd.edu>
November 20, 2012

Reads in DESeq input count table and attempts to cluster
genes which show similar trends of expression.    
"""
import numpy as np
from pandas import read_csv, Series
from pandas.tools.plotting import parallel_coordinates
from matplotlib import pyplot as plt
from scipy.cluster.vq import kmeans, whiten, vq

"""Main"""
#target="human"
target="tcruzi"

# Read in data
col_names = ['gene', 4, 6, 12, 20, 24, 48, 72]

replicates = read_csv('../data/%s_counts_replicates' % target,
                             index_col=0, skiprows=1, header=None,
                             names=col_names)

# Currently we don't have any control DEseq counts for T. cruzi so just
# use zero-filled dataframe.
if target == "tcruzi":
    controls = replicates.copy().clip(0, 0)
else:
    controls = read_csv('../data/%s_counts_controls' % target,
                               index_col=0, skiprows=1, header=None, 
                               names=col_names)
    controls = controls / controls.sum().astype('float')

# Normalize count data
#
# For now, a simple approach to normalizing data both between sample types
# and time periods is to divide each column by the total number of counts
# that column, resulting in a proportion of total expression.
replicates = replicates / replicates.sum().astype('float')

# Subtract the control counts from the infected counts to see look at the
# differences between infected (replicates) and control data sets.
data = replicates - controls

# filter out rows which don't show any significant change in transcription
def col_mask(data, col):
    return ((data[col] < data[col].mean() - (3 * data[col].std())) |
            (data[col] > data[col].mean() + (3 * data[col].std())))

# create a boolean mask with all rows that have no significant values
# set to False
bool_mask = np.any([col_mask(data, i) for i in data.columns], axis=0)
data = data[bool_mask]

# k-means clustering
# http://glowingpython.blogspot.com/2012/04/k-means-clustering-with-scipy.html
num_clusters = 4
data_matrix = data.as_matrix()
data_matrix = whiten(data_matrix)
centroids, distortion = kmeans(data_matrix, num_clusters)
idx, distortion = vq(data_matrix, centroids)

# add column with cluster grouping at end of data
data = data.join(Series(idx, index=data.index, name='cluster'))

# create parallel coordinates plot
plt.figure()
fig = parallel_coordinates(data, 'cluster', use_columns=True,
                           colors=('#556270', '#4ECDC4', '#C7F464', '#FF6B6B'))
fig.set_title("Gene Expression vs. Time")
fig.set_xlabel("Time (hours)")
fig.set_ylabel("Relative Expression")

plt.show()

