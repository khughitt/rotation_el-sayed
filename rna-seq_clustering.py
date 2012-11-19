#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
DESeq clustering
Keith Hughitt <khughitt@umd.edu>
November 12, 2012

Reads in DESeq input count table and attempts to cluster
genes which show similar trends of expression.

Example usage:

    count_table_clustering.py count_table
    
"""
import pandas
import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.vq import kmeans, whiten, vq

# Read in data
#replicates = pandas.read_csv('../data/combined_countable_replicates', index_col=0)
#controls = pandas.read_csv('../data/combined_countable_controls', index_col=0)
replicates = pandas.read_csv('../data/combined_countable_replicates',
                             index_col=0, skiprows=1, header=None,
                             names=['gene', 4, 6, 12, 20, 24, 48, 72])
controls = pandas.read_csv('../data/combined_countable_controls',
                           index_col=0, skiprows=1, header=None,
                           names=['gene', 4, 6, 12, 20, 24, 48, 72])

# Normalize count data
#
# @TODO talk to Yuan/check logs to see if any normalization has been done
# prior to count table generation...
#
# For now, a simple approach to normalizing data both between sample types
# and time periods is to divide each column by the total number of counts
# that column, resulting in a proportion of total expression.
replicates = replicates / replicates.sum().astype('float')
controls = controls / controls.sum().astype('float')

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
data_matrix = data.as_matrix()
data_matrix = whiten(data_matrix)
centroids, distortion = kmeans(data_matrix, 4)
idx, distortion = vq(data_matrix, centroids)

# add column with cluster grouping at end of data
#x = np.column_stack([x, idx])
data = data.join(pandas.Series(idx, index=data.index, name='cluster'))

plt.figure()
pandas.tools.plotting.parallel_coordinates(data, 'cluster', 
                                           colors=('#556270', '#4ECDC4', '#C7F464', '#FF6B6B'),
                                           use_columns=True)
#plt.legend().set_visible(False)
plt.show()


