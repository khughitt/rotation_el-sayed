#!/bin/bash
#
# Keith Hughitt <khughitt@umd.edu>
# 2012/11/20
#
# Parses combined DEseq count table and generates two separate tables:
# one containing a single replicate for each time period and another containing
# a single control for each time period.
#
# Some additional processing is applied to clean up the table and make it
# easier to work with in Python.
#

#
# Replicates used:
#
# 64 - 4 hours
# 66 - 6 hours
# 68 - 12 hours
# 54 - 20 hours
# 70 - 24 hours
# 56 - 48 hours
# 58 - 72 hours
#
# TODO: Add controls as they become available.
#
inputdir=../../data/DEseq_comparisons/hpgl0062vs64vs66vs68vs54vs59vs70vs56vs60vs72vs58vs61_tophatv2.0.3_deseq_esmer_oneloci
inputfile=hpgl0062vs64vs66vs68vs54vs59vs70vs56vs60vs72vs58vs61_tophatv2.0.3_deseq_esmer_oneloci.counttable.sorted
output=../../data/tcruzi_counts

# REPLICATES
# Grab columns of interest and rename col headers to include only hour portion
# remove 2nd, 3rd, and last two rows (non-genes)
awk '{print $1, $3, $4, $5, $6, $8, $9, $12}' \
    "${inputdir}/${inputfile}" |
sed '1 s/hr//g' |
sed '2,3d' | 
sed 's/[ ]\+/,/g'| 
sed 's/,$//g' | 
head --lines=-3 > ${output}_replicates

