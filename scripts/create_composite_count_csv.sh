#!/bin/bash
#
# Generates a composite count table CSV file containing one replicate from each
# available time period.
#
# Keith Hughitt <khughitt@umd.edu>
# 2012/11/12
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
# Controls used:
#
# 63 - 4 hours
# 65 - 6 hours
# 67 - 12 hours
# 53 - 20 hours
# 69 - 24 hours
# 55 - 48 hours
# 57 - 72 hours
#
basedir=../data/DEseq_comparisons
input=hpgl63vs64vs65vs66vs67vs68vs53vs54vs59vs69vs70vs55vs71vs56vs60vs72vs57vs58vs61_hg19_tophatv2.0.3_accepted_hits.counttable
output=../data/combined_countable

# REPLICATES
# Grab columns of interest and rename col headers to include only hour portion
# remove 2nd, 3rd, and last two rows (non-genes)
awk '{print $1, $3, $5, $7, $9, $12, $15, $19}' \
    "${basedir}/hpgl0053-0072_tophatv2.0.3_deseq_HFF_loose/${input}" |
sed '1 s/hrs_after_tc_infection1//g' |
sed '2,3d' | 
sed 's/ /,/g' |
head --lines=-3 > ${output}_replicates

# CONTROLS
# Grab columns of interest and rename col headers to include only hour portion
# remove 2nd, 3rd, and last two rows (non-genes)
awk '{print $1, $2, $4, $6, $8, $11, $13, $18}' \
    "${basedir}/hpgl0053-0072_tophatv2.0.3_deseq_HFF_loose/${input}" |
sed '1 s/hrs_before_tc_infection1//g' |
sed '1 s/hrs_before_tc_infection//g' |
sed '2,3d' | 
sed 's/ /,/g' |
head --lines=-3 > ${output}_controls

