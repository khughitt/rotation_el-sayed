#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Tryp2GO
Keith Hughitt <khughitt@umd.edu>

2012/12/09

Processes a TriTrypDB file containing Trypanosome gene ontology (GO) terms and
return a new table mapping gene id's to GO terms.
"""
import sys

input_file = "../data/TriTrypDB-4.2_TcruziEsmeraldo-LikeGene.txt"

# Parse TriTrypDB GO terms
current_id = None
mapping = []

for line in open(input_file).readlines():
    if line.startswith("Gene ID"):
        current_id = line.split(": ").pop().strip()
    elif line.startswith("GO:"):
        mapping.append([current_id] + line.split('\t')[0:5])

