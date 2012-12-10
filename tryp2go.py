#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Tryp2GO
Keith Hughitt <khughitt@umd.edu>

2012/12/09

Processes a TriTrypDB file containing Trypanosome gene ontology (GO) terms and
return a new table mapping gene id's to GO terms.

@TODO
-Handle annotated/predicted GO terms separately?
"""
import sys
import re
from goatools import obo_parser

# Gene ontology (1.2)
go_obo = "../data/gene_ontology.1_2.obo"

# Create a dict to map from GO names to identifiers
go_dag = obo_parser.GODag(go_obo)
go_terms = dict((v.name, k) for k, v in go_dag.items())

#TODO handle annotated vs. predicted differently?
input_file = "../data/TriTrypDB-4.2_TcruziEsmeraldo-LikeGene.txt"

# Parse TriTrypDB GO terms
current_id = None
mapping = []

regex = '(\w*) GO (\w*): ([\w, \ \/\[\]\-]*)'

for line in open(input_file).readlines():
    if line.startswith("Gene ID"):
        current_id = line.split(": ").pop().strip()

    elif line.startswith('Annotated GO'):
        # TYPE, CATEGORY, GO TERM NAMES
        groups = re.match(regex, line).groups()
        
        # make sure correct number of terms matched
        if len(groups) < 3:
            print("short line found: %s" % line)
            continue

        go_names = groups[2].strip()
        no_match = False

        # skip null entries (not helpful)
        if go_names == 'null':
            continue

        # skip very general (top-level) matches
        elif go_names in ['biological_process', 'cellular_component', 'molecular_function']:
            continue

        # first try to add entire line in case there are commas in the term itself
        if go_names in go_terms:
            mapping.append([current_id, go_terms[go_names], 
                            go_names, groups[0], groups[1]])
            continue

        # next, split into split comma-separated terms
        names = [x.strip() for x in go_names.split(',')]

        for name in names:
            if name in go_terms:
                mapping.append([current_id, go_terms[name], 
                                name, groups[0], groups[1]])
            else:
                no_match = True
        
        if no_match:
            print("No match found for: %s" % (go_names))
        

