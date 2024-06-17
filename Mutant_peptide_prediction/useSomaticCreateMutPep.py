#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

# Somatic mutations must contain Ensembl transcripts and protein mutation information.

# This script takes somatic mutation files in annovar format and peptide sequences for each transcript (obtained from peptide files downloaded from Ensembl) as input. It generates n-mer mutant peptides using the provided peptide sequences.

# The 'n-mer' length for mutant peptides can be customized based on user requirements.

# The script generates a final table containing mutant peptides, wild-type peptides, and additional information derived from somatic mutations.

# It creates an input file in the 'flurry' format for MHC_flurry 2.0 and a peptide file for NetMhcpan, which can be used for predicting peptide-MHC binding affinity.



"""
import os
import re
import sys

import numpy as np
import pandas as pd

neo_antigen_module = (
    '/work/szlab/Kun-Lin_pipeline/Allele_specific_model/Model_creation'
)

sys.path.append(neo_antigen_module)
from Create_mut_peptide_modules import *

## you can specifiy the peptide length (here we create 8mer to 11 mer)
nmer = list(range(8,15 ))

mut_file = sys.argv[1]

## you can choose some target genes or only you want to include all of the mutate genes
target_gene = ""

output_base = sys.argv[2]

### necessary files
data_source_folder = Path(
    '/work/szlab/Kun-Lin_pipeline/Allele_specific_model/mut_peptide_source'
)
## include random peptide
random_peptide_file = data_source_folder / "Pan_cancer_random_peptide.txt"
sequence_file = data_source_folder / "dog_gene_transcript_df_3.1.99.txt"
mhc_flurry_file = sys.argv[3]

final_meta_table = output_base + "/" + "Mut_peptide_meta_table.txt"
final_mut_peptide = output_base + "/" + "Mut_peptide.txt"

include_gene_list = target_gene.split(",")
# ["TP53", "BRAF","PIK3CA","KRAS","NRAS"]
random_peptide = pd.read_csv(random_peptide_file, sep=",")
seq_source = pd.read_csv(sequence_file, sep="\t")
# use Ensembl_transcript as the index and transform to dict to increase the searching peed
seq_source.index = seq_source["Ensembl_transcript"]
seq_source_dict = seq_source.to_dict()
mut_source = pd.read_csv(mut_file, sep="\t")
mut_source_col = list(mut_source.columns)
## if the data derved from TOSMIC, then only grep somatic mutations from Model_prediction
if "Model_prediction" in mut_source_col:
    mut_source = mut_source.loc[
        mut_source.Model_prediction.str.contains("Somatic", case=False)
    ]
## if the data derved from mutation prediction pipeline, then only grep somatic mutations from pipeline filtering
elif "Status" in mut_source_col:
    mut_source = mut_source.loc[mut_source.Status.str.contains("Somatic", case=False)]
else:
    mut_source = mut_source

## only grep SNV to create mutant peptides
mut_source = mut_source.loc[mut_source.Consequence == "nonsynonymous SNV"]
if include_gene_list != [""]:
    mut_source = mut_source[mut_source.Gene_name.isin(include_gene_list)]
mut_source = mut_source.to_dict("records")

total_sum = []
for each_nmer in nmer:
    for i, j in enumerate(mut_source):
        gene_name = mut_source[i]["Gene_name"]
        ensembl_transcript = mut_source[i]["Ensembl_transcripts"]
        Total_protein_change = mut_source[i]["Total_protein_change"]
        sample_name = mut_source[i]["Sample_name"]
        ## need to use re to extract mutation protein
        # x = re.search(r'(p.[A-Za-z])(\d+)([A-Za-z])',Total_protein_change) ## if pure annovar, use this one
        x = re.search(r"([A-Za-z])(\d+)([A-Za-z])", Total_protein_change)
        wildtype = x.group(1)
        mut_amino = x.group(3)  #'R'
        mut_position = int(x.group(2))  # int(1047)
        # check if the transcript have protein annotation and if it is a novel transcript (no peptide annotation), it will skip
        if ensembl_transcript in seq_source_dict["Sequence"].keys():
            protein_seq = seq_source_dict["Sequence"][ensembl_transcript]
            if protein_seq[int(mut_position) - 1] == wildtype:
                total_wild_mut_peptide = createSNVMutPeptide(
                    protein_seq,
                    each_nmer,
                    mut_position,
                    mut_amino,
                    gene_name,
                    ensembl_transcript,
                    sample_name,
                )
                total_sum += total_wild_mut_peptide
            else:
                print("Wildtype aa is not the same as the data shown")
final_sum = pd.DataFrame(total_sum)
final_sum.columns = [
    "Peptide",
    "Amino_acids",
    "Gene",
    "Ensembl_transcript",
    "Sample_name",
]
final_sum["Gene"] = final_sum["Gene"].apply(convertEmptyGene)

### We can also have a summarize table and check if mut peptides and wildtype peptides have same length
# mut_wt_sum_df = summarizeMutWTpeptideDf(final_sum)
# mut_wt_sum_df.to_csv(output_base + '/' + 'mut_wt_peptide_summary',
#                   sep ="\t", index = False ,header=True)


final_sum = pd.concat([final_sum, random_peptide])
final_sum.to_csv(final_meta_table, sep="\t", index=False, header=True)

final_sum.Peptide.to_csv(final_mut_peptide, sep="\n", index=False, header=None)


mhc_flurry_support_allele = pd.read_csv(mhc_flurry_file, header=None)
mhc_flurry_list = mhc_flurry_support_allele[
    mhc_flurry_support_allele[0].str.contains("DLA")
][0].to_list()
mhc_flurry_list = set(mhc_flurry_list)
if os.path.isdir(output_base + "/" + "Mhc_flurry_input")==False:
    os.mkdir(output_base + "/" + "Mhc_flurry_input")
    
for allele in mhc_flurry_list:
    #allele = allele.replace("*", "_").replace(":", "_")
    flurry_output = open(output_base + "/" + "Mhc_flurry_input" + "/" + allele + "_flurry_mut_peptide.txt", "w")
    flurry_output.write("allele,peptide\n")
    for peptide in final_sum.Peptide.to_list():
        flurry_output.write(allele + "," + peptide + "\n")

    flurry_output.close()
