#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import random
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from Neoantigene_modules import (
    convertEmptyGene,
    createSNVMutPeptide,
    summarizeMutWTpeptideDf,
)

neo_antigen_module = "/Volumes/Research/GitHub/NeoantigenIdentify/Module"
# r'C:\Users\abc73\Documents\GitHub\NeoantigenIdentify'
# '/Volumes/Research/GitHub/NeoantigenIdentify'
sys.path.append(neo_antigen_module)


### This function uses somatic mutation files to create to create two tables.
### the mutation results must contain "gene_name", "ensembl_transcripts", "total_protein_change", "sample_name" columns to create the tables, so it can come from pipeline results of TOSMIC or ml results of TOSMIC
## we can also extract certain genes only. ex:["TP53", "BRAF","PIK3CA","KRAS","NRAS"]
# Table A: table that contains mutant peptides column and wild-type peptides column and those peptide mutation information
# Table B: table containing  "Peptide","Amino_acids","Gene","Ensembl_transcript","Sample_name" table that can be used for peptide prediction with flurry or NetMhcpan
# necessary files: dog_gene_transcript_df_3.1.99.txt and mutation results


def create_mut_pep_sum(
    mut_file,
    data_source_folder,
    target_gene="",
    extract_data="all",
    includeRandom=False,
    nmer_list=list(range(8, 11)),
):
    np.random.seed(0)
    final_return = {}
    sequence_file = data_source_folder / "dog_gene_transcript_df_3.1.99.txt"

    include_gene_list = target_gene.split(",")

    seq_source = pd.read_csv(sequence_file, sep="\t")
    # use Ensembl_transcript as the index and transform to dict to increase the searching speed
    seq_source.index = seq_source["Ensembl_transcript"]
    seq_source_dict = seq_source.to_dict()
    mut_data = pd.read_csv(mut_file, sep="\t")
    ## in case the case of the column name is different, we change all to lowercase
    mut_data_col = [i.lower() for i in list(mut_data.columns)]
    mut_data.columns = mut_data_col

    ### check if the table contains all the necessary columns
    ## need gene_name, ensembl_transcript, total_protein_change, sample_name
    require_col = [
        "gene_name",
        "ensembl_transcripts",
        "total_protein_change",
        "sample_name",
    ]
    if len(intersection(require_col, mut_data_col)) != 4:
        miss_col = " ".join(set(require_col) - set(mut_data_col))
        print("missing column " + miss_col)
    else:
        ### consider whether we want somatic from pipeline only, model prediction, or use all data
        if extract_data == "p":
            mut_source = mut_data.loc[
                mut_data.status.str.contains("Somatic", case=False)
            ]
        elif extract_data == "m":
            mut_source = mut_data.loc[
                mut_data.model_prediction.str.contains("Somatic", case=False)
            ]
        elif extract_data == "all":  ## use all of the data
            mut_source = mut_data

        ## only grep SNV to create mutant peptides
        mut_source = mut_source.loc[
            mut_source.consequence.str.contains("nonsynonymous SNV")
        ]
        if include_gene_list != [""]:
            mut_source = mut_source[mut_source.gene_name.isin(include_gene_list)]
        mut_source = mut_source.to_dict("records")

        total_sum = []
        for each_nmer in nmer_list:
            for i, j in enumerate(mut_source):
                gene_name = mut_source[i]["gene_name"]
                ensembl_transcript = mut_source[i]["ensembl_transcripts"]
                total_protein_change = mut_source[i]["total_protein_change"]
                sample_name = mut_source[i]["sample_name"]
                # print(total_protein_change)

                ## need to use re to extract mutation protein
                # x = re.search(r'(p.[A-Za-z])(\d+)([A-Za-z])', total_protein_change) ## if pure Annovar, use this one
                x = re.search(r"([A-Za-z])(\d+)([A-Za-z])", total_protein_change)
                wildtype = x.group(1)
                mut_amino = x.group(3)  #'R'
                mut_position = int(x.group(2))  # int(1047)
                # check if the transcript has protein annotation and if it is a novel transcript (no peptide annotation), it will skip
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

        ## check if mutant peptides and wildtype peptides have the same length
        mut_wt_sum_df = summarizeMutWTpeptideDf(final_sum)

        if includeRandom:
            random_peptide_file = data_source_folder / "Pan_cancer_random_peptide.txt"
            random_peptide = pd.read_csv(random_peptide_file, sep=",")
            final_sum = pd.concat([final_sum, random_peptide])

        final_return["summary"] = mut_wt_sum_df
        final_return["meta"] = final_sum

    return final_return
