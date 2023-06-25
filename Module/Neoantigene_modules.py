#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 11:30:48 2021

The script contains the function that regular use for Neoantigen identification for Flurry or netMhcPan

@author: kun-linho
"""

import collections
import os
import random
import re
import sys
from copy import copy
from math import log, log2, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import distance
from six import StringIO

## the function that use to extract the ensembl gene and ensembl transcripts given the annovar annotated results
## ex: 'ENSCAFG00000023363:ENSCAFT00000000040:exon7:c.762dupA:p.G255fs,' and it will return ENSCAFG0000023363, ENSCAFT0000000040, G255fs
## even we have multiple annotation, it will join them together with ','


def extractAnnovarMutProtein(mut_info):
    try:
        Ensembl_gene = re.findall(r"(ENSCAFG)(\d+):", mut_info)
        Ensembl_trans = re.findall(r"(ENSCAFT)(\d+):", mut_info)

        ## the overall regex for annovar
        total_protein_mut = re.findall(r"p.([A-Z0-9a-z_.*-]*),", mut_info)
        total_ensembl_gene = ["".join(i) for i in Ensembl_gene]
        Total_protein_change = ["".join(i) for i in total_protein_mut]
        total_trans = ["".join(i) for i in Ensembl_trans]
        diff = abs(len(Total_protein_change) != len(total_trans))
        # in case that the annovar ensemble transcripts and the protein changes are not the same length, I add "No_Info_Provided" until they are the same
        while diff != 0:
            if len(Total_protein_change) > len(total_trans):
                total_trans += ["No_Info_Provided"]
                total_ensembl_gene += ["No_Info_Provided"]
            elif len(Total_protein_change) < len(total_trans):
                Total_protein_change += ["No_Info_Provided"]

            diff = abs(len(Total_protein_change) != len(total_trans))

        # combine into a big string
        Total_ensembl_gene = ",".join(total_ensembl_gene)
        Total_protein_change = ",".join(Total_protein_change)
        total_trans = ",".join(total_trans)

        final_return = pd.Series(
            [Total_ensembl_gene, total_trans, Total_protein_change]
        )

        return final_return
    except:
        print("There is an error in sample " + sample_name)


## this function will process annovar output and extract all the protein changes in the annovar file (not only one protein change)
## it will remove Ensembl gene that located in retro_gene list idenitfied in pan-cancer paper
## even the annovar contains sample name, it will overwrite the orignal one
## the final return table contains the following
## Line','Consequence','Gene_name','Chrom','Start','End','Ref','Alt','Sample_name','Ensembl_gene','Ensembl_transcripts','Total_protein_change'
## the reason I keep line is for future used for VAF calculation
def processAnnovar(annovar_gene_file, retro_gene_file, sample_name):
    retro_gene_list = pd.read_csv(retro_gene_file, sep="\n", header=None)
    retro_gene_list = retro_gene_list[0].to_list()
    annovar_gene_data = pd.read_csv(annovar_gene_file, sep="\t", header=None)
    annovar_output_col = [
        "Line",
        "Consequence",
        "Annovar_info",
        "Chrom",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Homo_hetro",
        "10",
        "11",
        "12",
        "13",
        "Gene_name",
        "Sample_name",
    ]
    annovar_gene_data.columns = annovar_output_col
    # annovar_gene_data = annovar_gene_data[annovar_gene_data[1]!='synonymous SNV']
    # annovar_gene_data.loc[:,'Sample_name']= sample_name
    # rename the first column for the future merge
    annovar_gene_data.loc[
        :, ["Ensembl_gene", "Ensembl_transcripts", "Total_protein_change"]
    ] = (annovar_gene_data["Annovar_info"].apply(extractAnnovarMutProtein).to_numpy())
    ## in case the Ensembl_transcripts and Total_protein_change are not in the same length, I add "No_Info_Provided" in the extractAnnovarMutProtein steps, so we need to exclude that information
    annovar_gene_data = annovar_gene_data.loc[
        (annovar_gene_data["Ensembl_transcripts"] != "No_Info_Provided")
        & (annovar_gene_data["Total_protein_change"] != "No_Info_Provided")
    ]
    target_annovar_info = annovar_gene_data.loc[
        :,
        [
            "Line",
            "Consequence",
            "Gene_name",
            "Chrom",
            "Start",
            "End",
            "Ref",
            "Alt",
            "Sample_name",
            "Ensembl_gene",
            "Ensembl_transcripts",
            "Total_protein_change",
        ],
    ]
    # target_annovar_info = annovar_gene_data[['Line',1,2,14,'Ensembl_gene','Ensembl_transcripts','Total_protein_change']]
    # target_annovar_info.columns = ['Line','Consequence','Gene_name','Sample_name','Ensembl_gene','Ensembl_transcripts','Total_protein_change']
    target_annovar_info.loc[:, "Gene_name"] = (
        target_annovar_info["Gene_name"].astype(str).apply(lambda x: x.split(","))
    )
    target_annovar_info.loc[:, "Ensembl_gene"] = (
        target_annovar_info["Ensembl_gene"].astype(str).apply(lambda x: x.split(","))
    )
    target_annovar_info.loc[:, "Ensembl_transcripts"] = (
        target_annovar_info["Ensembl_transcripts"]
        .astype(str)
        .apply(lambda x: x.split(","))
    )
    target_annovar_info.loc[:, "Total_protein_change"] = (
        target_annovar_info["Total_protein_change"]
        .astype(str)
        .apply(lambda x: x.split(","))
    )
    #### unnesting is the same as explode but just in case pandas version is < 1.3
    # target_annovar_info = unnesting(target_annovar_info,['Ensembl_transcripts','Total_protein_change'])
    target_annovar_info = target_annovar_info.explode(
        ["Gene_name", "Ensembl_gene", "Ensembl_transcripts", "Total_protein_change"]
    )
    target_annovar_info = target_annovar_info.loc[
        (target_annovar_info["Ensembl_transcripts"] != "No_Info_Provided")
        & (target_annovar_info["Total_protein_change"] != "No_Info_Provided")
    ]
    ### filter retrogene list with given ensembl id
    target_annovar_info = target_annovar_info[
        ~target_annovar_info.Ensembl_gene.isin(retro_gene_list)
    ]
    target_annovar_info.loc[:, "Gene_mut_info"] = (
        target_annovar_info["Gene_name"]
        + "_"
        + target_annovar_info["Total_protein_change"]
    )
    target_annovar_info.loc[:, "Transcript_mut_info"] = (
        target_annovar_info["Ensembl_transcripts"]
        + "_"
        + target_annovar_info["Total_protein_change"]
    )
    target_annovar_info = target_annovar_info.drop_duplicates()
    return target_annovar_info


## this function is used to reduce the memory usage for the DataFrame
def reduce_memory_usage(df, verbose=True):
    numerics = ["int8", "int16", "int32", "int64", "float16", "float32", "float64"]
    start_mem = df.memory_usage().sum() / 1024**2
    for col in df.columns:
        col_type = df[col].dtypes
        if col_type in numerics:
            c_min = df[col].min()
            c_max = df[col].max()
            if str(col_type)[:3] == "int":
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)
            else:
                if (
                    c_min > np.finfo(np.float16).min
                    and c_max < np.finfo(np.float16).max
                ):
                    df[col] = df[col].astype(np.float16)
                elif (
                    c_min > np.finfo(np.float32).min
                    and c_max < np.finfo(np.float32).max
                ):
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
    end_mem = df.memory_usage().sum() / 1024**2
    if verbose:
        print(
            "Mem. usage decreased to {:.2f} Mb ({:.1f}% reduction)".format(
                end_mem, 100 * (start_mem - end_mem) / start_mem
            )
        )
    return df


### the function will take a processed annovar output file as input to create nmer mutant peptide df
### the function will filter out mutations that are not SNV (ex: stop gain, fs)
def dfCreateSNVMut(df, protein_source_file, nmer):
    ## create source_proteins dict for creating mutant peptides
    seq_source = pd.read_csv(protein_source_file, sep="\t")
    # use Ensembl_transcript as the index and transform to dict to increase the searching peed
    seq_source.index = seq_source["Ensembl_transcript"]
    seq_source_dict = seq_source.to_dict()

    snv_df = df[df["Consequence"].str.contains("SNV")]
    mut_source = snv_df.to_dict("records")

    total_sum = []
    for i, j in enumerate(mut_source):
        gene_name = mut_source[i]["Gene_name"]
        ensembl_transcript = mut_source[i]["Ensembl_transcripts"]
        protein_change = mut_source[i]["Total_protein_change"]
        sample_name = mut_source[i]["Sample_name"]
        # print(protein_change)

        ## need to use re to extract mutation protein
        x = re.search(r"([A-Za-z])(\d+)([A-Za-z])", protein_change)
        wildtype = x.group(1)
        mut_amino = x.group(3)  #'R'
        mut_position = int(x.group(2))  # int(1047)
        # check if the transcript have protein annotation and if it is a novel transcript (no peptide annotation), it will skip
        if ensembl_transcript in seq_source_dict["Sequence"].keys():
            protein_seq = seq_source_dict["Sequence"][ensembl_transcript]
            total_wild_mut_peptide = createSNVMutPeptide(
                protein_seq,
                nmer,
                mut_position,
                mut_amino,
                gene_name,
                ensembl_transcript,
                sample_name,
            )
            total_sum += total_wild_mut_peptide

    final_sum = pd.DataFrame(total_sum)
    final_sum.columns = [
        "Peptide",
        "Amino_acids",
        "Gene",
        "Ensembl_transcripts",
        "Sample_name",
    ]
    final_sum["Gene"] = final_sum["Gene"].apply(convertEmptyGene)
    return final_sum


## the function is same as df.explode but in case pandas version < 1.3
## there is a extra function
def unnesting(df, explode):
    idx = df.index.repeat(df[explode[0]].str.len())
    df1 = pd.concat(
        [pd.DataFrame({x: np.concatenate(df[x].values)}) for x in explode], axis=1
    )
    df1.index = idx

    return df1.join(df.drop(explode, 1), how="left")


def createSNVMutPeptide(
    prot_seq, nmer, mut_pos, mut_amino, Gene, ensembl_transcript, sample_name
):
    ## consider mutation position, if it happens in the very beginning or very end of the peptide sequence

    finalpeptide = []

    ## if mutation occurs in the begining
    if mut_pos < nmer:
        potential_mut_peptide = (
            prot_seq[0 : mut_pos - 1]
            + mut_amino
            + prot_seq[int(mut_pos) : mut_pos + nmer]
        )
        wildtype_peptide = prot_seq[int(0) : mut_pos + nmer]

        for i in range(0, mut_pos):
            mutpeptide = potential_mut_peptide[i : i + nmer]
            wildtype = wildtype_peptide[i : i + nmer]
            # print(mutpeptide)
            finalpeptide.append(
                (
                    mutpeptide,
                    prot_seq[mut_pos - 1] + str(mut_pos) + mut_amino,
                    Gene,
                    ensembl_transcript,
                    sample_name,
                )
            )
            finalpeptide.append(
                (
                    wildtype,
                    prot_seq[mut_pos - 1] + str(mut_pos) + mut_amino + "_WildType",
                    Gene + "_WildType",
                    ensembl_transcript,
                    sample_name,
                )
            )

    ## if the mutation occurs in the end
    elif mut_pos + nmer - 1 > len(prot_seq):
        potential_mut_peptide = (
            prot_seq[int(mut_pos) - nmer : int(mut_pos) - 1]
            + mut_amino
            + prot_seq[int(mut_pos) : len(prot_seq)]
        )
        wildtype_peptide = prot_seq[mut_pos - nmer : len(prot_seq)]

        for i in range(0, len(prot_seq) - mut_pos + 1):
            mutpeptide = potential_mut_peptide[i : i + nmer]
            wildtype = wildtype_peptide[i : i + nmer]
            # print(mutpeptide)
            finalpeptide.append(
                (
                    mutpeptide,
                    prot_seq[mut_pos - 1] + str(mut_pos) + mut_amino,
                    Gene,
                    ensembl_transcript,
                    sample_name,
                )
            )
            finalpeptide.append(
                (
                    wildtype,
                    prot_seq[mut_pos - 1] + str(mut_pos) + mut_amino + "_WildType",
                    Gene + "_WildType",
                    ensembl_transcript,
                    sample_name,
                )
            )
    ## if mutation occurs in the middle
    else:
        potential_mut_peptide = (
            prot_seq[int(mut_pos) - nmer : int(mut_pos) - 1]
            + mut_amino
            + prot_seq[int(mut_pos) : int(mut_pos) + nmer - 1]
        )
        wildtype_peptide = prot_seq[int(mut_pos) - nmer : int(mut_pos) + nmer - 1]

        for i in range(0, len(potential_mut_peptide) - nmer + 1):
            mutpeptide = potential_mut_peptide[i : i + nmer]
            wildtype = wildtype_peptide[i : i + nmer]
            finalpeptide.append(
                (
                    mutpeptide,
                    prot_seq[mut_pos - 1] + str(mut_pos) + mut_amino,
                    Gene,
                    ensembl_transcript,
                    sample_name,
                )
            )
            finalpeptide.append(
                (
                    wildtype,
                    prot_seq[mut_pos - 1] + str(mut_pos) + mut_amino + "_WildType",
                    Gene + "_WildType",
                    ensembl_transcript,
                    sample_name,
                )
            )

    return finalpeptide


## this function will create input data for flurry prediction but need to provide what allele (after convert to flurry format (two-digs))
## if the Mhc flurry doesn't support sample allele, it will create a file tell you there is no mhc_flurry support alleles
def createFlurryInputSingleSample(
    mhc_flurry_support, sample_flurry_allele, mut_peptide, flurry_input_file
):
    # mhc_flurry_file = '/home/kh31516/kh31516_Lab_Share_script/IdentifiyNeoAntigene'
    mhc_flurry_list = pd.read_csv(mhc_flurry_support, header=None)
    mhc_flurry_list.columns = ["Allele_name"]
    mhc_flurry_support_list = set(mhc_flurry_list["Allele_name"])
    final_flurry_list = list(sample_flurry_allele & mhc_flurry_support_list)
    with open(flurry_input_file, "w") as w:
        if len(final_flurry_list) == 0:
            w.write("No Mhc_flurry supported alleles")
        else:
            w.write("allele,peptide\n")
            for allele in final_flurry_list:
                # file_name =allele.replace('*','_').replace(':','_')
                for peptide in mut_peptide:
                    w.write(allele + "," + peptide + "\n")


def createInputFlurry(
    mhc_flurry_support, mut_peptide, flurry_output_base, translate_to
):
    ### create MhcFlurry input files
    ## flurry accept two digits of allele name: ex DLA-88*05(two digits), while Yuans' allele names is DLA-88*005

    mhc_flurry_list = pd.read_csv(mhc_flurry_support, header=None)
    mhc_flurry_list.columns = ["Allele_name"]

    if translate_to.upper() == "DOG":
        final_mhc_flurry_list = list(
            mhc_flurry_list[mhc_flurry_list["Allele_name"].str.contains("DLA")][
                "Allele_name"
            ]
        )

    else:
        final_mhc_flurry_list = list(
            mhc_flurry_list[mhc_flurry_list["Allele_name"].str.contains("HLA")][
                "Allele_name"
            ]
        )

    final_mhc_flurry_list = set(final_mhc_flurry_list)

    ## seperate each allele as the output to do the prediction
    for allele in final_mhc_flurry_list:
        # file_name =allele.replace('*','_').replace(':','_')
        flurry_output = open(
            flurry_output_base + "/" + allele + "_" + "mutant_peptides", "w"
        )

        flurry_output.write("allele,peptide\n")
        for peptide in mut_peptide:
            flurry_output.write(allele + "," + peptide + "\n")

        flurry_output.close()


def createRandomPeptide(numberentries, notOverlapPeptide, nmer):
    notOverlapPeptide = set(notOverlapPeptide)
    aa_list = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]
    randomPeptide = []
    for i in range(numberentries):
        candidatePeptide = "".join(random.sample(aa_list, nmer))
        if candidatePeptide not in notOverlapPeptide:
            randomPeptide.append(candidatePeptide)
        # ["".join(random.sample(aa_list,9)) for i in range(numberentries) if random.sample(aa_list,9) not in notOverlapPeptide ]

    return randomPeptide


## the function take a dataframe created by iterating whole mutation data and create nmer peptide by the function createSNVMutPeptide
## it will return a df contains ['Mut_peptides','WT_peptides','Amino_acid','Gene','Ensembl_transcript','Sample_name']
def summarizeMutWTpeptideDf(mut_wt_Df):
    wildtype = mut_wt_Df[mut_wt_Df["Amino_acids"].str.contains("WildType")]
    mutant = mut_wt_Df[~mut_wt_Df["Amino_acids"].str.contains("WildType")]

    if len(wildtype.index) == len(mutant.index):
        print("WildType and Mutant are have same length")

        wildtype = wildtype.reset_index()
        mutant = mutant = mutant.reset_index()
        final_clean_table = pd.merge(
            mutant, wildtype, how="inner", left_index=True, right_index=True
        )
        # target_columns =['Allele_x', 'Peptide_x','Peptide_y' ,'Amino_acids_x', 'Gene_x', 'Ensembl_transcript_x','Sample_name_x']
        if "Allele" in mutant.columns:
            target_columns = [
                "Allele_x",
                "Peptide_x",
                "Peptide_y",
                "Amino_acids_x",
                "Gene_x",
                "Ensembl_transcript_x",
                "Sample_name_x",
            ]
            final_column_names = [
                "Allele",
                "Mut_peptides",
                "WT_peptides",
                "Amino_acid",
                "Gene",
                "Ensembl_transcript",
                "Sample_name",
            ]

        else:
            target_columns = [
                "Peptide_x",
                "Peptide_y",
                "Amino_acids_x",
                "Gene_x",
                "Ensembl_transcript_x",
                "Sample_name_x",
            ]
            final_column_names = [
                "Mut_peptides",
                "WT_peptides",
                "Amino_acid",
                "Gene",
                "Ensembl_transcript",
                "Sample_name",
            ]

        final_clean_table = final_clean_table[target_columns]
        final_clean_table.columns = final_column_names

        return final_clean_table
    else:
        print("Different Length in WildType and Mutant")


def converStringtoTwo(allele):
    ## change DLA-88*005 into DLA-88*05(two digits)
    final_allele = ""
    if ":" not in allele:
        return allele
    elif "New" in allele:
        return allele
    else:
        gene = allele.split("*")[0]
        group = allele.split("*")[1]
        last_dig = group.split(":")[1]
        if gene == "DLA-88":
            first_dig = group.split(":")[0]
            if first_dig[0] == "0" and len(first_dig) == 3:
                first_dig = first_dig[1:]
            final_allele = gene + "*" + first_dig + ":" + last_dig
        elif gene == "DLA-64":
            first_dig = group.split(":")[0]
            if len(last_dig) == 1:
                final_allele = gene + "*" + first_dig + ":" + "0" + last_dig
            else:
                final_allele = gene + "*" + first_dig + ":" + last_dig
        else:
            final_allele = allele

    return final_allele


def converStringtoThree(allele):
    final_allele = ""
    if ":" not in allele:
        return allele
    elif "New" in allele:
        return allele
    else:
        gene = allele.split("*")[0]
        group = allele.split("*")[1]
        last_dig = group.split(":")[1]

    if gene == "DLA-88":
        first_dig = group.split(":")[0]
        if len(first_dig) == 2:
            first_dig = "0" + first_dig
        elif len(first_dig) == 1:
            return allele

        final_allele = gene + "*" + first_dig + ":" + last_dig
    elif gene == "DLA-64":
        first_dig = group.split(":")[0]
        if len(last_dig) == 1:
            final_allele = gene + "*" + first_dig + ":" + "0" + last_dig
        else:
            final_allele = gene + "*" + first_dig + ":" + last_dig
    else:
        final_allele = allele

    return final_allele


def isNaN(num):
    return num != num


def convertEmptyGene(GeneName):
    if GeneName == "-" or GeneName == "-_WildType":
        GeneName = GeneName.replace("-", "NoGeneName")

    return GeneName


def transformBindingScores(x):
    y = x - 1
    value = 50000 ** (-y)
    return value


def check_gene_trancript_loc(info_list):
    # the function identify the Ensembl gene id and transcript id location
    gene_transcript_loc = ()
    gene_symbol = None
    for i in range(len(info_list)):
        if "gene:" in info_list[i]:
            gene_index = i
        elif "transcript:" in info_list[i]:
            trans_index = i
        elif "gene_symbol:" in info_list[i]:
            gene_symbol = i

    gene_transcript_loc = (gene_index, trans_index, gene_symbol)
    return gene_transcript_loc


def createPEPdatabse(pepGTFfile):
    ##### processinig human create a summary df containing transcript id, gene ensembl id, gene name, sequence ######
    ##### create a transcript peptide sequence database from PEP GTF files
    human_gene_transcript_summ = []
    gene_ensembl = ""
    transcript = ""
    gene_symbol = ""
    sequence = ""
    header = ""
    with open(pepGTFfile, "r") as f:
        total_file = f.read().split("\n")[:-1]
        for i, each_line in enumerate(total_file):
            if i != len(total_file) - 1:
                if each_line.startswith(">"):
                    each_sum = [header, transcript, gene_ensembl, gene_symbol, sequence]

                    if each_sum != ["", "", "", "", ""]:
                        human_gene_transcript_summ.append(
                            [header, transcript, gene_ensembl, gene_symbol, sequence]
                        )

                    info = each_line.split(" ")
                    header = info[0].split(">")[1]
                    index_info = check_gene_trancript_loc(info)
                    # print(index_info)
                    gene_ensembl = info[index_info[0]].split(":")[1].split(".")[0]
                    transcript = info[index_info[1]].split(":")[1].split(".")[0]

                    if index_info[2] != None:
                        gene_symbol = info[index_info[2]].split(":")[1]

                    else:
                        gene_symbol = "-"

                    sequence = ""
                else:
                    sequence += each_line.split("\n")[0]
            else:
                human_gene_transcript_summ.append(
                    [header, transcript, gene_ensembl, gene_symbol, sequence]
                )

    human_gene_transcript_df = pd.DataFrame(human_gene_transcript_summ)
    human_gene_transcript_df.columns = [
        "Header",
        "Ensembl_transcript",
        "Ensembl_gene",
        "Gene_name",
        "Sequence",
    ]
    human_gene_transcript_df.index = human_gene_transcript_df["Ensembl_transcript"]
    return human_gene_transcript_df


def createSwissportDb(Swissport_file):
    db_sum = []
    gene_symbol = ""
    sequence = ""
    header = ""
    with open(Swissport_file, "r") as f:
        total_file = f.read().split("\n")[:-1]
        for i, each_line in enumerate(total_file):
            if i != len(total_file) - 1:
                if each_line.startswith(">"):
                    each_sum = [header, gene_symbol, sequence]

                    if each_sum != ["", "", ""]:
                        db_sum.append([header, gene_symbol, sequence])

                    info = each_line.split(" ")
                    header = info[0].split(">")[1]
                    gene_info = each_line.split("GN=")
                    if len(gene_info) > 1:
                        gene_symbol = each_line.split("GN=")[1].split(" ")[0]
                    else:
                        gene_symbol = "NoGeneName"

                    sequence = ""
                else:
                    sequence += each_line.split("\n")[0]
            else:
                db_sum.append([header, gene_symbol, sequence])

        db_sum_df = pd.DataFrame(db_sum)
        db_sum_df.columns = ["Header", "Gene_name", "Sequence"]

    return db_sum_df


def createDictforHumanDogSearch(clean_translate_table):
    total_dict = {}
    #### create human_dog translation dict for future search
    human_dog_pos_dict = {}
    dog_human_pos_dict = {}
    human_aa_dict = {}
    dog_aa_dict = {}

    for i, j in enumerate(clean_translate_table):
        gene = j[1]
        human_pos = j[2]
        dog_pos = j[3]
        human_aa = j[4]
        dog_aa = j[5]

        if gene in dog_human_pos_dict.keys():
            human_dog_pos_dict[gene][human_pos] = dog_pos
            dog_human_pos_dict[gene][dog_pos] = human_pos
            human_aa_dict[gene][human_pos] = human_aa
            dog_aa_dict[gene][dog_pos] = dog_aa

        else:
            human_dog_pos_dict[gene] = {}
            human_dog_pos_dict[gene][human_pos] = dog_pos

            dog_human_pos_dict[gene] = {}
            dog_human_pos_dict[gene][dog_pos] = human_pos

            human_aa_dict[gene] = {}
            human_aa_dict[gene][human_pos] = human_aa

            dog_aa_dict[gene] = {}
            dog_aa_dict[gene][dog_pos] = dog_aa

    total_dict["human_dog_pos_dict"] = human_dog_pos_dict
    total_dict["dog_human_pos_dict"] = dog_human_pos_dict
    total_dict["human_aa_dict"] = human_aa_dict
    total_dict["dog_aa_dict"] = dog_aa_dict

    return total_dict


## current function only care about SNV and fs, and these two is the only we can do the dog_human comparison
## the consequence results (fs,snv) are derived from annovar annotation results, other annotation might not work
def identify_species_counterparts(
    gene_mut_info,
    human_dog_pos_dict,
    dog_human_pos_dict,
    human_aa_dict,
    dog_aa_dict,
    translate_to,
):
    gene_name = gene_mut_info.split("_")[0]
    mut_info = gene_mut_info.split("_")[1]

    ## if we want to translate to dog, then we need to use human_dog_pos_dict and vise versa
    if translate_to.upper() == "DOG":
        ref_dict = human_dog_pos_dict
        alt_dict = dog_human_pos_dict
        other_species_aa_dict = dog_aa_dict
    else:
        ref_dict = dog_human_pos_dict
        alt_dict = human_dog_pos_dict
        other_species_aa_dict = human_aa_dict

    pos = 0
    other_counterparts = " "

    ## extract mutation data
    ## consider three situation, SNV (stop_gain), fs, and other can't process (delines can't process because we don't know the downstream)

    if "fs" in mut_info:
        if re.search(r"([A-Za-z])(\d+)([A-Za-z]*)(fs)", mut_info):
            fs_info = re.search(r"([A-Za-z])(\d+)([A-Za-z]*)(fs)", mut_info)
            wt = fs_info.group(1)
            pos = int(fs_info.group(2))
            mut = ""
            # final_mut_info = wt+loc+fs
            situtation = "fs"
        else:
            other_counterparts = "No Counterparts"

    ## SNV or stop gain
    elif re.search(r"([A-Za-z])(\d+)([A-Z])", mut_info):
        SNV_info = re.search(r"([A-Za-z])(\d+)([A-Z])", mut_info)
        wt = SNV_info.group(1)
        pos = int(SNV_info.group(2))
        mut = SNV_info.group(3)

        situtation = "SNV"
    else:  ## if not SNV or fs types, just directly skip it (include delines)
        other_counterparts = "No Counterparts"

    if other_counterparts == " ":
        if gene_name in alt_dict.keys():
            # dog_target_gene_pos = set(translate_allign_table.loc[(translate_allign_table.Gene ==target_gene)]['QueryIdx'])
            ## if snv or fs
            if (
                (wt.upper() in common_amino_acid_value.keys())
                and (mut.upper() in common_amino_acid_value.keys())
            ) or ((wt.upper() in common_amino_acid_value.keys()) and mut == ""):
                if pos in ref_dict[gene_name].keys():
                    other_species_pos = ref_dict[gene_name][pos]

                    if other_species_pos in other_species_aa_dict[gene_name]:
                        other_species_wt = other_species_aa_dict[gene_name][
                            other_species_pos
                        ]
                        given_species_mut = mut
                        ## if mutation is synonymous SNV, then just replace the location and use counterparts WT for both aa
                        if wt == mut:
                            other_counterparts = (
                                gene_name
                                + "_"
                                + other_species_wt
                                + str(other_species_pos)
                                + other_species_wt
                            )

                        elif given_species_mut == other_species_wt:
                            other_counterparts = "No Mutation"

                        # human_info = translate_allign_table.loc[(translate_allign_table.Gene ==target_gene) & (translate_allign_table.QueryIdx == pos)][["RefIdx","RefAA"]]
                        # human_pos = list(human_info['RefIdx'])[0]
                        # human_aa = list(human_info['RefAA'])[0]
                        # print(human_pos)
                        else:
                            if situtation == "SNV":
                                ## if the mutation is synonymous SNV

                                other_counterparts = (
                                    gene_name
                                    + "_"
                                    + other_species_wt
                                    + str(other_species_pos)
                                    + given_species_mut
                                )
                            elif situtation == "fs":
                                other_counterparts = (
                                    gene_name
                                    + "_"
                                    + other_species_wt
                                    + str(other_species_pos)
                                    + "fs"
                                )

                    else:
                        other_counterparts = "No Counterparts"
                        #'Another species doesnt have the pos'
                else:
                    other_counterparts = "No Counterparts"
                    #'Current pos cannot align to another species'
            else:
                other_counterparts = "No Counterparts"
        else:
            other_counterparts = "No Counterparts"
            #'Another species has no '+gene_name+' in the databases'

    return other_counterparts


common_amino_acid_value = collections.OrderedDict(
    sorted(
        {
            "A": 0.05,
            "R": 0.05,
            "N": 0.05,
            "D": 0.05,
            "C": 0.05,
            "E": 0.05,
            "Q": 0.05,
            "G": 0.05,
            "H": 0.05,
            "I": 0.05,
            "L": 0.05,
            "K": 0.05,
            "M": 0.05,
            "F": 0.05,
            "P": 0.05,
            "S": 0.05,
            "T": 0.05,
            "W": 0.05,
            "Y": 0.05,
            "V": 0.05,
            "X": 0.05,
        }.items()
    )
)
