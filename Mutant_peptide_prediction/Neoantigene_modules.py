#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script contains the function that regular used for Neoantigen identification for Flurry or netMhcPan

"""

import numpy as np
import pandas as pd


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
