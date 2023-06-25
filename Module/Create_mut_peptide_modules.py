#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import collections
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

neo_antigen_module = "/Volumes/Research/GitHub/NeoantigenIdentify/Module"
# r'C:\Users\abc73\Documents\GitHub\NeoantigenIdentify'
# '/Volumes/Research/GitHub/NeoantigenIdentify'
sys.path.append(neo_antigen_module)


def convertEmptyGene(GeneName):
    if GeneName == "-" or GeneName == "-_WildType":
        GeneName = GeneName.replace("-", "NoGeneName")

    return GeneName


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
        print("There is an error in sample ")


## this function will process annovar output and extract all the protein changes in the annovar file (not only one protein change)
## it will remove Ensembl gene that located in retro_gene list idenitfied in pan-cancer paper
## even the annovar contains sample name, it will overwrite the orignal one
## the final return table contains the following
## Line','Consequence','Gene_name','Chrom','Start','End','Ref','Alt','Sample_name','Ensembl_gene','Ensembl_transcripts','Total_protein_change'
## the reason I keep line is for future used for VAF calculation
def processAnnovar(annovar_gene_file, sample_name):
    # retro_gene_list = pd.read_csv(retro_gene_file, sep="\n", header=None)[0].tolist()
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
    target_annovar_info = target_annovar_info.explode(
        ["Gene_name", "Ensembl_gene", "Ensembl_transcripts", "Total_protein_change"]
    )
    target_annovar_info = target_annovar_info.loc[
        (target_annovar_info["Ensembl_transcripts"] != "No_Info_Provided")
        & (target_annovar_info["Total_protein_change"] != "No_Info_Provided")
    ]
    # ### filter retrogene list with given ensembl id
    # target_annovar_info = target_annovar_info[
    #     ~target_annovar_info.Ensembl_gene.isin(retro_gene_list)
    # ]
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


### The function will create a list containing mutant peptides and wild-type peptides, Gene name, ensembl_transcript and smaple name
### A function used in create_mut_pep_sum function
def createSNVMutPeptide(
    prot_seq, nmer, mut_pos, mut_amino, Gene, ensembl_transcript, sample_name
):
    finalpeptide = []

    ## if mutation occurs in the begining
    if mut_pos < nmer:
        potential_mut_peptide = (
            prot_seq[0 : mut_pos - 1] + mut_amino + prot_seq[mut_pos : mut_pos + nmer]
        )
        wildtype_peptide = prot_seq[0 : mut_pos + nmer]

        for i in range(0, mut_pos):
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

    ## if the mutation occurs in the end
    elif mut_pos + nmer - 1 > len(prot_seq):
        potential_mut_peptide = (
            prot_seq[int(mut_pos) - nmer : int(mut_pos) - 1]
            + mut_amino
            + prot_seq[int(mut_pos) : len(prot_seq)]
        )
        wildtype_peptide = prot_seq[mut_pos - nmer : len(prot_seq)]

        for i in range(len(prot_seq) - mut_pos + 1):
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
    ## if mutation occurs in the middle
    else:
        potential_mut_peptide = (
            prot_seq[int(mut_pos) - nmer : int(mut_pos) - 1]
            + mut_amino
            + prot_seq[int(mut_pos) : int(mut_pos) + nmer - 1]
        )
        wildtype_peptide = prot_seq[int(mut_pos) - nmer : int(mut_pos) + nmer - 1]

        for i in range(len(potential_mut_peptide) - nmer + 1):
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


## the function is same as df.explode but in case pandas version < 1.3
## there is a extra function
def unnesting(df, explode):
    idx = df.index.repeat(df[explode[0]].str.len())
    df1 = pd.concat(
        [pd.DataFrame({x: np.concatenate(df[x].values)}) for x in explode], axis=1
    )
    df1.index = idx

    return df1.join(df.drop(explode, 1), how="left")


## the function take a dataframe created by iterating whole mutation data and create nmer peptide by the function createSNVMutPeptide
## it will return a df contains ['Mut_peptides','WT_peptides','Amino_acid','Gene','Ensembl_transcript','Sample_name']
def summarize_mut_wt_peptideDf(mut_wt_Df):
    wildtype = mut_wt_Df[mut_wt_Df["Amino_acids"].str.contains("WildType")]
    mutant = mut_wt_Df[~mut_wt_Df["Amino_acids"].str.contains("WildType")]

    if len(wildtype.index) == len(mutant.index):
        print("WildType and Mutant are have same length")

        wildtype = wildtype.reset_index()
        mutant = mutant = mutant.reset_index()
        final_clean_table = pd.merge(
            mutant, wildtype, how="inner", left_index=True, right_index=True
        )
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


### This function uses somatic mutation files (but only SNV) to create to create two tables .
### the mutation results must contain "gene_name", "ensembl_transcripts", "total_protein_change", "sample_name" columns to create the tables, so it can come from pipeline results of TOSMIC or ml results of TOSMIC
### extract_data argument can be all (no filtering), m (from machine learning model from TOSMIC) or p (from pipleine filtering from TOSMIC)
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
    if len(list(set(require_col) & set(mut_data_col))) != 4:
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
        mut_wt_sum_df = summarize_mut_wt_peptideDf(final_sum)

        if includeRandom:
            random_peptide_file = data_source_folder / "Pan_cancer_random_peptide.txt"
            random_peptide = pd.read_csv(random_peptide_file, sep=",")
            final_sum = pd.concat([final_sum, random_peptide])

        final_return["summary"] = mut_wt_sum_df
        final_return["meta"] = final_sum

    return final_return


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
