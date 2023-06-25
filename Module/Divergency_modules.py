#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 14:34:11 2021

@author: kun-linho
"""
import collections
import random
import sys
from copy import copy
from math import log, log2, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import distance
from six import StringIO

# sys.path.append('/Volumes/Research/GitHub/NeoAntigeneModelCreation')
#     #r'C:\Users\abc73\Documents\GitHub\NeoAntigeneModelCreation')
# from Blosum62_matrix import *

random.seed(1000)
np.random.seed(1000)


def kl_divergence(p, q):
    return sum(p * np.log2(p / q))


# return sum(p[i] * log2(p[i]/q[i]) for i in range(len(p)))


# a function take the peptide_dataframe as input and give back the frequency of each position in the mer
def numberLetter(peptide_df, allele, nmer):
    target_peptide = peptide_df[peptide_df.allele == allele].peptide.drop_duplicates()
    summary = []
    # pd.DataFrame()
    for i in range(nmer):
        each_letter = target_peptide.str[i]
        # df = pd.DataFrame(target_peptide.str[i])
        group_df = each_letter.groupby(each_letter).size()
        # group_df = group_df.sort_index()
        current_index = list(group_df.keys())
        # some amino acid might not show up, so fill up with 0 in the Series
        diff = list(set(peptide_list).symmetric_difference(set(current_index)))
        if len(diff) != 0:
            for each_diff in diff:
                group_df = group_df.append(pd.Series({each_diff: 0}))

        # reset_index(name='position'+str(i+1)).set_index('peptide')
        # sort_df = group_df.sort_values(by='counts_position'+str(i+1),ascending= False)
        group_df = group_df.sort_index()
        # print(group_df)
        summary.append(group_df.rename("position" + str(i + 1)))
        # pd.concat([summary,group_df.rename('position'+str(i+1))],axis=1)
    summary = pd.concat(summary, axis=1)
    count_summary = summary.copy()
    summary = summary.transform(lambda x: x / len(target_peptide.index))
    col_name = ["ratio_" + "position_" + str(i + 1) for i in range(nmer)]
    summary.columns = col_name
    final_summary = pd.concat([count_summary, summary], axis=1)

    # summary = summary.fillna(0)

    return final_summary


def js_divergence(p, q):
    # we don't use this function in this script because scipy package provide the faster speed.
    # but the result is the same.
    m = 0.5 * (p + q)
    value1 = kl_divergence(p, m)
    value2 = kl_divergence(q, m)
    return 0.5 * value1 + 0.5 * value2


## random.permutation
def createRandomMatrixRowColumn(actualDataFrame, column):
    # this function create random matrix with shuffle column and row. This step is the most computational expensive step
    # np.shuffle is computational expensive, so I create a list first and then reassemble into a df and then transpose
    pos = np.random.permutation(column)
    summary = []
    for i in range(len(pos)):
        orig = actualDataFrame[pos[i]]
        # newcopy = orig.copy()
        summary.append(np.random.permutation(orig))

    newdf = pd.DataFrame(summary)
    newdf = newdf.transpose()
    newdf.columns = column
    return newdf


def createRandomMatrix(reference_matrix, column, common_amino_acid_value):
    # this function create random matrix with shuffle column and row
    # but the shuffle row select the same index of the panda. eg. each column select the index for its' row (1,5,9,17,10....same row order)
    # so even it shuffle the rows, each column has the same shuffle order for the row
    sampler = np.random.permutation(reference_matrix.shape[1])
    newcol = reference_matrix.take(sampler, axis=1)
    sampler = np.random.permutation(reference_matrix.shape[0])
    newcolrow = newcol.take(sampler, axis=0)
    newcolrow.index = common_amino_acid_value.keys()
    newcolrow.columns = column

    return newcolrow


def mapDictory(amino_acid_dict, amino_acid_arr):
    # a function that calculate how many amino acid in each column and divide by total number of amino acid
    amino_acid_list = amino_acid_arr.tolist().copy()
    new_amino_acid_dict = amino_acid_dict.copy()

    for eachaa in amino_acid_list:
        new_amino_acid_dict[eachaa.upper()] += 1

    total = sum(new_amino_acid_dict.values())

    for key in new_amino_acid_dict.keys():
        new_amino_acid_dict[key] = new_amino_acid_dict[key] / total

    return new_amino_acid_dict


def peptideListIntoDataFrame(peptide_list):
    # the input file is eg: AARRSER. so the function will create a df for each amino acid and then put into a df.
    # then it will create a probablility matrix for each position (column), here we add a pseuodo count 0.05 to prevent log0 situation.
    # the output is a df with probability matrix for each position (column)
    # use list to append the dataframe row and then tranposes will be much faster than append column
    total_line = []
    for i, j in enumerate(peptide_list):
        each_line = []
        for k in peptide_list[i]:
            each_line.append(k)
        total_line.append(each_line)
    # pd.read_csv(fileName, header = None)
    file = pd.DataFrame(total_line)
    nmer = len(file.columns)
    column = ["Position" + str(i + 1) for i in range(nmer)]
    newdf = file
    newdf.columns = column

    dfsummary = []
    for col in column:
        summary = mapDictory(common_amino_acid_value, newdf[col])
        dfsummary.append(summary.values())
        # dfsummary[col]= summary.values()

    dfsummary = pd.DataFrame(dfsummary).transpose()
    index = list(common_amino_acid_value.keys())

    dfsummary.index = index
    dfsummary.columns = column
    return dfsummary


def processfileIntoDataFrame(fileName):
    # the input file is eg: AARRSER. so the function will create a df for each amino acid and then put into a df.
    # then it will create a probablility matrix for each position (column), here we add a pseuodo count 0.05 to prevent log0 situation.
    # the output is a df with probability matrix for each position (column)
    # use list to append the dataframe row and then tranposes will be much faster than append column
    if isinstance(fileName, pd.DataFrame):
        file = fileName
    else:
        with open(fileName, "r") as f:
            peptide_list = f.read().split("\n")[:-1]
            dfsummary = peptideListIntoDataFrame(peptide_list)

    return dfsummary


# the function take orignal 20 aa frequency matrix and sum the prob within the same group
# ex: it will sum KER prob into a KER group rather than K,E,R seperately
def group_aa(freq_df, grouping_method):
    summary = []
    for each_group in grouping_method:
        if len(each_group) == 1:
            each_group_sum = freq_df.loc[each_group]
            summary.append(each_group_sum)
        else:
            each_letter_sum = pd.Series([], dtype="float64")
            for letter in each_group:
                each_letter_sum = each_letter_sum.add(freq_df.loc[letter], fill_value=0)

            summary.append(each_letter_sum)

    newdf = pd.DataFrame(summary)
    newdf.index = grouping_method
    return newdf


def superTypeDistance(dfA, dfB):
    total_pos = list(dfA.columns)
    each_pos_dis = 0
    for each_pos in total_pos:
        numerator = np.dot(dfA[each_pos], dfB[each_pos])
        denominator = sqrt(sum(dfA[each_pos] ** 2)) * sqrt(sum(dfB[each_pos] ** 2))
        distance = 1 - (numerator / denominator)
        each_pos_dis += distance
        return round(distance, 5)


def calculateDivergency(target_matrix, reference_matrix):
    # if grouping_method!=None:
    #     target_matrix = group_aa(target_matrix,grouping_method)
    #     reference_matrix = group_aa(reference_matrix,grouping_method)

    predictVSActualValue = {}
    # js_summary =[]
    # for col in column:
    # js_divergencevalue = js_divergence(target_matrix[col].values, reference_matrix[col].values)
    js_distance = distance.jensenshannon(target_matrix, reference_matrix)
    js_divergenceValue = js_distance**2
    kl_divergenceValue = kl_divergence(target_matrix.values, reference_matrix.values)
    mean_js_divergenceValue = js_divergenceValue.mean()
    mean_kl_divergenceValue = kl_divergenceValue.mean()
    # output.write(str(kullback_Leibler_divergence)+"\t")

    js_divergenceValue = np.append(js_divergenceValue, mean_js_divergenceValue)
    kl_divergenceValue = np.append(kl_divergenceValue, mean_kl_divergenceValue)
    predictVSActualValue["js_summary"] = js_divergenceValue
    predictVSActualValue["kl_summary"] = kl_divergenceValue

    return predictVSActualValue


def calculatePermutationPvalue(target_matrix, reference_matrix, iteration, allele):
    nmer = len(target_matrix.columns)
    iteration = int(iteration)
    column = reference_matrix.columns.values.tolist()

    aa_index = list(reference_matrix.index.tolist())

    predictVSActualValue = calculateDivergency(target_matrix, reference_matrix)

    js_summary = predictVSActualValue["js_summary"][0:-1]
    kl_summary = predictVSActualValue["kl_summary"][0:-1]

    js_count_value = np.zeros(len(js_summary))
    kl_count_value = np.zeros(len(kl_summary))

    for i in range(int(iteration)):
        permutationMatrix = createRandomMatrixRowColumn(reference_matrix, column)
        permutationMatrix.index = aa_index

        ## the order might need to modify, should use target_matrix to permut or refer to permut ?
        randomMatrixSummary = calculateDivergency(target_matrix, permutationMatrix)

        randomjsSummary = randomMatrixSummary["js_summary"][0:-1]

        randomklSummary = randomMatrixSummary["kl_summary"][0:-1]

        # outrandomjsSummary = randomjsSummary.tolist()
        # outrandomklSummary = randomklSummary.tolist()

        js_count_value = np.where(
            js_summary < randomjsSummary, js_count_value, js_count_value + 1
        )
        kl_count_value = np.where(
            kl_summary < randomklSummary, kl_count_value, kl_count_value + 1
        )

    js_count_value = np.append(js_count_value, np.mean(js_count_value))
    kl_count_value = np.append(kl_count_value, np.mean(kl_count_value))

    df_columns = ["position" + str(i + 1) for i in range(nmer)]
    df_columns.append("average")

    # addition_columns=['Average','Methods','percentage','allele']
    # df_columns.extend(column)
    # df_columns.extend(addition_columns)

    # df_count = pd.DataFrame([js_count_value,kl_count_value],columns=df_columns)
    # df_count['methods'] = ['js','kl']
    # df_count['allele'] = allele

    # df_count.to_csv(output_folder+'/'+allele+'_count_summary.txt',
    #                  sep="\t",index=False)

    js_pvalue = js_count_value / iteration
    kl_pvalue = kl_count_value / iteration
    # previous already calculated mean, so no need to have mean
    # js_pvalue = np.append(js_pvalue,np.mean(js_pvalue))
    # kl_pvalue = np.append(kl_pvalue,np.mean(kl_pvalue))

    df_pvalue = pd.DataFrame([js_pvalue, kl_pvalue], columns=df_columns)
    df_pvalue["methods"] = ["js", "kl"]
    df_pvalue["allele"] = allele

    # df_pvalue.to_csv(output_folder+'/'+allele+'_pvalue_summary.txt',
    #                  sep="\t",index=False)

    print("successfully calculate permutation p value")
    return df_pvalue


COMMON_AMINO_ACIDS = collections.OrderedDict(
    sorted(
        {
            "A": "Alanine",
            "R": "Arginine",
            "N": "Asparagine",
            "D": "Aspartic Acid",
            "C": "Cysteine",
            "E": "Glutamic Acid",
            "Q": "Glutamine",
            "G": "Glycine",
            "H": "Histidine",
            "I": "Isoleucine",
            "L": "Leucine",
            "K": "Lysine",
            "M": "Methionine",
            "F": "Phenylalanine",
            "P": "Proline",
            "S": "Serine",
            "T": "Threonine",
            "W": "Tryptophan",
            "Y": "Tyrosine",
            "V": "Valine",
            "X": "Unknown",
        }.items()
    )
)

COMMON_AMINO_ACIDS_WITH_UNKNOWN = copy(COMMON_AMINO_ACIDS)
COMMON_AMINO_ACIDS_WITH_UNKNOWN["X"] = "Unknown"

AMINO_ACID_INDEX = dict(
    (letter, i) for (i, letter) in enumerate(COMMON_AMINO_ACIDS_WITH_UNKNOWN)
)

AMINO_ACIDS = list(COMMON_AMINO_ACIDS_WITH_UNKNOWN.keys())

BLOSUM62_MATRIX = (
    pd.read_csv(
        StringIO(
            """
  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0  0
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3  0
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  0
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  0
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1  0
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  0
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3  0
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3  0
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1  0
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1  0
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1  0
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2  0
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0  0 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3  0
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1  0
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4  0
   X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1
   """
        ),
        sep="\s+",
    )
    .loc[AMINO_ACIDS, AMINO_ACIDS]
    .astype("int8")
)
assert (BLOSUM62_MATRIX == BLOSUM62_MATRIX.T).all().all()


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
