#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 12:25:09 2021

@author: kun-linho
"""
import collections
import io
import os
import random
from copy import copy
from pathlib import Path

import keras
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
import tensorflow as tf
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.layers import Dense, Dropout
from keras.models import Sequential, load_model
from keras.wrappers.scikit_learn import KerasClassifier
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from six import StringIO
from sklearn.metrics import (classification_report, confusion_matrix,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.utils import resample
from tensorflow.keras import backend


### we can make the length as long as 15 or shorter, like Mhcflurry 2.0,
def makeSameLength(each_peptide):
    first = each_peptide
    middle = each_peptide.center(15, "X")
    end = ""

    while len(first) < 15:
        first += "X"
        end += "X"
    end = end + each_peptide[-len(each_peptide) :]
    return first + middle + end


### the random peptide is created from pan_cancer mutations, total are 100k random peptides, so that we don't have to create it again
def importRandomPeptide(peptideData):
    with open(peptideData, "r") as f:
        data = f.read().split("\n")[:-1]
    total_list = []
    for line in data:
        small_list = []
        for each in line:
            small_list.append(each)
        total_list.append(small_list)

    return total_list


def extractdata(dataset):
    total_data = (
        dataset.iloc[:, 3:12].apply(lambda x: x.astype(str).str.upper()).values.tolist()
    )

    return total_data


def extractbinder(dataset):
    binderData = dataset[dataset["Status"] == 1].iloc[:, 3:12].values.tolist()
    return binderData


def encodeWithBLOSUM62(amino_acids):
    ## a function that returns the blosum62 vector given a certain aa
    return list(BLOSUM62_MATRIX[amino_acids].values)


def blosumEncoding(data):
    ## a function that creates amino acid sequence that encode by blosum62 matrix
    total_data_row = []
    for each in data:
        eachrow = []
        for aa in each:
            eachrow = eachrow + encodeWithBLOSUM62(aa)

        total_data_row.append(eachrow)

    return pd.DataFrame(total_data_row)


def build_classifier():
    classifier = Sequential()  # use the sequential layer
    ## init = kernel_initializer
    classifier.add(
        Dense(units=80, kernel_initializer="uniform", activation="tanh", input_dim=189)
    )
    classifier.add(Dropout(rate=0.5))
    classifier.add(
        Dense(units=80, kernel_initializer="uniform", activation="tanh", input_dim=189)
    )
    classifier.add(Dropout(rate=0.5))
    classifier.add(
        Dense(units=80, kernel_initializer="uniform", activation="tanh", input_dim=189)
    )
    ## if we deal with more than 2 categories, the activation function needs to use softmax
    classifier.add(Dense(units=1, kernel_initializer="uniform", activation="sigmoid"))
    opt = keras.optimizers.RMSprop(learning_rate=0.001)
    # Compiling the ANN
    classifier.compile(optimizer=opt, loss="binary_crossentropy", metrics=["accuracy"])

    return classifier


def createRandomPeptide(BinderData, NumberData, withBinder):
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
    # aa_list_new= np.array(aa_list)
    total_list = set()
    if withBinder:
        total_list = set(BinderData)
    for _ in range(NumberData):
        each_create = "".join(random.sample(aa_list, 9))
        if each_create not in total_list:
            total_list.add(each_create)

    total_list = list(total_list)

    return total_list


def check_legit(peptide):
    for i in peptide:
        if i.upper() not in COMMON_AMINO_ACIDS.keys():
            return False
    return True


def writeOuput(outputdf, cutOff, outputName):
    final_list = outputdf.Sequence[outputdf.Total_rank < cutOff].values.tolist()
    with open(outputName, "w") as f:
        for each in final_list:
            f.write(each + "\n")


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
