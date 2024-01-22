#!/usr/bin/env python3
"""
This script utilizes the Artificial Neural Network (ANN) model that we've developed to predict the binding potential of peptides. We have the flexibility to set a probability cutoff, enabling us to determine whether the peptide is likely to bind or not."
"""

import os
import sys

from scipy import stats

sys.path.append(r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\Model_creation")
from keras.models import load_model

from Blosum62_Model_modules import *

base = r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\example_data"
os.chdir(base)
model_name = "DLA-50101_best_10batch_model.h5"
peptide_file = base + "/" + "Mut_peptide.txt"
meta_table_file = base + '/' + 'Mut_peptide_meta_table.txt'
new_model = keras.models.load_model(base + "/" + model_name)

### create random peptide and turn into DataFrame
peptide_data = pd.read_csv(peptide_file,header=None)
meta_table = pd.read_csv(meta_table_file, sep="\t")
final_peptide = [split_word(makeSameLength(each_peptide)) for each_peptide in peptide_data[0]]

predictpeptide = blosumEncoding(final_peptide).values

### Test the model performance with random dataset
test_new_pred = new_model.predict(predictpeptide)
total_rank = (stats.rankdata(-test_new_pred, "dense") / len(test_new_pred)) * 100


outputdf = pd.DataFrame({"Peptide": peptide_data[0], 
                         "Prob":test_new_pred.flatten(),
                         "Rank":total_rank
                         })

final_table = outputdf.merge(meta_table,on='Peptide',how = 'left')
final_table = final_table[final_table['Gene']!='random']
final_table.to_csv("Peptide_predict_results.txt",sep ="\t", index = False)