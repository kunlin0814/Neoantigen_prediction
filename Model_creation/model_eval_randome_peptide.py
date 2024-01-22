#!/usr/bin/env python3
"""
This script is used to test the Artificial Neural Network (ANN) model created for motif identification in experimental data. 
The purpose is to assess whether the model can accurately identify motifs similar to those observed in the actual experiment data.
We can select a target rank (e.g., top 1% binding) to identify potential binders. 
By comparing the motifs between the predicted peptides and the actual peptides, we can evaluate the model's performance in motif recognition and binder prediction.
random_peptide_file contains 100k random peptides 
"""

import sys

from scipy import stats

sys.path.append(r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\Model_creation")
from keras.models import load_model

from Blosum62_Model_modules import *

top_binder = 0.1
model_name = "DLA-50101_best_10batch_model.h5"
base = r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\example_data"
random_peptide_file = base + "/" + "Random100TwithoutBinder.txt"
new_model = keras.models.load_model(base + "/" + model_name)

### create random peptide and turn into DataFrame
randomPeptide = pd.read_csv(random_peptide_file,sep="\t",header=None)
final_peptide = [split_word(makeSameLength(each_peptide)) for each_peptide in randomPeptide[0]]
predictRandom = blosumEncoding(final_peptide).values

### Test the model performance with random dataset
test_new_pred = new_model.predict(predictRandom)
total_rank = (stats.rankdata(-test_new_pred, "dense") / len(test_new_pred)) * 100
outputdf = pd.DataFrame({"Sequence": randomPeptide[0], "Total_rank": total_rank,"Prob":test_new_pred.flatten()})

outputdf.to_csv("Peptide_predict_results.txt",sep ="\t", index = False)