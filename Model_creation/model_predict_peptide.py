#!/usr/bin/env python3
"""
This script utilizes the Artificial Neural Network (ANN) model that we've developed to predict the binding potential of peptides. We have the flexibility to set a probability cutoff, enabling us to determine whether the peptide is likely to bind or not."
"""

import sys

from scipy import stats

sys.path.append(r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\Model_creation")
from keras.models import load_model

from Blosum62_Model_modules import *

top_binder = 0.1
model_name = "DLA-50101_best_10batch_model.h5"
base = r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\example_data"
peptide_file = base + "/" + "Random100TwithoutBinder.txt"
new_model = keras.models.load_model(base + "/" + model_name)

### create random peptide and turn into DataFrame
randomPeptide = importRandomPeptide(peptide_file)
predictRandom = blosumEncoding(randomPeptide).values

### Test the model performance with random dataset
test_new_pred = new_model.predict(predictRandom)
total_rank = (stats.rankdata(-test_new_pred, "dense") / len(test_new_pred)) * 100
resultPeptide = np.array(["".join(randomPeptide[i]) for i in range(len(randomPeptide))])

outputdf = pd.DataFrame({"Sequence": resultPeptide, "Total_rank": total_rank})

writeOuput(
    outputdf,
    cutOff=top_binder,
    outputName=base
    + "/"
    + "Blosum_"
    + top_binder
    + "_DLA88-50101_MHC_model_redict_result.txt",
)
