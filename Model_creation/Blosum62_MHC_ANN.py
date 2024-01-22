#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will take a peptide binding data to create an ANN model.
We can use this model to make an prediction and select the target rank as our binder (ex: top 1% binding as our target binder) and compare the motif between predicted peptides and actual peptides
The ANN model is created based on the parameters found on the MHC_flurry 2.0

The model can take peptides with length from 8-15 as training input data
"""
import os
import sys

sys.path.append(r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\Model_creation")
#'/Volumes/Research/GitHub/NeoAntigeneModelCreation/')

from Blosum62_Model_modules import *

test_ratio = 0.2
val_size = 0.1
base = r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\example_data"
os.chdir(base)
dataset = pd.read_csv(base + "/" + "Final_DLA88-50101_training_data.txt", sep="\t")
peptides = dataset['peptide'].apply(lambda x:x.upper())
final_peptide = [split_word(makeSameLength(each_peptide)) for each_peptide in peptides]
total_data = extractdata(dataset)
y = dataset.iloc[:, 2].values
X = blosumEncoding(final_peptide).values

# from sklearn.preprocessing import LabelEncoder, OneHotEncoder
# Splitting the dataset into the Training set and Test set
from sklearn.model_selection import train_test_split

# create train, test and validation dataset
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=test_ratio, random_state=0, shuffle=True, stratify=y
)
X_train, X_val, y_train, y_val = train_test_split(
    X_train, y_train, test_size=val_size, random_state=0, shuffle=True
)
# Dense object will take care to initialize the random number close to 0 ( first ANN step)
ANN_classifier = build_classifier()
# add earlystopping to prevent overfitting
es = EarlyStopping(monitor="val_loss", mode="auto", verbose=0, patience=50)
# check the model performance and save the best model
mc = ModelCheckpoint(
    base + "/" + "DLA-50101_best_10batch_model.h5",
    monitor="val_loss",
    mode="min",
    save_best_only=True,
)
history = ANN_classifier.fit(
    X_train,
    y_train,
    validation_data=(X_val, y_val),
    verbose=0,
    batch_size=10,
    epochs=300,
    callbacks=[mc, es],
)

## we can visualize the training process
plt.plot(history.history["loss"], label="train")
plt.plot(history.history["val_loss"], label="validation")
plt.legend()
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.show()
history.history["loss"]


y_pred = ANN_classifier.predict(X_test)
# use 0.5,0.7, 09 to as the cutoff to decide the instance is a positive instance or negative instance
## the output give the probability and we apply the cut off to decide if the peptide is a binder or not

thresheld = [0.5, 0.7, 0.9]
y_pred_value = (y_pred > thresheld).tolist()

y_pred_list = []
for i in thresheld:
    y_pred_value = (y_pred > i).tolist()
    y_pred_list.append(y_pred_value)

# Making the Confusion Matrix
from sklearn.metrics import confusion_matrix

cm_list = []
for i in range(len(y_pred_list)):
    cm = confusion_matrix(y_test, y_pred_list[i], labels=[1, 0])
    cm_list.append(cm)
    print("", "True", "False", sep="\t")
    print("True", cm[0, 0], cm[0, 1], sep="\t")
    print("Flase", cm[1, 0], cm[1, 1], sep="\t")
    print("-" * 20)
    print("accuracy:", sklearn.metrics.accuracy_score(y_pred_list[i], y_test))


test_loss, test_acc = ANN_classifier.evaluate(X_test, y_test, verbose=2)
print(test_acc)


"""
K fold cross validation, 
We can test the average of the accuracy and std of the accuracy 
we need to combine keras to scikit_learn 
"""


from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import (StratifiedKFold, cross_val_score,
                                     train_test_split)

kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=7)
classifier = KerasClassifier(build_fn=build_classifier, batch_size=50, epochs=200)
accuracies = cross_val_score(
    estimator=classifier, X=X_train, y=y_train, cv=kfold, n_jobs=-1, scoring="accuracy"
)
mean = accuracies.mean()
variance = accuracies.std()


# ## Use GridSearchCV to find the best parameter
# from sklearn.model_selection import GridSearchCV

# classifier = KerasClassifier(build_fn=build_classifier)
# # create dict to tune the hyper-parameter to find the best parameter
# parameters = {
#     "batch_size": [10, 25, 30],  # we can choose batch_size we want
#     "nb_epoch": [100, 200, 300],
# }  # tune the archticture of the ann

# grid_search = GridSearchCV(
#     estimator=classifier,  # classifier = ANN
#     param_grid=parameters,  # test all combination
#     scoring="accuracy",
#     cv=5,  # k-fold validation
# )

# # Identify the best parameters for the model creation
# # overwrite the original grid_search
# grid_search = grid_search.fit(
#     X_train, y_train
# )  # fit the grid_search to the training set
# best_parameters = grid_search.best_params_
# best_accuracy = grid_search.best_score_
