# Neoantigen Prediction

This repository contains code and models for predicting peptide-MHC binding affinity, an important step in identifying potential neoantigens for cancer immunotherapy. Neoantigens are tumor-specific antigens encoded by mutated genes that can be recognized by the immune system. Accurately predicting peptide-MHC binding is crucial for neoantigen identification.

This repository contains two main directories for predicting neoantigens:

## Mutant Peptide Prediction

The `Mutant_peptide_prediction` directory is designed to work with mutation data. It can derived from the TOSMIC package or your own mutation data. It generates mutant peptides from the given mutation data and provides the flexibility to predict the binding affinity of these peptides using one of the following methods:

1. MHC_flurry 2.0: A state-of-the-art tool for predicting peptide-MHC binding affinity.
2. NetMhcpan4.1: Another widely used tool for predicting peptide-MHC binding.

Additionally, we have developed our own custom model for peptide-MHC binding prediction, which can also be utilized in this directory.

### Usage

1. Prepare the mutation data in the required format. Ensure that the mutation data includes the following columns: Consequence, Sample_name, Gene_name, Ensembl_transcripts, and Total_protein_change. If the mutation data comes from TOSMIC, the script will create mutant peptides from somatic mutations only.

2. Run the `useSomaticCreateMutPep.py` script to generate mutant peptides from the mutation data. The script will generate several output files:

   - `Mut_peptide.txt`: Contains the generated mutant peptides and includes 100k random peptides for evaluation purposes.
   - `Mut_peptide_meta_table.txt`: Provides additional information about the candidate mutant peptides and includes 100k random peptides for evaluation purposes.
   - `allele + flurry_mut_peptide.txt`: This file serves as the input for MHC_flurry 2.0 prediction.

   **Note:**
   a. Currently, the script supports "nonsynonymous SNV" mutations only for creating mutant peptides.
   b. MHC_flurry supports neoantigen prediction for Dog, Human, and Mouse alleles. In the `useSomaticCreateMutPep.py` script, we extract Dog alleles (currently 75 alleles), and each allele will have its own input file for MHC_flurry 2.0 prediction.

3. Utilize the `allele + flurry_mut_peptide.txt` file as input for MHC_flurry 2.0 prediction to assess peptide-MHC binding affinity.

4. For NetMhcpan prediction, use the `Mut_peptide.txt` file containing the generated mutant peptides.

5. If you need to check which mutation types correspond to the candidate mutant peptides, refer to the `Mut_peptide_meta_table.txt` file.

Please ensure that you have all the necessary dependencies and data files in place before running the script. Detailed instructions and additional information can be found in the provided documentation for further guidance on usage and setup.

## Model Creation Directories

The `Model creation directories` contain code to create our own model based on the Artificial Neural Network (ANN) architecture, implemented using the Keras library. This model is designed for predicting peptide-MHC binding affinity. To enhance the model's prediction accuracy, we utilize the Blosum62 matrix to encode amino acids during the training process.

### Usage

1. Prepare the training data in the appropriate format (refer to the example file for guidance).

2. Train the ANN model using the `Blosum62_MHC_ANN.py` scripts in the respective directories. These scripts implement the Blosum62 matrix encoding to improve the model's performance.

3. Evaluate the performance of the trained model using `model_eval_random_peptide.py`. This evaluation script allows you to assess the model's predictive capabilities on random peptide samples. Additionally, you can use this model to make predictions and select the target rank, such as choosing the top 1% of binding predictions as your target binder. You can then compare the motifs between the predicted peptides and the actual peptides.

The ANN model is created based on the parameters found in MHC_flurry 2.0, ensuring it leverages valuable insights from the state-of-the-art tool to enhance its performance.

Feel free to explore each directory for detailed instructions on data preparation, model training, and prediction procedures. Additional information can be found within the directories to aid you in utilizing and customizing the model according to your specific requirements.

### Requirements

Keras 2.11.0
TensorFlow 2.11.0
NumPy 1.23.5
Pandas >= 1.3
Scikit-learn 1.0
Matplotlib 3.6.2
