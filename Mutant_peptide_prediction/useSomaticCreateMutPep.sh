#!/bin/bash       
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (4)
#SBATCH --mem=50G                   # Job memory limit (10 GB)
#SBATCH --time=10:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds

scripts='/work/szlab/Kun-Lin_pipeline/Allele_specific_model/Mutant_peptide_prediction'
mut_file='/work/szlab/Kun-Lin_pipeline/Allele_specific_model/example_data/Total_mutation.txt'
output_base='/scratch/kh31516/GR_neoantigen/source'

source activate py38

# file_pattern= sys.argv[1]
# final_df_out = sys.argv[2]
# add_bioproject= sys.argv[3]

python $scripts/useSomaticCreateMutPep.py \
${mut_file} \
${output_base}


# mut_file = sys.argv[1]
# output_base = sys.argv[2]
#r"C:\Users\abc73\Documents\GitHub\Neoantigen_prediction\example_data"
