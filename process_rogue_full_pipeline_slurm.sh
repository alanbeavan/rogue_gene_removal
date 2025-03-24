#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20g
#SBATCH --time=3-00:00:00
#SBATCH --array 15-15
#go to 576 for full run

echo "Running on `hostname`"
cd ${SLURM_SUBMIT_DIR}
family=`sed "${SLURM_ARRAY_TASK_ID}q;d" remove_rogue_post_gum/all_families.txt`
echo $family
echo ../rogue_remove_pipeline/pre_remove/ortho_alignments/${family}.fa
echo ../rogue_remove_pipeline/pre_remove/true_orthologs/${family}.nwk
#rogue_remove_pipeline/pre_remove/true_orthologs

python3 process_rogue_full_pipeline.py  rogue_remove_pipeline/pre_remove/ortho_alignments/${family}.fa rogue_remove_pipeline/pre_remove/true_orthologs/${family}.nwk 6 rogue_remove_pipeline/${family}



sleep 40
echo "Finished job now"
