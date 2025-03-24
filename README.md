Allows the identification of sequences that significantly increase the number of gene losses needed for reconciliation under a model of gene duplication and loss, the implication being that these are either a) out paralogs, b) a sequecing contamination or c) horizontally transferred.

files:
remove_rogue.py - the main program that sequentially removes sequences from the gene tree and reconciles it with the species tree, providing a table that shows how many "extra" gene loss events are needed to explain the presence of each gene.
example - a dummy dataset
process_rogue_full_pipeline.py - A pipeline to sequentially remove rogue genes until none more exist that require over a certain number of extra loss events to explain
process_rogue_full_pipeline_slurm.sh - a slurm job script to run proces_rogue_full_pipeline.py for a set of gene families as a job array

In future, I might make this work for groups of several genes, to identify horizontally transferred groups of 2, 3, 4 or more genes.
