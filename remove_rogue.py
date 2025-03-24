#!/usr/bin/env python3.6
"""
remove_rogue.py

This is a heuristic approach to remove a tip from a gene tree if that tip is
either horrizontally transferred or a sequnecing contaminent.

To do this, for each gene duplication in a processed generax family, get the
clade in which it duplicated, then sequentially remove each descendent tip and
again take the clade of duplication.

If the clade of duplication differs significantly, we say that the sequence
causing this is horrizontally transferred or a contamient, so remove it, which
will mean that when inferreing ancestral gene counts, it will make more sense.

The question is how significantly is significantly - I think it depends on how
many losses needs to be inferred for that species to have the gene. I guess
that means I have to calculate the number of losses needed to explain each
situation.

So I need to measure the cost (number of lossese) trees with each tip removed.
I can probably do this by interacting with a secondary program (ALE for
example).

OK - so generax can do this with --strategy EVAL. - keep removing the tip and
run generax with EVAL, saving each result to a table, then if it improves by,
say, 50% (but more than 2 losses) (from each duplication), then the sequence
is removed.

So I will need to set up a families file each time.

You get lines like this
Node_Phaeodactylum.tricornutum_Pseudonitzschia.multistriata_0 S=66 SL=48 D=0 T=0 TL=0
from which I can count the number of SL= to get the number of losses
"""

import ete3
import glob
import os
from random import randint
import shutil
import subprocess
import sys
import my_module as mod

def get_args():
    """Get user arguments."""
    if not len(sys.argv) == 5:
        print("USAGE: python3 remove_rogue.py alignment gene_tree species_tree output_prefix")
        sys.exit()
    return sys.argv[1:]

def set_up_dummy_tree_and_seqs(t, leaf, tip_names, seqs, gr_dir):
    """Write the tree for the dummy generax analysis."""
    to_prune = [x for x in tip_names if x != leaf.name]
    pruned = t.copy()
    pruned.prune(to_prune)

    pruned.write(format = 1, outfile = gr_dir + "/tree.nwk")
    with open(gr_dir + "/seqs.fa", "w", encoding = "utf8") as out:
        for key, value in seqs.items():
            if key in to_prune:
                out.write(">" + key + "\n" + value + "\n")
            #else:
                #print(key + " gone")
def set_up_dummy_families(gr_dir):
    """set up the dummy generax "families" file."""
    with open(gr_dir + "/families.txt", "w", encoding = "utf8") as out:
        out.write("[FAMILIES]\n")
        out.write("-temp\n")
        out.write("alignment = " + gr_dir + "/seqs.fa\n")
        out.write("starting_gene_tree = " + gr_dir + "/tree.nwk\n")
        out.write("subst_model = LG\n")

def main():
    """Do the things."""
    ali_file, gt_file, st_file, output = get_args()
    unique = ""
    for i in range(15):
        unique = unique + str(randint(0,9))
    gr_dir = "temp_remove_rogue_" + unique
    results = {}
    seqs = mod.read_fasta(ali_file)
    t = ete3.Tree(gt_file, format = 1)
    #iterate though each tip and set up the generax --strategy EVAL
    #without said tip using subprocess or something
    tip_names = []
    for leaf in t:
        tip_names.append(leaf.name)
    
    #First, do the reconciliation without removing anything
    if not os.path.exists(gr_dir):
        os.makedirs(gr_dir)
    with open(gr_dir + "/families.txt", "w", encoding = "utf8") as out:
        out.write("[FAMILIES]\n")
        out.write("-temp\n")
        out.write("alignment = " + ali_file + "\n")
        out.write("starting_gene_tree = " + gt_file + "\n")
        out.write("subst_model = LG\n")
    
    
    print("generax -r UndatedDL -f " + gr_dir + "/families.txt -s" + st_file + " --geneSearchStrategy EVAL -p" + gr_dir + "/GeneRax")
    i = 0
    worked = 0
    while i < 10:
        ans = subprocess.call(["generax",
                               "-r", "UndatedDL",
                               "-f", gr_dir + "/families.txt",
                               "-s", st_file,
                               "--geneSearchStrategy", "EVAL",
                               "-p", gr_dir + "/GeneRax"])
        if ans == 0:
            worked = 1
            break
        i += 1
    if not worked:
        print("Generax didn't work for some reason - the temporary files are still there so go check them out.")
        print("The directory is " + gr_dir)
        sys.exit()
    
    for line in open(gr_dir + "/GeneRax/reconciliations/temp_eventCounts.txt", "r").readlines():
        if "SL:" in line:
            results["original"] = int(line.split(":")[1])
            break




    for leaf in t:
        #write the dummy tree
        set_up_dummy_tree_and_seqs(t, leaf, tip_names, seqs, gr_dir)
        #and write the families file
        set_up_dummy_families(gr_dir)

        #run generax
        i = 0
        worked = 0
        while i < 10:
            ans = subprocess.call(["generax",
                                   "-r", "UndatedDL",
                                   "-f", gr_dir + "/families.txt",
                                   "-s", st_file,
                                   "--geneSearchStrategy", "EVAL",
                                   "-p", gr_dir + "/GeneRax"])
            print(ans)
            if ans == 0:
                worked = 1
                break
            i += 1
        if not worked:
            print("Generax didn't work for some reason - the temporary files are still there so go check them out.")
            print("The directory is " + gr_dir)
            sys.exit()

        #It bloomin works!
        #then just go through the lines in the per-species-event-counts file
        #and count the number of SL=., saving each to a dictionary
        for line in open(gr_dir + "/GeneRax/reconciliations/temp_eventCounts.txt", "r").readlines():
            if "SL:" in line:
                results[leaf.name] = int(line.split(":")[1])
                break
        
        #Then just delete the files and crack on
        shutil.rmtree(gr_dir + "/GeneRax")
        print(results)
    
    
    #calculate the reduction in number of losses
    diff_table = {}
    proportion_table = {}
    orig = results["original"]
    for key, value in results.items():
        if key != "original":
            diff_table[key] = orig - value
            proportion_table[key] = round((orig - value) / orig, 3)


    #write results
    with open(output + "_loss_reduction.csv", "w", encoding = "utf8") as out:
        out.write("gene,loss_reduction,proportion_reduction\n")
        for key, value in diff_table.items():
            out.write(key + "," + str(value) + "," + str(proportion_table[key]) + "\n")
    
    #remove the temporary directory
    shutil.rmtree(gr_dir)
    print("job done")


if __name__ == "__main__":
    main()
