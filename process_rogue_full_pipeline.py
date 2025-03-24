#!/usr/bin/env python3.6
"""
for a given gene family
    i.  make dummy version with terminal duplication collapsed
    ii. run remove_rogue
    iii. remove the highest scoring gene according to remove_rogue.py if that
        is greater than say 5
    iv. jump to ii until none score over 5
    v. un_collapse and add to reformatted, then I can do the EVAL generax
"""

import ete3
import os
import subprocess
import sys
from random import randint
import my_module as mod

def get_args():
    """
    Get user arguments.
    family_alignment
    family_tree
    threshold
    prefix
    """
    if not len(sys.argv) == 5:
        print("Usage: python3 process_rogue_full_pipeline.py family_alignment family_tree loss_threshold prefix")
        sys.exit()
    return sys.argv[1:]

def collapse_same_sp(tree):
    """
    Collapse terminal pairs of tips belonging to the same species.

    Keep doing this until no more remain

    Return the tree and the sequences removed with the one that represents
    them as a dictionary.
    """
    replacements = {}
    to_remove = []
    leaves = []
    for leaf in tree:
        leaves.append(leaf.name)
    unfinished = 1
    while unfinished:
        done = 0
        for node in tree.traverse():
            if len(node) == 2:#this assumes bifurcation
                #then both descendants are terminal tips
                sp1 = node.children[0].name.split("_")[0]
                sp2 = node.children[1].name.split("_")[0]
                if sp1 == sp2:
                    #then they are the same species and should be collapsed
                    if node.children[1].name in replacements:
                        replacements[node.children[0].name] = node.children[1].name + ";" + replacements[node.children[1].name]
                    else:
                        replacements[node.children[0].name] = node.children[1].name
                    to_remove.append(node.children[1].name)
                    node.remove_child(node.children[1])
                    node.delete()
                    done = 1
                    break
        if not done:
           unfinished = 0
    return tree, replacements, to_remove

def main():
    """Do the things."""
    alignment, treefile, threshold, prefix = get_args()
    threshold = int(threshold)
    unique = str(randint(0,9))
    os.makedirs(prefix, exist_ok=True)
    
    #here, I still need to collapse single species clades first.
    tree = ete3.Tree(treefile, format = 1)
    tree, id_table, to_remove = collapse_same_sp(tree)
    tree.write(outfile = prefix + "/tree.nwk", format = 1)
    seqs = mod.read_fasta(alignment)
    with open(prefix + "/seqs.fa", "w", encoding = "utf8") as out:
        for key, value in seqs.items():
            if key not in to_remove:
                out.write(">" + key + "\n" + value + "\n")
    cwd = os.getcwd()
    done = 0
    #temp for debugging
    removed = []
    while not done:
        #run rogue removal
        #remove the sequence with the highest number of extra losses above the threshold
        #if there isn't 1, set done to 1
        #python3 remove_rogue.py alignment gene_tree species_tree output_prefix
        ans = subprocess.call(["python3", "../../../../remove_rogue/remove_rogue.py",
                               "seqs.fa", "tree.nwk",
                               "../../no_gum.nwk", "temp"], cwd=cwd + "/" + prefix)
        print(ans)
        if ans:
            print("remove_rogue has stopped at some point - go check it out")
            sys.exit()

        #read results
        current = 0
        current_line = ""
        to_remove = ""
        for line in mod.get_file_data(cwd + "/" + prefix + "/temp_loss_reduction.csv")[1:]:
            fields = line.split(",")
            if int(fields[1]) > current and int(fields[1]) >= threshold:
                to_remove = fields[0]
                current = int(fields[1])
                current_line = line
        if len(to_remove) == 0:
            done = 1
            break
        removed.append(to_remove)
        if to_remove in id_table:
            if ";" in id_table[to_remove]:
                removed.extend(";".split(id_table[to_remove]))
            else:
                removed.append(id_table[to_remove])

        #now remove the sequence in question
        seqs = mod.read_fasta(cwd + "/" + prefix + "/seqs.fa")
        with open(cwd + "/" + prefix + "/seqs.fa", "w", encoding = "utf8") as out:
            for key, value in seqs.items():
                if key != to_remove:
                    out.write(">" + key + "\n" + value + "\n")

        t = ete3.Tree(cwd + "/" + prefix + "/tree.nwk", format = 1)
        to_prune = []
        for leaf in t:
            if to_remove not in leaf.name:
                to_prune.append(leaf.name)
        t.prune(to_prune)    
        t.write(format = 1, outfile = cwd + "/" + prefix + "/tree.nwk")

    with open(cwd + "/" + prefix + "/removed_sequences.txt", "w", encoding = "utf8") as out:
        out.write("\n".join(removed) + "\n")

    #now I just need the full pruned tree with the terminal duplications
    #included and the alignment with all those seqs too.
    tree = ete3.Tree(treefile, format = 1)
    to_prune = []
    for leaf in tree:
        if leaf.name not in removed:
            to_prune.append(leaf.name)
    tree.prune(to_prune)
    tree.write(format = 1, outfile = prefix + "/pruned_tree.nwk")
    og_ali = mod.read_fasta(alignment)
    with open(prefix + "/pruned_alignment.fa", "w", encoding = "utf8") as out:
        for key, value in og_ali.items():
            if key not in removed:
                out.write(">" + key + "\n" + value + "\n")


if __name__ == "__main__":
    main()
