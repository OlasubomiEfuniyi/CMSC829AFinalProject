# CMSC829AFinalProject
This repository contains the data set and scripts I used for my CMSC829A final project. Checkout the final paper named "CMSC829A_Final_Paper.pdf"

* **{exon/intron}_fasta_files** contains selected exon fasta files.
* **{exon/intron}_fastme_files** contains the files produces by running FastME on the corresponding data set.
* **{exon/intron}_iqtree_fastme_files** contains the files produced by running IQ-TREE using the FastME tree as the start tree.
* **{exon/intron}_iqtree_mp_files** contains the files produced by running IQ-TREE using the default Maximum Parsimony tree as the start tree.
* **{exon/intron}_phylip_files** contains the phylip version of the files in {exon/intron}_fasta_files.
* **exon_mp_files** contains the files produced by running MPBoot on the exon data set.
* **exon_iqtree_computed_mp_files** contains the files produced by running IQ-TREE using the MPBoot tree as the start tree.
* **images** contains the pictures of generated trees.
* **analysis.py** is a script used to calculate tree reliability and how many bipartitions two trees have in common.
* **gap_counter.py** is a script used to count the number of gaps in each sequence in a Multiple Sequence Alignment.
* **max_length_sequences_finder.py** is a script used to find the file that provides the longest MSA aligned sequences from the downloaded data. This is what is saved in {exon/intron}_fasta_files.

# Dependencies
*analysis.py* depends on dendropy.
