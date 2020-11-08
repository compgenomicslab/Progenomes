# Progenomes
Python package for Jorge´s JS visualization


usage:

Usage: main.py [mode]

Modes:
-B for blast analysis using a query sequence

  usage: main.py -B <protein_sequence.fasta> <identity_cut-off>


-E to retrieve all sequence from query Eggnog

  usage: main.py -E <Eggnog>


For blast analysis its neccesary to get acces to Progenomes blast dabase at: 
/home/giner/Progenomes_fr11/blast_e5DB

File outputs will be save at /tmp folder:

genes_coordinataes.csv (neighbourhood gene coordinates to plot in JS, not sorted respect to center gene in positive orientation)
genes_coordinataes_sorted.csv (neighbourhood gene coordinates sorted respect to center gene in positive orientation)
protein_sequence.fasta (query fasta protein of blast)
result_blast.tsv (result of blast in table format, HEADER: query_sequence_header;subject_sequence_name;identity;e-value;bit score;coverage;subject seq protein)
tree.nwx (EGGnog tree in newick format)
