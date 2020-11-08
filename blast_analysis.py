import os, sys
from Bio.Blast.Applications import NcbiblastnCommandline



def blast(queryFile,blastdb):



    blastn_cline = NcbiblastnCommandline(cmd='blastp',
                                        query=queryFile,
                                        db=blastdb,
                                        outfmt='"6 qseqid sseqid pident evalue bitscore qcovs sseq"',
                                        evalue='0.00001',
                                        num_alignments = '5000',
                                        num_threads='12'
                                        )	

    stdout, stderr = blastn_cline()
    result = stdout.split("\n")
    return stdout