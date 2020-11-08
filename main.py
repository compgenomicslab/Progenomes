
import sys

from mongo_client import *
from filter_blast import *
from blast_analysis import *

from retrieve_neighbours_data import *
from order_gene_list_by_strand import *
from getting_trees import *



#blast data base containing progenomes representatives protein sequences.
#add here correct path to blast db.
blastdb = "~/data2/Eggnog5_taxID_members/eggnog5_members_annotations_trees/blast_e5DB/e5.proteomes"


def analysis(protein_sequence,identity_cutoff,blastdb = blastdb):

    """
        Esta seccion genera una lista con aquellos eggnos que contienen secuencias con genes ortologos a la secuencia introducida
        hasta el cut-off de identidad indicado.
        datos mostrados para cada eggnog:
        -eggnog
        -numero de secuencias con las que tienen un alto nivel de homologia en el eggnog
        -numero total de secuencias del eggnog
        -identidad media que tiene la secuencia query con aquellas secuencias con la que muestra homologia en el eggnog
        -nombre de los genes que codifican para estas proteinas
        -KEGG del eggnog
        -Descripcion del eggnog
    """

    #remove temporary files
    tmp_folder = "tmp/"
    filelist = [ f for f in os.listdir(tmp_folder) if f.endswith(".csv") or f.endswith(".nwx")]
    for f in filelist:
        os.remove(os.path.join(tmp_folder, f))

    protein_sequence_file = open(protein_sequence,"r")
    protein_sequence_file = protein_sequence_file.read()
    print(protein_sequence_file)
    protein_file =open("tmp/protein_sequence.fasta","w")
    protein_file.write(protein_sequence_file)
    protein_file.close()
    
    
    protein_name = protein_sequence.split("\n")[0].replace(">","")

    #hacemos blast con la protei_NA_ query
    result_blast = blast("tmp/protein_sequence.fasta",blastdb)
    

    result_blast_file = open("tmp/result_blast.tsv","w")
    result_blast_file.write(result_blast)
    
    # #a_NA_lizamos los resultados con 
    # #parameters
    #identity_cutoff = 90.0

    target_genes_list = retrieve_OGs("tmp/result_blast.tsv",identity_cutoff)
    genes_dict_processed, eggnog_list, taxID_list= retrieve_gene_information(target_genes_list,coll_contigs,coll_annotations)
    result = printing_results(eggnog_list, genes_dict_processed, coll_members)




    return result




def make_tree(eggnog):
    
    #remove temporary files
    tmp_folder = "tmp/"
    filelist = [ f for f in os.listdir(tmp_folder) if f.endswith(".csv") or f.endswith(".nwx")]
    for f in filelist:
        os.remove(os.path.join(tmp_folder, f))


    try:
        eggnog = eggnog.split("@")[0]
    except:
        print("running!")    
    #retrieve corrdinates of neighbours genes from a eggnog 
    generate_genes_table_for_drawn(eggnog,"tmp/genes_coordinates.csv")
    
    #sort genes by strand to draw all query genes by positive strand
    table = order_gene_list_by_strand("tmp/genes_coordinates.csv")

    #get the eggnog tree
    output_file_for_tree = "tmp/tree.nwx"
    get_eggnog_tree(eggnog,output_file_for_tree)






#Si quieres hacer busquedas mediante blast y luego seleccionar los eggnogs que contienen homologos de esa proteina
# utiliza la opcion -B, si quieres buscar y graficar todos los integrantes de un EGGnog us la opcion -E


if len(sys.argv) > 1:
    if sys.argv[1] == "-B":


        protein_sequence = sys.argv[2]  #archivo fasta que contiene la proteina query
        identity_cutoff = sys.argv[3]  # cut-off de identidad para seleccionar los resultados del BLAST
        result = analysis(protein_sequence,identity_cutoff)

        eggnog = result[0][0] #por defecto de la lista que obtiene de eggnogs (result) se queda con el primero.
        make_tree(eggnog)

    elif sys.argv[1] == "-E":

        eggnog = sys.argv[2]
        make_tree(eggnog)

    else:
        print("Usage: main.py mode\n")
        print("-B for blast analysis using a query sequence\n")
        print("usage: main.py -B <protein_sequence.fasta> <identity_cut-off>\n")
        print("-E to retrieve all sequence from query Eggnog\n")
        print("usage: main.py -E <Eggnog (containing  TaxID@eggnog)>\n")

else:        

    print("Usage: main.py mode\n")
    print("-B for blast analysis using a query sequence\n")
    print("usage: main.py -B <protein_sequence.fasta> <identity_cut-off>\n")
    print("-E to retrieve all sequence from query Eggnog\n")
    print("usage: main.py -E <Eggnog>\n")
