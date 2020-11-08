from __future__ import print_function
import sys

from collections import Counter
from mongo_client import *



eggnog_list = []



def retrieve_gene_information(target_genes_list,coll_contigs,coll_annotations):
    """
        output files:
        genes_dict_processed: diccio_NA_rio con los genes y su informacion
        eggnog_list: lista con todos los eggnogs que tienen estos genes asociados
        taxID_list: lista con todos los taxID asociados a los eggnogs

    """
    
    #primero obtenemos el strand y el codigo del gen o_i 
    #necesario para acceder a la base de datos annotation
    genes_dict_processed = {}
    eggnog_list = []
    taxID_list = []
    for target_gene in target_genes_list:        
        
        ident = target_gene.split("@")[1]
        target_gene = target_gene.split("@")[0]
        

        gene_information_from_contigs = coll_contigs.find({"o":target_gene})
        for element in gene_information_from_contigs:
            strand = element["str"]
            orf_ID_for_annotation = element["o_i"]
            

        
            gene_information_from_annotation = coll_annotations.find({"o":orf_ID_for_annotation})
            for item in gene_information_from_annotation:
                #kegg = item["Kegg"]
                gene_name = item["g_n"]
                description = item["d"]
                eggnog = item["enog"]
                kegg = item["Kegg"]

                
                
                for egg in eggnog.split(","):
                    taxID = egg.split("@")[1]
                    taxID_list.append(taxID)
                   
                    eggnog_list.append(egg)

                genes_information = "{}#{}#{}#{}#{}#{}#{}".format(target_gene,ident,strand,gene_name,eggnog,description,kegg)
                
                
                genes_dict_processed[target_gene] = genes_information

                

    return genes_dict_processed, eggnog_list, taxID_list




def calculate_sequences_with_higher_identity(blast_file):
    """
    procesa el archivo de blast generado por blast_a_NA_lysis.py
    y calcula el percentil 10 de lass secuencias con mas identidad
    """
    
    identity_list = []

    for line in blast_file:
        #query,target,ident,evalue,x,cov,seq = map(str, line.split("\t"))
        try:
            fields = line.split("\t")
            target = fields[1]
            ident = float(fields[2])

        except:
            continue
        
        identity_list.append(target)

    identity_list.sort(reverse = True)

    number_sequences = len(identity_list)

    ident_range_selected = 0.1*number_sequences

    percentil_10 = (round(ident_range_selected))

    return percentil_10




def retrieve_OGs(blast_file, identity_cutoff):
    blast_file = open(blast_file,"r")

    target_genes_list = []
    for line in blast_file:
        #query,target,ident,evalue,x,cov,seq = map(str, line.split("\t"))
        
        try:
            
            fields = line.split("\t")
            target = fields[1]
            ident = fields[2]

        except:
            continue
        
        if float(ident) > float(identity_cutoff):
            target_genes_list.append(target+"@"+str(ident))



    return target_genes_list


    




def printing_results(eggnog_list, genes_dict_processed,coll_members):
    
    result_list = []
    Count = Counter(eggnog_list)
    for k, v in Count.items():
        ident_list = []
        eggnog = k
        count = v
        taxID = int(eggnog.split("@")[1])

        if taxID != 1 and taxID != 2:
            OG_members= coll_members.find({"e":eggnog.split("@")[0]})
            
            sequence_count = 0
            for item in OG_members:

                if item == None:
                    continue
                else:   
                    #print (eggnog,item)
                    sequence = item["s"]
                    
                    sequence_count +=1

                    try: 
                        fields = genes_dict_processed[sequence].split("#")
                        ident = fields[1]
                        ident_list.append(float(ident))
                        desc = fields[5] 
                        kegg = fields[6]  
                        gene_name = fields[3]
                    
                        
                    except:
                        continue
                    
                ident_average = sum(ident_list)/len(ident_list)
                ident_average = "%0.2f" % (ident_average)

        else:
            continue

        result = [eggnog, count, sequence_count,ident_average,gene_name,kegg, desc]
        result_list.append(result)
    result_list = sorted(result_list, key=lambda x: x[3], reverse=True)
    return result_list





if __name__ == "__main__":
    blast_file = open(sys.argv[1],"r")
    target_genes_list = retrieve_OGs(blast_file)
    genes_dict_processed, eggnog_list, taxID_list= retrieve_gene_information(target_genes_list,coll_contigs,coll_annotations)


    printing_results(eggnog_list, genes_dict_processed, coll_members)
