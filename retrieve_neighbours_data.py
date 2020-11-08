from __future__ import division, print_function
import json
from pymongo import MongoClient
from collections import Counter
from mongo_client import *
import os,sys
import time
import re

 
###Este script no analiza los eggnogs de nivel 1 y 2 porque son muy grandes y no se pueden graficar,
###para analizar dichos eggnogs hay que cambiar la opcion 2_group a "True" de la funcion 
###retrieve_eggnods_from_input_sequence



kegg_pathways = open("KEGGs_pathways.txt","r")


def make_kegg_dict(kegg_pathways):
    """generamos el diccionario con los kegg pathways y sus descripciones"""
    global kegg_dict
    kegg_dict = {}


    for line in kegg_pathways:
       fields = line.strip("\n").split("\t")
       kegg= fields[0]
       description= " ".join(fields[1::]).rstrip(" ")
       kegg_dict[kegg]= description

    return kegg_dict


def get_kegg_description(kegg):

        description = kegg_dict[kegg]
        return description    


def retrieve_OG_genes (eggnog, coll_members):
    
    """ 
    obtiene los genes que pertenecen a un determinado OG 
    """
    
    gene_list = []
    get_eggnog = coll_members.find({"e":eggnog})
    
    for seq in get_eggnog:
        eggnog_OG = seq["e"]
        eggnog_taxID = seq["l"]
        get_genes_from_eggnog = coll_members.find({"e":eggnog_OG})
        
        for gene in get_genes_from_eggnog:
            orf = gene["s"]
            if orf not in gene_list:
                gene_list.append(orf)


    return [gene_list, eggnog_OG, eggnog_taxID]



def retrieve_neighbours_data(genes, number_neigh, coll_contig_clusters):
    """ 
    retrieve neighbours genes for every gene of the OG
    number_neigh = number of neighbours genes to analyze
    """
    
    neigh_dict = {}

    gene_list = [i for i in range(-number_neigh,number_neigh+1,1)]

    seqs = ""
    for orf in genes:
        try:
            fields = orf.split("_")
            contig = fields[0]
            query_orf = fields[1]
        except:
            continue


        if "gene" in orf:        
            sequence_function = coll_contig_clusters.find({"contig":contig})

            for element in sequence_function:
                seqs = element['o']
                
       
        gene_ordered = []
        

        for gene_pos in gene_list:        
            try:
                prefix,number = re.search('([^0-9]+)(\d+)$', query_orf).groups() # busca letras y numernos en el gen y si los hay los separa en letra por un lado y numero por otro
                gene_number = int(number)+int(gene_pos)
                gene = prefix+str(gene_number)

            except:
                #asi evitamos que genes letra_num_letra como sr2032c rompan el ciclo. 

                try:
                    gene = int(query_orf)+int(gene_pos)
                except:
                    continue

            if "gene" in orf:
            
                if gene in seqs:
                    gene = contig+"_"+gene
                    gene_ordered.append(gene)
                else:
                    gene_ordered.append("NIcontig") # no esta presente en el contig
            else:
                length_query_orf = len(str(query_orf))
                length_gene = len(str(gene))
                if length_gene != length_query_orf:
                    gene = str(gene).zfill(length_query_orf)
                gene = contig+"_"+str(gene)
                gene_ordered.append(gene)


        #if "NIcontig" not in gene_ordered:        # Asi si alguno de los genes no tiene neighbourhood lo descartamos
            neigh_dict[orf]=gene_ordered

    return neigh_dict




def retrieve_gene_information_for_neighbours(neighbours_genes_dict, coll_contigs, coll_annotations):
    """
        obtenemos: 
        strand, 
        eggnog description, 
        keggs, 
        gene_name

    para el diccionario (neighbours_genes_dict) que contiene 
    los genes neighbourhood de todas las orfs de un OG
    """
    neighbours_genes_dict_processed = {}
    for gene, neighbours in neighbours_genes_dict.items():
  
    
    
    #primero obtenemos el strand y el codigo del gen o_i 
    #necesario para acceder a la base de datos annotation
        neighbour_list = []
        for orf in neighbours:
            gene_information = "_NA_"
            gene_information_from_contigs = coll_contigs.find({"o":orf}) #hay genes que no estan organizados en +-1 sino +-5 o +-10
            for element in gene_information_from_contigs:
      
                strand = element["str"]
                orf_ID_for_annotation = element["o_i"]


               
                gene_information_from_annotation = coll_annotations.find({"o":orf_ID_for_annotation})


                for item in gene_information_from_annotation:

                    eggnog = "/".join(item["enog"].split(","))
                    kegg = item["Kegg"]
                    if kegg == "":
                        kegg == "Not Kegg"
                    else:                   
                        kegg = "/".join(item["Kegg"].split(","))
                  
                    gene_name = item["g_n"]
                    description = item["d"].replace(","," ")

                    gene_information = "{}#{}#{}#{}#{}#{}".format(orf,strand,gene_name,kegg,eggnog,description)

                    neighbour_list.append(gene_information)

               
            if gene_information == "_NA_":            #neighbour_list.append("na") # hay genes que no estan en la base de datos de progenomes, no se por que?
                neighbour_list.append(gene_information)  


        neighbours_genes_dict_processed[gene]=neighbour_list
    return neighbours_genes_dict_processed



def retrieve_eggnods_from_input_sequence (gene,coll_contigs,coll_annotations, group_two = "False"):

    """
    Retrieve the list of eggnogs assgined to input gene sequence
    """
    
    eggnog_list = []

    gene_information_from_contigs = coll_contigs.find({"o":gene})
    for element in gene_information_from_contigs:
        strand = element["str"]
        orf_ID_for_annotation = element["o_i"]

        gene_information_from_annotation = coll_annotations.find({"o":orf_ID_for_annotation})
        for item in gene_information_from_annotation:
            eggnogs = item["enog"]
            for egg in eggnogs.split(","):
                fields = egg.split("@")
                egg = fields[0]
                taxID = str(fields[1])
                
                if group_two == "False":
                    if taxID != "1" and taxID != "2":  # por defecto no sacamos informacion de los eggnog nivel 2 ya que son muy grandes y no se pueden graficar
                        eggnog_list.append(egg)
                
                if group_two == "True":
                    if taxID != "1" :  # por defecto no sacamos informacion de los eggnog nivel 2 ya que son muy grandes y no se pueden graficar
                        eggnog_list.append(egg)


    return eggnog_list

def function_for_drawn(neighbours_genes_dict_processed,save_path):
        """con esta funcion podemos plotear el OG que queramos en un grafico"""

        output_file = open(save_path,"w")
        output_file.write("Unnamed: 0,Gene_id,Start,End,Strand,Kegg,EGGNOG,Function\n")
        for k,v in neighbours_genes_dict_processed.items():
            gene_from_eggnog = k
            list_of_neigh = v
     
            for neigh in list_of_neigh:
                if neigh == "_NA_":
                    output_file.write("1,"+"_NA_"+","+"100,200,"+"+"+","+"_NA_"+","+"_NA_"+","+"_NA_"+"\n") #Esto es o bien porque el gen no tiene datos en proveniente de progeneomes o el gen no sigue una serie +-1, sino +-5 o +-10
               
                else:
                 
                    gene,strand,gene_name,kegg,OG,OG_description = neigh.split("#")
                    output_file.write("1,"+gene+","+"100,200,"+strand+","+kegg+","+OG+","+OG_description+"\n")






"""
Future addition
In case we want to search using a gene name we can retrieve 
the eggnog using the retrieve_eggnods_from_input_sequence function.
Future addition
"""
#gene = "525903.Taci_0012"
#eggnog_list = retrieve_eggnods_from_input_sequence(gene,coll_contigs,coll_annotations)






def generate_genes_table_for_drawn(eggnog,save_file):


    

        
    gene_list, eggnog_OG, eggnog_taxID = retrieve_OG_genes(eggnog, coll_members)

    if len(gene_list) <= 200:

        neighbours_genes_dict = retrieve_neighbours_data(gene_list, 2,coll_contig_clusters)
        neighbours_genes_dict_processed = retrieve_gene_information_for_neighbours(neighbours_genes_dict, coll_contigs, coll_annotations)


        function_for_drawn(neighbours_genes_dict_processed,save_file)



#This is only for testing
eggnog = "41F15"
save_file="results.tsv"

generate_genes_table_for_drawn(eggnog,save_file)
