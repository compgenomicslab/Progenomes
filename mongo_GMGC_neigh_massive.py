from __future__ import division
import json
from pymongo import MongoClient
from collections import Counter
import os,sys
import time
import re

current_directory = os.path.dirname(os.path.abspath(__file__))
save_folder = current_directory+"/"
output = open(save_folder+"neigh_results.txt","w")



# Llamamos a la base de datos de Mongo
db = None
if not db:
	client = MongoClient()
	client = MongoClient('localhost', 27017)
	db = client.freeze11
	coll_members = db.members
	coll_contig_clusters = db.contig_clusters
	coll_trees = db.trees
	coll_contigs = db.contigs # now include correct orf
	coll_gene_description = db.gene_description
	coll_annotations = db.annotations
## Entrada de datos en el script






### Parseamos el fichero con las descripciones de los KEGGs pathways
### y generamos un diccionario para almacenarlos.
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




def retrieve_kegg(gene,taxID):
	""" obtiene el kegg y el eggnog de una lista de orfs de la coleccion members """
	
	selected = []
	kegg = "NA"
	for n in gene:
		freeze11_function = coll_members.find({"s":n,"l":taxID})

		try:
			for element in freeze11_function:
				kegg = element["k"]
		except:
			kegg = "NA"

		selected.append(kegg)
	return selected


def get_strand(list_of_genes):
	"""busca el strand para la lista de genes neighbours y el query"""
	list_of_strand = []
	for gene in list_of_genes:
		if gene != "NA":
			sequence_function = coll_contigs.find({"o":gene})
			for element in sequence_function:
				strand = element["str"]
				list_of_strand.append(strand)



	return list_of_strand

def retrieve_neighbours_data(sequences):
	""" obtiene una lista con los cuatro genes -2-1+1+2 que rodean en el contig al orf que estamos analizando, si la lista 
	de los cuatro genes contiene algun gen missing no metemos dicha lista de genes en el analisis"""
	neigh_dict = {}


	for sequence in sequences.split(","):
		if "_" in sequence:
			fields = sequence.split("_")
			contig = fields[0]
			query_gene = fields[1]
			
			try:
				sequence_function = coll_contig_clusters.find({"contig":contig})

				for element in sequence_function:

					seqs = element['o']

				gene_list =[-2,-1,1,2]

				gene_ordered = []


				for gene_pos in gene_list:	
					gene = query_gene
					try:
						prefix,number = re.search('([^0-9]+)(\d+)$', gene).groups() # busca letras y numernos en el gen y si los hay los separa en letra por un lado y numero por otro
						gene = int(number)+int(gene_pos)
						#print gene
						gene = prefix+str(gene)
						#print contig,gene
					except:
						gene = int(gene)+int(gene_pos)

					if gene in seqs:
						gene = contig+"_"+gene
						gene_ordered.append(gene)
					else:
						gene_ordered.append("NIcontig") # no esta presente en el contig

				if 	"NIcontig" not in gene_ordered:	# Asi si alguno de los genes no tiene neighbourhood lo descartamos
					neigh_dict[query_gene]=gene_ordered
			except:
				continue

	return neigh_dict


def retrieve_REFSEQ_description(seq):
	orf_features = coll_gene_description.find({"s":seq})
	for feature in orf_features:
		description = feature['d']

	return description.replace(",",".")



def function_for_drawn(keggs_for_draw,OG,strand_for_draw,gene_list):
	"""con esta funcion podemos plotear el OG que queramos en un grafico"""

	keggs_for_draw = [keggs_for_draw[0],keggs_for_draw[1],OG,keggs_for_draw[2],keggs_for_draw[3]]

	kegg_not_allowed = ["01100","01110","01120","01130","01200","01210","01212","01230","01220"]
	count_for_draw = 0
	position = 0
	for gene in gene_list:
		description = retrieve_REFSEQ_description(gene) 
		OG = keggs_for_draw[int(position)]
		strand = strand_for_draw[int(position)]
		for kegg in kegg_not_allowed:
			OG = OG.replace(kegg+",","").replace(","+kegg,"")
		OG = OG.replace(",","/")
		print("1,"+gene+","+"100,200,"+strand+","+OG+","+description)
		position +=1


def create_strand_for_drawn(neighbours_list,query_gene):
	"""generamos una lista con los strands del gen query y sus cuatro neighbours"""

	strand_for_draw = []
	gene_list = ",".join(neighbours_list).split(",")
	contig = gene_list[0].split("_")[0]
	gene_list = [gene_list[0],gene_list[1],contig+"_"+query_gene,gene_list[2],gene_list[3]]
	strands = get_strand(gene_list)
	for orientation in strands:
		strand_for_draw.append(orientation)

	return 	[strand_for_draw,gene_list]





def main(file,taxID): #entra el archivo con los OG Kegg seqID 
	for line in file:
		kegg_list = []
		query_kegg_filtered = []
		kegg_not_allowed = ["01100","01110","01120","01130","01200","01210","01212","01230","01220"]

		gmgc_orf_dict = {}
		number_neigh = 0
		neigh_with_keggs = 0
		analysed_orfs = 0
		keggs_for_draw_list=[] #obtener los stats de los genes con Kegss
		unigenes_functions=[]
		neigh_with_keggs_real = 0

		fields = line.split("\t")
		OG = fields[0]

		sequences = fields[2]
		analysed_orfs = len(sequences.split(','))

		kegg_query = fields[1] 

		# hacemos esto para descartar los keggs que pertenecen al grupo 1.0 de keggs pathways que son generales y no estan el kegg dict	
		for kegg in kegg_query.split(","):
			if kegg not in kegg_not_allowed:
				query_kegg_filtered.append(kegg)

		kegg_query = ",".join(query_kegg_filtered)


		gene_ordered = retrieve_neighbours_data(sequences) #diccionario key gene0 y el resto son -2-1+1+2

		for k,v in gene_ordered.items():
			neighbours_list = v	
			query_gene = k
			neighbours =retrieve_kegg(neighbours_list,taxID) # lista que contiene los kegggs de cada neighbourhood

			################################ modulo para generar un output para graficar un eggnog concreto
			# strand_for_draw = []
			# gene_list = ",".join(neighbours_list).split(",")
			# contig = gene_list[0].split("_")[0]
			# gene_list = [gene_list[0],gene_list[1],contig+"_"+query_gene,gene_list[2],gene_list[3]]
			# strands = get_strand(gene_list)
			# for orientation in strands:
			# 	strand_for_draw.append(orientation)

			#print strand_for_draw
			################################


			keggs_for_draw = []


			for kegg in neighbours:
				#print kegg
				number_neigh +=1 #sumamos 1 por cada neighbourhood con kegg

				if kegg != "NA":
					neigh_with_keggs_real +=1 # numero real de genes neigh que tienen un kegg y no una asignacion NA

				keggs_for_draw.append(kegg)


				for sub_kegg in kegg.split(","): #depuro la lista de kegg para quedarme solo con los kegg pathways

					if sub_kegg in kegg_dict :
						unigenes_functions.append(sub_kegg)




			################ print keggs to use in graphic script
			#strand_for_draw = create_strand_for_drawn(neighbours_list,query_gene)[0]
			#gene_list = create_strand_for_drawn(neighbours_list,query_gene)[1]
			#function_for_drawn(keggs_for_draw,OG,strand_for_draw,gene_list)


					#for neigh visulization
			keggs_for_draw_list.append(keggs_for_draw)

			kegg_depured_list =[]
			for k in keggs_for_draw_list:
				if k.count("NA") != 4:
					kegg_depured_list.append(k)


			neigh_orf_with_keggs = len(kegg_depured_list)
			#print kegg_depured_list




		query_kegg_list=kegg_query.split(",")


		kegg_description_dict = {}

		Count = Counter(unigenes_functions)
		Count = {key:val for key, val in Count.items() if val != 1}		## Eliminamos los Keggs que estan una sola vez
		for k,v in Count.items():
			count = int(v)
			neigh_with_keggs += count # numero total de keggs en el analisis menos los keggs que han aparecido una sola vez


		#print neigh_with_keggs
		subject_kegg_list = []
		subject_kegg_dict = {} #almacenamos todos los kegg que pasen el porcentaje para compararlos con los query_kegg_list
		for k,v in Count.items(): # v es el numero de veces que ha aparecido en ese cluster por vencidad dicha funcion
			KEGG_count = int(v)
			percentage = float(v/neigh_with_keggs)*100 # calculamos el porcentaje sobre el numero de genes neighbours with keggs
			percentage = ("{0:.2f}".format(percentage))

			if float(percentage) >= 25:

				kegg = k
				kegg_description = get_kegg_description(kegg)
				kegg_description_dict[kegg]=kegg_description
				subject_kegg_dict[kegg]=percentage
				subject_kegg_list.append(kegg)


		hits_kegg_list = [] # lista que almacena los hits que coinciden con los iniciales del query
		hit_kegg_percentage = [] # guarda el porcentage con el que paso dicho kegg, porcentaje de abundancia entre los keggs.
		hits = 0
		hit_kegg_count_per_orf = {}
		hit_kegg_count_per_orf_negative = {}
		for k,v in subject_kegg_dict.items():
			kegg = k
			percentage = v
			k_p = k+"@"+str(percentage)
			hit_kegg_percentage.append(k_p)
			num_query_kegg = len(kegg_query.split(","))
			num_subject_kegg = len(subject_kegg_list)


			if kegg in query_kegg_list:
				hits +=1
				hits_kegg = kegg+"_"+str(percentage)
				hits_kegg_list.append(hits_kegg)


				kegg_count = 0
				for k in kegg_depured_list:
					for orf in k:
						if kegg in orf:
							kegg_count +=1
				#print kegg_count
				hit_kegg_count_per_orf[kegg]=kegg_count

			if kegg not in 	query_kegg_list:

				kegg_count = 0
				for k in kegg_depured_list:
					for orf in k:
						if kegg in orf:
							kegg_count +=1
				#print kegg_count
				hit_kegg_count_per_orf_negative[kegg]=kegg_count


		passed_kegg_average =0
		percentage_kegg_sum= 0
		passed_kegg_count = 0
		try:
			for kegg,count in hit_kegg_count_per_orf.iteritems():

				try:
					kegg_n = ("{0:.2f}".format(count/neigh_orf_with_keggs)) # proporcion del numero de orf que tenian el kegg en sus neighbourhood. cuanto mayor mejor
				except:
					kegg_n = 0

				percentage_kegg_sum = percentage_kegg_sum+float(kegg_n)

				passed_kegg_count += 1
			passed_kegg_average = percentage_kegg_sum/passed_kegg_count

		except:
			passed_kegg_average=0


		passed_kegg_average_negative =0
		percentage_kegg_sum= 0
		passed_kegg_count = 0
		try:
			for kegg,count in hit_kegg_count_per_orf_negative.iteritems():

				try:
					kegg_n = ("{0:.2f}".format(count/neigh_orf_with_keggs)) # proporcion del numero de orf que tenian el kegg en sus neighbourhood. cuanto mayor mejor
				except:
					kegg_n = 0
				percentage_kegg_sum = percentage_kegg_sum+float(kegg_n)
				passed_kegg_count += 1
			passed_kegg_average_negative = percentage_kegg_sum/passed_kegg_count

		except:
			passed_kegg_average_negative=0






		try:
			kegg_positives = ("{0:.2f}".format(int(hits)/num_query_kegg)) # cuanto mas cercano a 1 mejor.
		except:
			kegg_positives = 0

		try:
			kegg_accuracy = ("{0:.2f}".format(hits/len(subject_kegg_list))) # porcentaje de hits sobre todos los keggs predichos. Cuanto menor mejor
		except:
			kegg_accuracy = 0


		try:
			kegg_proportion = ("{0:.2f}".format(neigh_with_keggs/number_neigh)) # esta midiendo la dispersion de los keggs, si es alto quiere decir que la mayporia de los keggs estaban agrupados en cluster de mas de un kegg.
			#kegg_proportion = ("{0:.2f}".format(neigh_with_keggs_real/number_neigh)) # kegg totales que entraron en la lista de Counter dividido por el numero total de neigh en el analisis

		except:
			kegg_proportion = 0


		try:

			for kegg in query_kegg_list:
                                kegg_description = get_kegg_description(kegg)
                                kegg_description_dict[kegg]=kegg_description


			description_list = []
			for kegg, description in kegg_description_dict.iteritems():
				k_d= str(kegg)+"@"+str(description)
				description_list.append(k_d)


		except:
			description_list=[]




		# neigh_with_keggs == number of keggs that appears more than once in a OG analysis
		#anadir las descriciones al final y porcentajes de los kegg
		#print ("#gmgc"+"\t"+"query_keggs"+"\t"+"subject_keggs"+"\t"+"keggs_positives"+"\t"+"keggs_accuracy"+"\t"+"analysed_orfs"+"\t"+"neigh_genes"+"\t"+"neigh_with_keggs"+"\t"+"kegg_proportion"+"\t"+"presence_of_kegg"+"\t"+"hit_kegg_percentage"+"\t"+"kegg_description")
		output.write(
		
			str(OG)+"\t"+
			",".join(query_kegg_list)+"\t"+
			",".join(subject_kegg_list)+"\t"+
			str(kegg_positives)+"\t"+
			str(kegg_accuracy)+"\t"+
			str(analysed_orfs)+"\t"+
			str(number_neigh)+"\t"+
			str(neigh_with_keggs)+"\t"+
			str(kegg_proportion)+"\t"+
			str(passed_kegg_average)+"\t"+
			str(passed_kegg_average_negative)+"\t"+
			",".join(hit_kegg_percentage)+"\t"+
			";".join(description_list)+"\n"
			)








make_kegg_dict (kegg_pathways)

input_file = sys.argv[1]
taxID = input_file.split(".")[0].split("_")[0] #taxID_level
file = open(input_file, "r")

output_folder= current_directory+"/output_analysis/"
output = open(output_folder+taxID+".txt","w")

print(taxID)


main(file,taxID)
