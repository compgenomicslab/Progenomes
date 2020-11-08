import numpy as np
import pandas as pd 
import sys



global output_list
output_list = []


def output_list_of_images(df):
    """obtenemos los genes que estan en hebra negativa """

    gene_list = list(df.loc[:,"Gene_id"])
    strand_list = list(df.loc[:,"Strand"])
    for n in range(2,len(gene_list),5):
        gene = gene_list[n].replace(".","_")
        strand = strand_list[n]
        if strand == "-":
        
        	output_list.append(n)


      


def inverse_operon_for_negative_strand (df):
	""" si la ORF es negativa cambia la orientacion de su operon para graficar 
	siempre las ORF en sentido positivo """	


	df1 = df
	for n in output_list: # para cada posicion en el index q contiene un ORF en strand -
		
		df2 = df.reindex([n+2,n+1,n,n-1,n-2])

		for i in [n+2,n+1,n,n-1,n-2]:
			strand = df2.loc[i,"Strand"]
			if strand == "+":
				df2.loc[i,"Strand"]="-"
			else:
				df2.loc[i,"Strand"]="+"	

		#print(df.index[[n+2,n+1,n,n-1,n-2]])
		df3 = df1.drop(df.index[[n+2,n+1,n,n-1,n-2]])	

		frames = [df2, df3]
	
		df_final = pd.concat(frames) # unimos df2 con df3
		df1 = df_final


	return df1




#input = sys.argv[1] # csv with neighbourhood


def order_gene_list_by_strand(gene_coordinates_file):
	df = pd.read_csv(gene_coordinates_file) 
	output_list_of_images(df) 
	df = inverse_operon_for_negative_strand(df)
	df.to_csv("tmp/genes_coordinates_sorted.csv",index=False)





