from pymongo import MongoClient
from mongo_client import *



def get_eggnog_tree(eggnog,path_to_save):
    
    output_file = open(path_to_save,"w")

    get_tree_information = coll_trees.find({"e":eggnog})        
    for element in get_tree_information:
        tree = element["tree"]
    
    output_file.write(tree)


#eggnog = "41F15"
#get_eggnog_tree(eggnog,path_to_save)
