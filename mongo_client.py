from pymongo import MongoClient
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
        coll_annotations = db.annotation