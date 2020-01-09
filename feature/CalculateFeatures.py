# coding=utf-8
import os
import sys
from Bio import SeqIO
import numpy as np

from feature.utils import ExtractProteins
from feature.utils import CodonDistance
from feature.utils import ClassifyProduct
from feature.utils import HTHGap
from feature.utils import MGE
from feature.utils import CodonDistance


def CalculateFeatures(fileName):
    gbkFile = fileName + ".gbff"
    phasterFile = fileName + ".phage"
    islandviewerFile = fileName + ".tsv"
    featureFile = fileName + ".txt"

    fo = open(featureFile, "w")
    fo.write("accession,length,product_category,HTH_gap,island_type,codon_distance\n")

    # calculate features
    faa, protein_id_list = ExtractProteins.extract(gbkFile)
    hth_proteins = HTHGap.extractHTH(faa)
    phaster_locations = {}
    if os.path.isfile(phasterFile):
        phaster_locations = MGE.PHASTERLocation(phasterFile, gbkFile)
    gene_dict, Dict_genome_CDS = CodonDistance.distance(gbkFile)  # calculate codon distance

    print("Calculating features...")
    for record in SeqIO.parse(faa, 'fasta'):
        des = record.description
        info, location, product = des.split('\t')
        accession, LOCUS = info.split('|')
        location = location.replace("<", "")
        location = location.replace(">", "")
        if "join" in location:
            location = location.split(",")[-1].strip()[:-1]  # join{[0:434](-), [6799346:6799785](-)}

        if LOCUS in phaster_locations.keys():
            phaster_location = phaster_locations[LOCUS]
        else:
            phaster_location = []

        length = len(record.seq)  # feature 1, length
        category = ClassifyProduct.classify(product)  # feature 2, product type
        HTH_gap = HTHGap.caculateHTHgap(protein_id_list, hth_proteins, accession)  # feature 3, Aca gap
        island_type = MGE.identify(islandviewerFile, accession, location[1:-4], phaster_location, product)  # feature 4, MGE type
        Acr_frequence = CodonDistance.calculate_frequence(gene_dict[accession])
        codon_distance = np.linalg.norm(Dict_genome_CDS['genome'] - Acr_frequence)  # feature 5, codon distance

        fo.write("%s\n" % ",".join([str(item) for item in
                                    [accession, length, category, HTH_gap, island_type, codon_distance]]))

    fo.close()


if __name__ == "__main__":
    _, fileName = sys.argv
    CalculateFeatures(fileName)
