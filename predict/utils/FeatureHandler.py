import pandas as pd
import numpy as np


def createDataset(file):
    row_data = {'Length': [],
                'Product': [],
                'Downstream-HTH': [],
                'MGE': [],
                'Codon': []}
    protein_ids = []

    with open(file) as f:
        f.readline()
        for line in f:
            temp = line.strip().split(',')  # accession, length, category, HTH_gap, island_type, codon_distance
            row_data['Length'].append(float(temp[1]))
            row_data['Product'].append(float(temp[2]))
            row_data['Downstream-HTH'].append(float(temp[3]))
            row_data['MGE'].append(float(temp[4]))
            row_data['Codon'].append(float(temp[5]))

            protein_ids.append(temp[0])

    return pd.DataFrame(row_data), protein_ids


def getFeatures(file, predict_proteins):
    island = {"1": "prophage", "2": "normal gene island", "0": "out of GIs"}

    features = {}
    with open(file) as f:
        f.readline()
        for line in f:
            accession, length, category, HTH_gap, island_type, codon_distance = line.strip().split(',')
            if accession not in predict_proteins:
                continue
            if accession not in features.keys():
                features[accession] = {}
            features[accession]["hth"] = HTH_gap
            features[accession]["mge"] = island[island_type]
            features[accession]["codonDistance"] = codon_distance

    return features


if __name__=='__main__':
    file = '/media/pudongkai/sda3/dataSet/predict_valid/old/CP011374.1/testDataset.arff'
    create_data(file)
