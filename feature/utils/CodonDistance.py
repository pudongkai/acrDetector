# coding=utf-8
import re
import numpy as np
from Bio import SeqIO


def distance(gbk):
    # caculate Acr distance
    handle = SeqIO.parse(gbk, "genbank")
    gene_dict = {}
    Dict_genome_CDS = {}
    genome_cds_seq = ''
    for record in handle:
        for feature in record.features:
            if feature.type == 'CDS' and feature.qualifiers.get('protein_id') != None:
                Dict_genome_CDS[feature.qualifiers.get('protein_id')[0]] = calculate_frequence(
                    str(feature.extract(record.seq)))
                genome_cds_seq = genome_cds_seq + str(feature.extract(record.seq))
                gene_dict[feature.qualifiers.get('protein_id')[0]] = str(feature.extract(record.seq))
    Dict_genome_CDS['genome'] = calculate_frequence(genome_cds_seq)
    return gene_dict, Dict_genome_CDS


def calculate_frequence(target_seq):
    #################################################
    # code by ZengZhi
    #################################################
    try:
        code_l = re.findall(r'(\w\w\w)', target_seq)
    except Exception as e:
        # print type(target_seq)
        print(e)
    sum_num = len(code_l)
    CodonsDict = dict([('TTT', 0), ('TTC', 0), ('TTA', 0), ('TTG', 0), ('CTT', 0), ('CTC', 0), ('CTA', 0),
                     ('CTG', 0), ('ATT', 0), ('ATC', 0),
                     ('ATA', 0), ('ATG', 0), ('GTT', 0), ('GTC', 0), ('GTA', 0), ('GTG', 0), ('TAT', 0),
                     ('TAC', 0), ('TAA', 0), ('TAG', 0),
                     ('CAT', 0), ('CAC', 0), ('CAA', 0), ('CAG', 0), ('AAT', 0), ('AAC', 0), ('AAA', 0),
                     ('AAG', 0), ('GAT', 0), ('GAC', 0),
                     ('GAA', 0), ('GAG', 0), ('TCT', 0), ('TCC', 0), ('TCA', 0), ('TCG', 0), ('CCT', 0),
                     ('CCC', 0), ('CCA', 0), ('CCG', 0),
                     ('ACT', 0), ('ACC', 0), ('ACA', 0), ('ACG', 0), ('GCT', 0), ('GCC', 0), ('GCA', 0),
                     ('GCG', 0), ('TGT', 0), ('TGC', 0),
                     ('TGA', 0), ('TGG', 0), ('CGT', 0), ('CGC', 0), ('CGA', 0), ('CGG', 0), ('AGT', 0),
                     ('AGC', 0), ('AGA', 0), ('AGG', 0),
                     ('GGT', 0), ('GGC', 0), ('GGA', 0), ('GGG', 0)])
    for item in CodonsDict.keys():
        if item in code_l:
            CodonsDict[item] = code_l.count(item)
    tensor = []
    for i in sorted(CodonsDict.keys()):
        tensor.append(float('%.4f' % (CodonsDict[i]/float(sum_num))))
    return np.array(tensor)
