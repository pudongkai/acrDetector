# coding=utf-8
import json, re, os
from Bio import SeqIO


def PHASTERLocation(phaster, gbk):
    #########################################################
    # extract prophage location from PHASTER's result
    # return as a dict: {LOCUS: [[start, end], [start, end], ]}
    #########################################################
    locusNum = 0
    for record in SeqIO.parse(gbk, 'genbank'):
        locusNum += 1
        locus = record.id

    regex = re.compile('\\s+')
    phage_location = {}
    with open(phaster, 'r') as f:
        content = f.read().strip().replace('\n', r'\n').split(r'\n')
        if "error" in content:
            return phage_location
        for item in content:
            if '--------------' in item:
                index = content.index(item)
        for item in content[index + 1:]:
            temp = regex.split(item.strip())
            if len(temp) != 17:
                continue
            info = temp[4].split(',')
            location = info[-1].split(':')[-1]
            if locusNum > 1:
                locus = info[0]
                phage_location = addDict(phage_location, locus, list(map(int, location.split('-'))))
            else:
                phage_location = addDict(phage_location, locus, list(map(int, location.split('-'))))

    return phage_location


def islandProteins(island):
    ####################################################
    # Extract proteins on island from IslandViewer result
    # return proteins as a list
    ####################################################
    result = []
    if not os.path.isfile(island):
        print(island, 'is not exist')
        return []
    with open(island, 'r') as f:
        f.readline()
        for line in f:
            protein = line.split('\t')[4]
            result.append(protein)
    return result


def onPhage(Acr_location, phaster_location, product):
    ###############################################
    # Determine whether a Acr is on phage
    # return True or False
    #	 Acr_location:  location_start?:location_end
    #	 phage_location:  [location_start?-location_end]
    ###############################################
    Acr_start = min(map(int, Acr_location.split(':')))
    Acr_end = max(map(int, Acr_location.split(':')))
    for phage_location in phaster_location:
        phage_start = min(phage_location)
        phage_end = max(phage_location)
        if len(set(range(Acr_start, Acr_end)) & set(range(phage_start, phage_end))) > 10 or "phage" in product:
            return True
    return False


def identify(island, protein_tag, protein_location, phaster_location, product):
    #########################################################
    # input:
    # 	PHASTER and islandViewer result, protein id and protein location
    # detect island type: phage(1), island(2) or none(0)
    #########################################################

    # if Acr in prophage
    protein_of_islands = islandProteins(island)
    if onPhage(protein_location, phaster_location, product):
        # print protein_tag, ' is on phage.'
        island_type = 1
    # if Acr in island
    elif protein_tag in protein_of_islands:
        island_type = 2
    else:
        island_type = 0

    return island_type


def addDict(Dict, key, value):
    if key not in Dict.keys():
        Dict[key] = [value, ]
    else:
        Dict[key].append(value)
    return Dict
