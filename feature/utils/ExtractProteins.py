# coding=utf-8
from Bio import SeqIO


def extract(gbkFile):
    faa = gbkFile.replace(".gbff", ".faa")

    proteins = []
    fo = open(faa, "w")

    handle = SeqIO.parse(gbkFile, 'genbank')
    for record in handle:
        LOCUS = record.id
        # record.annotations
        for (index, feature) in enumerate(record.features):
            if feature.type == 'gene' or feature.type == 'CDS':
                qualifier = feature.qualifiers
                location = str(feature.location)

                if feature.type == 'CDS' and 'translation' in qualifier:
                    qualifier = feature.qualifiers
                    protein_id = qualifier['protein_id'][0]
                    translation = qualifier['translation'][0]
                    product = qualifier['product'][0]
                    # print protein_id
                    fo.write('>%s|%s\t%s\t%s\n%s\n' % (protein_id, LOCUS, location, product, translation))
                    proteins.append(protein_id)
    fo.close()
    return faa, proteins
