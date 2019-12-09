import os
from Bio import SeqIO
from . import FeatureHandler

def save(dataSetFile, predict_proteins):
    path = os.path.split(dataSetFile)[0]
    fo = open(path + "/result.txt", "w")

    features = FeatureHandler.getFeatures(dataSetFile, predict_proteins)

    gbkFile = dataSetFile.replace(".txt", ".gbff")
    handle = SeqIO.parse(gbkFile, 'genbank')
    for record in handle:
        locus = record.id
        genome = record.annotations["source"]
        fo.write("%s\n" % genome)
        fo.write("GenBank: %s\n" % locus)
        fo.write("Potential Acr(s): %d\n\n" % len(predict_proteins))
        fo.write("%s\n" % "id,protein_id,length,Aca_gap,codon_distance,mge,product")
        line_index = 1
        for (index, feature) in enumerate(record.features):
            if feature.type == 'gene' or feature.type == 'CDS':
                qualifier = feature.qualifiers
                if feature.type == 'CDS' and 'translation' in qualifier:
                    qualifier = feature.qualifiers
                    proteinId = qualifier['protein_id'][0]
                    product = qualifier['product'][0]
                    seq = qualifier['translation'][0]
                    length = len(seq)

                    if proteinId not in predict_proteins:
                        continue

                    mge = features[proteinId]["mge"]
                    hth = features[proteinId]["hth"]
                    codonDistance = features[proteinId]["codonDistance"]

                    data = [str(line_index), proteinId, str(length), hth, codonDistance, mge, product]
                    fo.write("%s\n" % ",".join(data))

                    line_index += 1

    fo.close()
