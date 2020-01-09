import joblib
import os
import sys

from predict.utils import FeatureHandler
from predict.utils import ResultSaver

model = "predict/model/train_model.m"


def predict(dataSetFile):
    dataSetFile = dataSetFile + ".txt"

    print("Predicting Acr...")

    X, protein_ids = FeatureHandler.createDataset(dataSetFile)

    dtc = joblib.load(model)
    y_predict = dtc.predict(X)

    # get potential Acr Accession
    predict_proteins = set()
    for i in range(0, len(y_predict)):
        p = y_predict[i]  # predict value
        if p == '1':
            predict_proteins.add(protein_ids[i])

    # save into file
    ResultSaver.save(dataSetFile, predict_proteins)

    return predict_proteins







if __name__ == "__main__":
    _, dataSetFile = sys.argv
    # dataSetFile = "/home/pudongkai/Downloads/e.coli/CP014225.1.txt"
    predict(dataSetFile)

