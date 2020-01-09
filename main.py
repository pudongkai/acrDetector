import sys
from feature import CalculateFeatures
from predict import PredictAcr


def main(fileName):
    print("Calculating features...")
    CalculateFeatures.CalculateFeatures(fileName)
    print("Predicting Acr...")
    PredictAcr.predict(fileName)


if __name__=="__main__":
    _, fileName = sys.argv
    main(fileName)
