import sys
from feature import CalculateFeatures
from predict import PredictAcr


def main(fileName):
    CalculateFeatures.CalculateFeatures(fileName)
    PredictAcr.predict(fileName)


if __name__=="__main__":
    _, fileName = sys.argv
    main(fileName)
