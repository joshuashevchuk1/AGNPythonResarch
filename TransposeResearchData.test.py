import pandas as pd


def transposeData():
    df = pd.read_csv('file.csv')
    transposedData = df.T
    return transposedData


def saveData(data):
    data.to_csv(str(data_frame[n]['DirName'] + "_" + "time_GTN" + ".csv"), sep='\t')

saveData(transposeData())

