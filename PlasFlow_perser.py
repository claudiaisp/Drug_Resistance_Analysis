import glob
from difflib import SequenceMatcher

def plasflow_reader(file):
    f = open (file, 'r')
    lines = f.readlines()
    header = lines[0].split('\t')
    static_array = []
    del lines[0]
    for line in lines:
        if float(line.split('\t')[2].split('_')[5]) > 80 and float(line.split('\t')[3]) > 500:
            if line.split('\t')[5].find("unclassified") != -1 :
                value = -1
                array = [line.split('\t')[2], line.split('\t')[3], line.split('\t')[5], value]
                static_array.append(array)
            else:
                pos = header.index(line.split('\t')[5])
                array = [line.split('\t')[2], line.split('\t')[3], line.split('\t')[5], line.split('\t')[pos]]
                static_array.append(array)
        else:
            pass
    return static_array

def plasflow_writer(array):
    f = open("plasmids_filtered.txt", "a")
    for line in array:
        f.write('{0},{1},{2},{3},{4}'.format(line[0], line[1], line[2], line[3], '\n'))





if __name__ == '__main__':
    prediction = glob.glob('/*.txt')
    prediction.sort()
    #print(prediction)
    file = "test.plasflow_predictions.txt"
    for file in prediction:
        plasflow_writer(plasflow_reader(file))