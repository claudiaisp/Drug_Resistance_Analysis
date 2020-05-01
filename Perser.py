import glob
import numpy as np

def read_file (file):
    f = open(file, 'r')
    lines = f.readlines()
    auxiliary_array = []
    array = []
    del lines[0]
    for line in lines:
        array = [line.split(',')[0],line.split(',')[2].split(" ")[0].replace('"', ""), str(" ".join(line.split(',')[2].split(" ")[1:7]))]
        #print(array)
        if float(line.split(',')[0].split("_")[5]) > 80 and float(line.split(',')[0].split("_")[3]) > 500 and float(line.split(',')[-5]) > 95:
            auxiliary_array.append(array)


        else:
            pass
    f.close()
    return auxiliary_array

def write_file(auxiliary_array):
    #print(auxiliary_array)
    f = open("plasmids_filtered.txt", "a")
    for i in auxiliary_array:
        #str(auxiliary_array).replace("[","").replace("]","")
        f.write('{0},{1},{2}{3}'.format(i[0].replace('"', ""), i[1], i[2], '\n'))
    f.close()
    return None

def duplicates (name1):
    f2 = open(name1, 'r')
    lines = f2.readlines()
    array = np.unique(lines)
    f3 = open("plasmids_filtered_unique.txt", "a")
    for a in array:
        f3.write('{0},{1},{2}{3}'.format(a.split(",")[0], a.split(",")[1], a.split(",")[2], '\n'))
    return None

if __name__ == '__main__':
    name = "Data.csv"
    name1 = "plasmids_filtered.txt"
    write_file(read_file(name))
    duplicates(name1)






