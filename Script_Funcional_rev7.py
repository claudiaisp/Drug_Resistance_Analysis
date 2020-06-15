# ------------------------------------LIBRARY START-------------------------------------------#
import glob
from difflib import SequenceMatcher
import numpy as np
#import plotly.graph_objects as go
# ------------------------------------LIBRARY END---------------------------------------------#
# ----------------------------------FUNCTIONS START-------------------------------------------#

# Open the AB_Genes.txt document,
# Reads the lines,
# Split the lines, using \t
# Selects the AntibioticClass column, determined as variable (pos)


def ab_genes_reader(file, pos):
    f = open(file, 'r')
    lines = f.readlines()
    del lines[0]  # Headers
    auxiliary_array = []
    for line in lines:
        # Splitting TSV format, and returning key position
        if line.split('\t')[pos].upper() in auxiliary_array:
            pass
        else:
            auxiliary_array.append(line.split('\t')[pos].upper())
    f.close()
    return auxiliary_array

# Formats the result in a manner that can be compared
# Replace \t per None __ - per None and MAJ
# Re-writes AntibioticClass


def ab_genes_writer(auxiliary_array):
    # String formatting CAPS and blanks/-
    f = open("AB_Genes_mod.txt", "a")
    for var1 in auxiliary_array:
        var1 = var1.replace(' ', '').replace('-', '').upper()
        f.write('{0}{1}'.format(var1, '\n'))
    f.close()


def genes_writer(auxiliary_array):
    auxiliary_array_2 = []
    for val in auxiliary_array:
        if val in auxiliary_array_2:
            pass
        else:
            auxiliary_array_2.append(val.split('\n')[0])
    f = open("AB_Genes_genes.txt", "a")
    for val in auxiliary_array_2:
        f.write('{0}{1}'.format(val, '\n'))
    f.close()
    return auxiliary_array_2

# Open file = Isolate
# Select the variables of intrest [Identifier, Gene, Class, Coverage, Id]
# Read lines
# Compares with the AB_Genes file (comparator function)
# And filters by ID && Coverage
#


def amr_reader(file):
    key_vars = [5, 13, 16, 10]
    static_auxiliary_array = [file[4:11]]
    f = open(file, 'r')
    lines = f.readlines()
    del lines[0]
    for line in lines:
        auxiliary_array = []
        var_3 = False
        var_1 = comparator(line.split('\t')[10])
        var_2 = gene_validation(line.split('\t')[5])
        if var_1 and var_2 and float(line.split('\t')[15]) > 60 and float(line.split('\t')[16]) > 95:
            var_1 = False
            var_2 = False
            auxiliary_array.append(line.split('\t')[5])
            auxiliary_array.append(float(line.split('\t')[3]) - float(line.split('\t')[2]))
            auxiliary_array.append(line.split('\t')[16])
            auxiliary_array.append(line.split('\t')[10])
            var_3 = True
        else:
            pass
        if var_3:
            static_auxiliary_array.append(auxiliary_array)
            var_3 = False
    f.close()
    return static_auxiliary_array

# Re-writes the format in case of an Isolate having more than 1 AntibioticClass
# Replace the list format and the \ to create CSV.


def amr_writer(auxiliary_array):
    f = open("amr_data.txt", "a")
    if len(auxiliary_array) > 1:
        for i in range(1, len(auxiliary_array)):
            auxiliary_string = str(auxiliary_array[i]).replace('\'',"").replace('[',"").replace(']',"")
            f.write('{0},{1}{2}'.format(auxiliary_array[0], auxiliary_string, '\n'))
    f.close()
    return None

# Open file = Isolate
# Select the variables of interest [String, Coverage_1, Coverage_2, Id]
# Read lines
# Formats MAJ
# Compares with the AB_Genes file (comparator function)
# And filters by ID && Coverage

def megares_reader(file):
    static_auxiliary_array = [file[8:15]]
    f = open(file, 'r')
    lines = f.readlines()
    formatting_array_genes = []
    for line in lines:
        auxiliary_array = line.split('\t')
        node = auxiliary_array[0]
        auxiliary_array_2 = auxiliary_array[1].split('|')
        auxiliary_array_3 = []
        # if the comparison is not needed, comment code line 100.
        # The following comparison has to be made with the antibiotic Class and the Gene.
        var_1 = True
        var_1 = comparator(auxiliary_array_2[2].upper())
        if auxiliary_array_2[-1] == 'RequiresSNPConfirmation':
            var_2 = gene_validation(auxiliary_array_2[-2].upper())
        else:
            var_2 = gene_validation(auxiliary_array_2[-1].upper())
        if auxiliary_array_2[1] == 'Drugs' and float(auxiliary_array[3])/float(auxiliary_array[2]) > 0.6 \
            and float(auxiliary_array[4]) > 95 and var_1 and var_2:
            auxiliary_array_3.append(node)
            auxiliary_array_3.append(auxiliary_array_2[4].upper())
            auxiliary_array_3.append(float(auxiliary_array[8])-float(auxiliary_array[7]))
            auxiliary_array_3.append(float(auxiliary_array[4].split('\n')[0]))
            auxiliary_array_3.append(auxiliary_array_2[2].upper())
            # Python automatically erase \n when conversion to float entity
            var_1 = False
            var_2 = False
            formatting_array_genes.append(auxiliary_array_3)
    f.close()
    # Check in order to avoid appending void arrays.
    if len(formatting_array_genes) > 0:
        static_auxiliary_array.append(formatting_array_genes)
    else:
        pass
    return static_auxiliary_array

# Re-writes the format in case of an Isolate having more than 1 AntibioticClass
# Replace the list format and the \ to create CSV.


def megares_writer(auxiliary_array):
    f = open("megares_data.txt", "a")
    if len(auxiliary_array) > 1:
        for i in auxiliary_array[1]:
            auxiliary_string = str(i).replace('\'', '')
            f.write('{0}, {1}{2}'.format(auxiliary_array[0], auxiliary_string[1:-2], '\n'))
    f.close()
    return None

# Uses the SequenceMatcher.ratio function to compare the similarities of the strings
# If is higher that the threshold determined as 0.9 is close enough


def comparator(input_data):
    comparation_array = []
    for line in _ANTIBIOTIC:
        comparation_array.append(SequenceMatcher(None, input_data, line.split('\n')[0]).ratio())
    if max(comparation_array) > 0.9:
        return True
    else:
        return False


def gene_validation(input_data):
    comparation_array = []
    for line in _GENE:
        # modify input data to len size of _GENE
        if line.find('*') > 0:
            var_1 = input_data.upper()[:line.index('*')]
            var_1_true = line.split('\n')[0][:-1]
            comparation_array.append(SequenceMatcher(None, var_1, var_1_true).ratio())
        else:
            var_1 = input_data.upper()[:len(line)]
            var_1_true = line.split('\n')[0]
            comparation_array.append(SequenceMatcher(None, var_1, var_1_true).ratio())
    if max(comparation_array) > 0.9:
        return True
    else:
        return False


def read_plasmids (file):
    f = open(file, 'r')
    lines = f.readlines()
    auxiliary_array = []
    array = []
    for line in lines:
        array = [line.split('\t')[0], line.split('\t')[2].split(' ')[0], str(" ".join(line.split('\t')[2].split(" ")[1:])), line.split('\t')[3], line.split('\t')[0].split("_")[3]]
        if float(line.split('\t')[0].split("_")[5]) > 80 and float(line.split('\t')[0].split("_")[3]) > 500 and float(line.split('\t')[3]) > 95:
            auxiliary_array.append(array)

        else:
            pass
    f.close()
    return auxiliary_array


def write_plasmids(auxiliary_array, name):
    f = open(name, "a")
    for i in auxiliary_array:
        print('{0},{1},{2}'.format(file[23:30], str(i).replace("[","").replace("]",""), '\n'))
        f.write('{0},{1},{2}'.format(file[23:30], str(i).replace("[","").replace("]",""), '\n'))
    f.close()
    return None


def duplicates(name1):
    f2 = open(name1, 'r')
    lines = f2.readlines()
    array = np.unique(lines)
    f3 = open(name1, "a")
    for i in array:
        f3.write('{0}'.format(i ,'\n'))
    return None


def adjust_results(file):
    f = open(file, 'r')
    lines = f.readlines()
    del lines[0]
    static_array = []
    for line in lines:
        try:
            pos = line.split('\t')[4].split(' ').index('plasmid')
            array = [file, line.split('\t')[1], line.split('\t')[7], line.split('\t')[2], line.split('\t')[3],
                     line.split('\t')[4].split(' ')[pos:]]
            static_array.append(array)
        except:
            pass

    fout = open('plasmid_name.txt', "a")
    for i in static_array:
        fout.write('{0},{1}'.format(str(i).replace('[', '').replace(']', '').replace("'", '').replace('"', ''), '\n'))
    fout.close()


# -----------------------------------FUNCTIONS END--------------------------------------------#
# -----------------------------------SCRIPT STARTS--------------------------------------------#


if __name__ == '__main__':
    # Clean initial files, in case to be modified, all files are empty for formatting
    # CAP vars indicate const vars.

    _GENE = genes_writer(ab_genes_reader('AB_Genes.txt', pos=3))
    _ANTIBIOTIC = ab_genes_reader('AB_Genes.txt', pos=2)

    names = glob.glob('AMR/*.fasta')
    names.sort()
    for file in names:
        amr_writer(amr_reader(file))

    names = glob.glob('MEGARES/*.fasta')
    names.sort()
    for file in names:
        megares_writer(megares_reader(file))

    names = glob.glob('Brooks_et_al/*.fasta')
    names.sort()
    for file in names:
        write_plasmids(read_plasmids(file),"Brooks_filtered.txt")
    duplicates('Brooks_filtered.txt')

    names = glob.glob('PLSDB/*.fasta')
    names.sort()
    for file in names:
        write_plasmids(read_plasmids(file), "'PLSDB_filtered.txt'")
    duplicates('PLSDB_filtered.txt')

    names = glob.glob('Filtered_results_*')
    for file in names:
        adjust_results(file)

'''
    # AMR/MEGARES file writing, uniques and reformatted values.
    # [ISOLATE_1, [information from 1st gene]
    # [ISOLATE_1, [information from 2st gene]
    # [ISOLATE_2, [information from 1st gene]
    # ...
    # [ISOLATE_n, [information from nst gene]
    # CASE OF NOT METTING CRETERIA:
    # Whole isolate do not met criteria: [ISOLATE_x]
    # Some gene from an ISOLATE_x does NOT met the criteria: ISOLATE_1, [],[Information from gene that met criteria],..,[]]
    names = glob.glob('AMR/*.fasta')
    names.sort()
    for file in names:
        amr_writer(amr_reader(file))
    names = glob.glob('MEGARES/*.fasta')
    names.sort()
    for file in names:
        megares_writer(megares_reader(file))

# -------------------------------------SCRIPT END--------------------------------------------#
# -------------------------------------TEST START--------------------------------------------#
    # Headers creations for table
    headers = headers()
    headers.insert(0, 'ISOLATES')

    # ISOLATE value creation for table
    values_isolate_amr = isolate_values('amr_data.txt')
    values_isolate_megares = isolate_values('megares_data.txt')

    # Merging both nested dictionaries,the order is kept. amr_dict - megares_dict
    values_isolate_amr.update(values_isolate_megares)
    #print(values_isolate_amr)


    # Grouping data for Results
    isolates = list(values_isolate_amr.keys())

    # Generating data points for the table, keeping the format shown below.
    # Accessing into the nested dictionary,
    # For key in values_isolate_amr:  <-- Access to Isolate files
    #   for info in values_isolate_amr[key]: <-- Access to key values from ISOLATED files
    # To obtain access to certain values of certain [antibiotic] of the isolate file number [N]:
    #   values_isolate_amr[N].get(antibiotic)
    #
    #print(headers)
    values_for_table = []
    for values in range(len(antibiotic_unique)):
        auxiliary_array = []
        for key in values_isolate_amr:
            auxiliary_array.append(values_isolate_amr[key][antibiotic_unique[values].upper()])
        values_for_table.append(auxiliary_array)
    values_for_table.insert(0, [isol for isol in isolates])
    #print(values_for_table)
    for i in range(len(values_for_table)):
        print(len(values_for_table[i]))
    #print('-')
    #print(len(values_for_table))

    # Table format, value position,
    #   ISOLATES ID     Header_1                Header_2                ...         Header_N
    #   ISOLATE_1       antibiotic1_isolate1    antibiotic2_isolate1  ...         antibioticN_isolate1
    #   ISOLATE_2       antibiotic1_isolate2    antibiotic2_isolate2  ...         antibioticN_isolate2
    #   ISOLATE_3       antibiotic1_isolate3    antibiotic2_isolate3  ...         antibioticN_isolate3
    #      ...              ...                     ...                 ...            ....
    #   ISOLATE_N       antibiotic1_isolateN    antibiotic2_isolateN              antibioticN_isolateN
    # ----------->
    # [[ISOLATES], [ANTIBIOTIC_1 FOR ALL ISOLATES], [ANTIBIOTIC_2 FOR ALL ISOLATES],  (continuation below)
    #   [ANTIBIOTIC_3 FOR ALL ISOLATES], ... , [ANTIBIOTIC_N FOR ALL ISOLATES]]
# -------------------------------------TEST END--------------------------------------------#
'''