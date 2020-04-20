# ------------------------------------LIBRARY START-------------------------------------------#
import glob
from difflib import SequenceMatcher
import plotly.graph_objects as go
import numpy as np
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
    auxiliar_array = []
    for line in lines:
        # Splitting TSV format, and returning key position (CLASS)
        if line.split('\t')[pos] in auxiliar_array:
            pass
        else:
            auxiliar_array.append(line.split('\t')[pos])
    f.close()
    return auxiliar_array

# Formats the result in a manner that can be compared
# Replace \t per None __ - per None and MAJ
# Re-writes AntibioticClass
def ab_genes_writer(auxiliary_array):
    # String formatting CAPS and blanks/-
    for var1 in auxiliary_array:
        var1 = var1.replace(' ', '').replace('-', '').upper()
        f = open("AB_Genes_mod.txt", "a")
        f.write('{0}{1}'.format(var1, '\n'))
        f.close()

# Open file = Isolate
# Select the variables of intrest [Identifier, Gene, Class, Coverage, Id]
# Read lines
# Compares with the AB_Genes file (comparator function)
# And filters by ID && Coverage
#
def amr_reader(file):
    key_vars = [1, 5, 10, 15, 16]
    static_auxiliary_array = [file[4:10]]
    f = open(file, 'r')
    lines = f.readlines()
    del lines[0]
    for line in lines:
        # if the comparison is not needed, comment code line 57.
        var_1 = True
        var_1 = comparator(line.split('\t')[10])
        if var_1 and float(line.split('\t')[15]) > 60 and float(line.split('\t')[16]) > 95:
            static_auxiliary_array.append([line.split('\t')[var] for var in key_vars])
            var_1 = False
        else:
            pass
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
#
def megares_reader(file):
    static_auxiliary_array = [file[8:15]]
    string_values = [4, 2]
    f = open(file, 'r')
    lines = f.readlines()
    formatting_array_genes = []
    for line in lines:
        auxiliary_array = line.split('\t')
        auxiliary_array_2 = auxiliary_array[0].split('\n')[0].split('|')
        auxiliary_array_3 = []
        # if the comparison is not needed, comment code line 100.
        var_1 = True
        var_1 = comparator(auxiliary_array_2[2].upper())
        if auxiliary_array_2[1] == 'Drugs' and float(auxiliary_array[2])/float(auxiliary_array[3]) > 0.6 \
            and float(auxiliary_array[4]) > 95 and var_1:
            # Python automatically erase \n when conversion to float entity
            auxiliary_array_3.append(auxiliary_array[1])
            for val in string_values:
                auxiliary_array_3.append(auxiliary_array_2[val].upper())
            auxiliary_array_3.append(str(float(auxiliary_array[2])/float(auxiliary_array[3])*100))
            auxiliary_array_3.append(auxiliary_array[4].split('\n')[0])
            var_1 = False
            formatting_array_genes.append(auxiliary_array_3)
    # Check in order to avoid appending void arrays.
    if len(formatting_array_genes) > 0:
        static_auxiliary_array.append(formatting_array_genes)
    else:
        pass
    f.close()
    return static_auxiliary_array

# Re-writes the format in case of an Isolate having more than 1 AntibioticClass
# Replace the list format and the \ to create CSV.
def megares_writer(auxiliary_array):
    f = open("megares_data.txt", "a")
    if len(auxiliary_array)>1:
        for i in auxiliary_array[1]:
            auxiliary_string = str(i).replace('\'', '' )
            f.write('{0},{1}{2}'.format(auxiliary_array[0], auxiliary_string[1:-2], '\n'))
    f.close()
    return None

# Uses the SequenceMatcher.ratio function to compare the similarities of the strings
# If is higher that the threshold determined as 0.9 is close enough
def comparator(input_data):
    f = open('AB_Genes_mod.txt', 'r')
    lines = f.readlines()
    comparation_array = []
    for line in lines:
        comparation_array.append(SequenceMatcher(None, input_data, line.split('\n')[0]).ratio())
        #print(max(comparation_array))
    if max(comparation_array) > 0.9:
        f.close()
        return True
    else:
        f.close()
        return False


def reformatting_file(file):
    file = open(file, "r+")
    file.truncate(0)
    file.close()
    return None


def headers():
    f = open('AB_Genes_mod.txt', 'r+')
    lines = f.readlines()
    auxiliary_array = []
    for line in lines:
            auxiliary_array.append(line.split('\n')[0])
    f.close()
    return auxiliary_array


def isolate_values(file):
    f_1 = open(file, 'r')
    lines = f_1.readlines()
    file_dictionary = {}
    for line in lines:

        # Resetting global loop variables
        auxiliary_array = []
        auxiliary_dictionary = {}

        # This parsing can be done, due it si a custom file. the original formatting is already known.
        # [ISOLATE_1, [information from 1st gene],[information from 2nd gene],...[information from nth gene]]
        # Each gen cell, contains a total of 5 data keys. The total cell positions, is multiple of 5.
        string_modifier = line.split('\n')[0]
        string_modifier = string_modifier.replace('\'', '').replace('[', '')\
            .replace(']', '').replace(' ', '').split(',')

        # Data parser -- Only Isolates
        isolate_key_values = string_modifier[0]

        # Data parser -- Gene
        values_gene = string_modifier[1:]

        # Rearranging items into sub-lists
        for position in range(0, len(values_gene)-5, 5):
            auxiliary_array.append(values_gene[position:position+5])

        # Arranging dictionary keys from AB_Genes_mod file for simplicity
        f_2 = open('AB_Genes_mod.txt', 'r+')
        anti_lines = f_2.readlines()
        for anti in anti_lines:
            auxiliary_dictionary[anti.split('\n')[0]] = []
        f_2.close()

        # Arranging gene data as values for the dictionary if the criteria is met
        # Relevant data:
        # Find key position of auxiliary_array[gene][2] in auxiliary_dictionary
        # Append associated gene name identifier
        for gene in range(len(auxiliary_array)):
            comparator_array = []
            for position in range(len(list(auxiliary_dictionary.keys()))):
                comparator_array.append(SequenceMatcher(None,
                                                        auxiliary_array[gene][2],
                                                        list(auxiliary_dictionary.keys())[position]).ratio())
            match = comparator_array.index(max(comparator_array))
            if comparator_array[match] > 0.9:
                key = list(auxiliary_dictionary.keys())[match]
                if key in auxiliary_dictionary:
                    auxiliary_dictionary[key].append(auxiliary_array[gene][1])
                else:
                    auxiliary_dictionary[key] = auxiliary_array[gene][1]

        # Nesting dictionaries in order to obtain the isolate file as a key, and its dictionary as a value.
        # file_dictionary = {ISOL_1 : auxiliary_dictionary, ISOL_2 : auxiliary_dictionary...}
        file_dictionary[isolate_key_values] = auxiliary_dictionary

    f_1.close()
    return file_dictionary


# -----------------------------------FUNCTIONS END--------------------------------------------#
# -----------------------------------SCRIPT STARTS--------------------------------------------#


if __name__ == '__main__':
    # Clean initial files, in case to be modified, all files are empty for formatting
    reformatting_file('AB_Genes_mod.txt')
    reformatting_file('amr_data.txt')
    reformatting_file('megares_data.txt')
    #
    ab_genes_writer(ab_genes_reader('AB_Genes.txt', pos=2))
    antibiotic_unique = ab_genes_reader('AB_Genes.txt', pos=2)
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
        #print(len(values_for_table[i]))
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

