# Drug_Resistance_Analysis
Final degree code

# Introduction
This code can be used for the identification of drug resitance genes of interest within assemblies already aligned with a drug resistance database.
The databases used for this code are: AMRFinder and Megares.
The assemblies used for this code are confidential and should be classified and stored depending on the DB used:
  1. AMRFinder assemblies should be kept in a folder named AMR and extension of the file .fasta
  2. MEGARES assemblies should be kept in a folder named MEGARES and extension of the file .fasta

The outputs have been obtained using: Data_Obtention_Script.sh


# Description of the code workflow
1. Read AB_Genes and stores the information into a TSV file
    1.1 AB_Genes is a document with the names of the antibiotics tested        experimentally.
        Input (TSV):[AntibioticAbreviation AntibioticName  AntibioticClass]
        Output (TSV): [AntibioticClass1 AntibioticClass2 ....]
        
2. Read filtered AMRFinder DB files
    2.1 Filtered by ID>95%, COV>60 and AntibioticClass present in AB_Genes.
    2.2 Formats the result so it can be comparable
        Input (TSV): Column selection,filtration and formatting needed
        Output(CSV):[Isolate_Id,AntibioticClass1,Gene1]
 
 3. Read filtered Megares DB files
    3.1 Filtered by ID>95% (In bash script) and AntibioticClass present in AB_Genes.
    3.2 Formats the result so it can be comparable
        Input (TSV): Column selection,filtration and formatting needed
        Output (CSV): [Isolate_Id,AntibioticClass1,Gene1]
        
4. Comparator
        
        
        
        
