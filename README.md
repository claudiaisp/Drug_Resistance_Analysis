# Drug_Resistance_Analysis
Final degree code

# Introduction
This code aims to be used for the identification of drug resitance genes of interest within assemblies already aligned with a drug resistance database.
The databases this code treats are: AMRFinder, MEGARes, PLSDB and Brooks et al.
The assemblies used for this code are confidential and should be classified and stored depending on the DB used:
  1. AMRFinder assemblies should be kept in a folder named AMR and extension of the file .fasta
  2. MEGARES assemblies should be kept in a folder named MEGARES and extension of the file .fasta
  3. PLSDB assemblies should be kept in a folder named PLSDB and extension of the file .fasta
  4. Brooks et al. assemblies should be kept in a folder named Brooks_et_al and extension of the file .fasta

The outputs have been obtained using: Data_Obtention_Script.sh

# Requirements
import glob
from difflib import SequenceMatcher
import numpy as np
        
        
        
        
