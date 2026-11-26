# -*- coding: Utf-8 -*-
import os
import sys


WorkingDirectory = os.getcwd()

with open("iteration_number.txt", "r") as f:
    line = f.readline()
ite_calcul = int(line)

if ite_calcul == 1:
    os.makedirs("Extracted_files")
os.chdir("Extracted_files")
string = "Extracted_files" + str(ite_calcul)
os.makedirs(string)

string2 = "/Extracted_files/Extracted_files" + str(ite_calcul)
sys.path.insert(0, WorkingDirectory + string2) 

from extraction_surface import extraction_var_surf
extraction_var_surf(WorkingDirectory + '/projet_MIG.odb', 4.0,4.0,0.40,0.020, WorkingDirectory + string2)
