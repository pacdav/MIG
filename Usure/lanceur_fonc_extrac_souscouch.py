# -*- coding: Utf-8 -*-
import os
import sys
WorkingDirectory = os.getcwd()
# sys.path.insert(0, WorkingDirectory + '//fretting')  

with open("iteration_number.txt", "r") as f:
    line = f.readline()
ite_calcul = int(line)
string2 = "/Extracted_files/Extracted_files" + str(ite_calcul)
sys.path.insert(0, WorkingDirectory + string2) 

from extrraction_sous_couche import extraction_var_sous
extraction_var_sous(WorkingDirectory + '/projet_MIG.odb', 4.0,0.40,0.020,WorkingDirectory +string2)
