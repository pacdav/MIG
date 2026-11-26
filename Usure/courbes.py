import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import shutil


wd = os.getcwd()
nb_iteration = len(os.listdir("Extracted_files"))

if os.path.exists(wd + "/Evolution_maillage"):
    shutil.rmtree(wd + "/Evolution_maillage")
if os.path.exists(wd + "/Profils de pression"):
    shutil.rmtree(wd + "/Profils de pression")

os.makedirs(wd + "/Evolution_maillage")
os.makedirs(wd + "/Profils de pression")

for ite in range(nb_iteration + 1):
    os.chdir(wd)
    lines = []
    if ite==0:
        with open("Fichiers avant modification/projet_MIG.inp", 'r') as f:
            for line in f:
                lines.append(line)
    
    else:
        with open(f"Extracted_files/Extracted_files{ite}/projet_MIG_{ite}.inp", 'r') as f:
            for line in f:
                lines.append(line) 

    os.chdir("Evolution_maillage")

    liste_aux = []
    for i in range(24288, 30684):
        liste_aux.append(list(map(float, lines[i][:-1].split(','))))
    coordonnees_plan = np.array(liste_aux)

    liste_aux = []
    for i in range(20, 9355):
        liste_aux.append(list(map(float, lines[i][:-1].split(','))))
    coordonnees_cylindre = np.array(liste_aux)    
    keys_plan = coordonnees_plan[:, 0]
    coord_plan_all = coordonnees_plan[:, 1:3]
    keys_cylindre = coordonnees_cylindre[:, 0]
    coord_cylindre_all = coordonnees_cylindre[:, 1:3]
    df_plan = pd.DataFrame(coord_plan_all, columns=['x', 'y'])
    df_plan.reset_index(inplace=True)
    df_plan['index'] += 1
    df_plan.set_index('index', inplace=True)

    # Ne pas modifier les valeurs de la ligne suivante !!

    df2_plan = df_plan[(df_plan.y >= -0.05) & (df_plan.x <= 0.75) & (df_plan.x >= -0.75)]
    df_cylindre = pd.DataFrame(coord_cylindre_all, columns=['x', 'y'])
    df_cylindre.reset_index(inplace=True)
    df_cylindre['index'] += 1
    df_cylindre.set_index('index', inplace=True)
    df2_cylindre = df_cylindre[(df_cylindre.y <= 0.05) & (df_cylindre.x <= 0.75) & (df_cylindre.x >= -0.75)]
    plt.plot(df2_plan.x, df2_plan.y, '.', label='Plan')
    plt.plot(df2_cylindre.x, df2_cylindre.y, '.', label='Cylindre')
    plt.legend()
    plt.title('Usure au niveau des mailles centrales')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(f"Usure_iteration_{ite}")
    plt.clf()

ite = np.linspace(1, nb_iteration, nb_iteration)
ite = ite.astype('uint8')
pression = []
step = 0.02
nb_points = int(4 / step) + 1
X = np.linspace(-2, 2, nb_points)
for i in ite :
    os.chdir(wd)
    with open(f'Extracted_files/Extracted_files{i}/res_pression_plan.dat', 'r') as fichier :  
        all_pressures = []
        for line in fichier :
            all_pressures.append(float(line))
        pression.append(max(all_pressures))

    plt.clf()
    os.chdir("Profils de pression")
    plt.plot(X, all_pressures, label=f"Itération {i}")
    plt.xlabel('x')
    plt.ylabel('Pression')
    plt.title('Pression')
    plt.legend()
    plt.savefig(f'profil_pression_{i}')

plt.clf()
os.chdir(wd)

with open('Extracted_files/Extracted_files1/res_pression_plan.dat', 'r') as fichier :  
    all_pressures = []
    for line in fichier :
        all_pressures.append(float(line))
    plt.plot(X, all_pressures, label=f"Itération 1")

os.chdir(wd)
for i in range(1, nb_iteration // 10 + 1):
    with open(f'Extracted_files/Extracted_files{i * 10}/res_pression_plan.dat', 'r') as fichier :  
        all_pressures = []
        for line in fichier :
            all_pressures.append(float(line))
    plt.plot(X, all_pressures, label=f"Itération {i * 10}")

os.chdir("Profils de pression")
plt.legend()
plt.xlabel('x')
plt.ylabel('Pression')
plt.title('Evolution du profil de pression')
plt.savefig("pression_evolution")

plt.clf()
plt.plot(ite, pression)
plt.title("Pression maximale en fonction de l'itération")
plt.xlabel('Itération')
plt.ylabel('Pression maximale')
plt.savefig("pression_max")
