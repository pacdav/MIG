import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
import os

iteration_calcul = 200

iteration_max = 0

wd = os.getcwd()


# Calcul : 

# boucle principale d'usure

for ite_calcul in range(1, iteration_calcul + 1):

#il faut un dossier vide ou il n'y a pas eu de calcul précédent, supprimé le dossier précédent

    if ite_calcul == 1 and os.path.exists("Extracted_files") :
        raise FileExistsError("Déplacez le dossier Extracted_files dans un autre dossier")
# on lit le nombre de cycle par saut de calcul, beta pour la wearbox

    with open(wd + "/iteration_number.txt", "w") as f:
        f.write(str(ite_calcul))
# on lance le calcul et les extractions associées

    os.system("echo yes | abaqus cpus=4 job=projet_MIG interactive")
    os.system("echo yes | abaqus cae nogui=lanceur_fonc_extrac_surf.py")
    os.system("echo yes | abaqus cae nogui=lanceur_fonc_extrac_souscouch.py")

    
# on va récupérer les neouds en contact sur le plan et le cylindre, 
# ici le code a en dur leur position dans l'INP ...    
    lines = []
    with open("projet_MIG.inp", 'r') as f:
        for line in f:
            lines.append(line)
    string_noeuds_contact_cylindre = lines[39730:39743]
    string_noeuds_contact_plan = lines[39746:39759]
    noeuds_contact_cylindre = []
    noeuds_contact_plan = []

    for chaine in string_noeuds_contact_cylindre:
        list_int = list(map(int, chaine[1:-1].split(', ')))
        noeuds_contact_cylindre += list_int
    for chaine in string_noeuds_contact_plan:
        list_int = list(map(int, chaine[1:-1].split(', ')))
        noeuds_contact_plan += list_int

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
    df2_plan = df_plan[(df_plan.y >= -0.4) & (df_plan.x <= 2) & (df_plan.x >= -2)]

    df_cylindre = pd.DataFrame(coord_cylindre_all, columns=['x', 'y'])
    df_cylindre.reset_index(inplace=True)
    df_cylindre['index'] += 1
    df_cylindre.set_index('index', inplace=True)
    
    df2_cylindre = df_cylindre[(df_cylindre.y <= 0.45) & (df_cylindre.x <= 2) & (df_cylindre.x >= -2)]

    step = 1/20
# on récupère le cisaillement du pan et du cylindre
    cis = []
    fichier = open(f"Extracted_files/Extracted_files{ite_calcul}/res_cisaillement_plan.dat", 'r')
    for line in fichier :
        cis.append(list(map(float, line.split())))
    q_plan = np.array(cis)[:, 2:]
    if iteration_max == 0 :
        iteration_max = len(cis[0])
    fichier.close()

    cis = []
    fichier = open(f"Extracted_files/Extracted_files{ite_calcul}/res_cisaillement_cyl.dat", 'r')
    for line in fichier :
        cis.append(list(map(float, line.split())))
    q_cylindre = np.array(cis)[:, 2:]
    fichier.close()

# on récupère le déplacement sur la slave surface, ici le plan

    dep = []
    fichier = open(f"Extracted_files/Extracted_files{ite_calcul}/res_cdisp_plan.dat", 'r')
    for line in fichier :
        dep.append(list(map(float, line.split())))
    dep = np.array(dep)
    d_plan = dep[:, 2:] - dep[:,1:iteration_max-1]
    fichier.close()
# le coef d'usure en dur ...
    alpha = 10*10**(-8)
    beta = 200
    alpha *= beta
# on interpole le déplacement sur le cylindre 
    liste_abs = [df_cylindre.x[noeud] for noeud in noeuds_contact_cylindre]
    list_of_functions = []
    for ite in range(iteration_max-2):
        liste_ord = - d_plan[:, ite]
        list_of_functions.append(interpolate.interp1d(liste_abs, liste_ord))

# on calcul lenergie dissipé sur le cylindre
    h_cylindre = {}
    j = 0
    liste_de_x = []
    for number in noeuds_contact_cylindre:
        x = df_cylindre.x[number]
        sum = 0
        for i in range(iteration_max-2):
            f = list_of_functions[i]
            sum += f(x) * q_cylindre[j, i]
        sum *= alpha
        h_cylindre[x] = sum
        liste_de_x.append(x)
        j += 1
# on déplace la position des noeuds
    for x in liste_de_x:
        h = h_cylindre[x]
        df_aux = df2_cylindre.y[abs(df2_cylindre.x - x) < 0.00001].sort_values(ascending = True)
        corresponding_nodes = df_aux.index
        for i in range(20):
            number = corresponding_nodes[i]
            df_cylindre['y'][number] = df_cylindre['y'][number] - h * (1 - step * i)

    df3_cylindre = df_cylindre[(df_cylindre.y <= 0.4) & (df_cylindre.x <= 2) & (df_cylindre.x >= -2)]
    
# on calcul l'énergie sur le plan
    h_plan = {}
    j = 0
    for number in noeuds_contact_plan:
        x = df_plan.x[number]
        sum = 0
        for i in range(iteration_max-2):
            sum += d_plan[j, i] * q_plan[j, i]
        sum *= alpha
        h_plan[x] = sum
        j += 1
# on applique le déplacement sur les noeud
    for x in df2_plan.x.unique():
        h = h_plan[x]
        df_aux = df2_plan.y[abs(df2_plan.x - x) < 0.00001].sort_values(ascending = False)
        corresponding_nodes = df_aux.index
        for i in range(20):
            number = corresponding_nodes[i]
            df_plan['y'][number] = df_plan['y'][number] - h * (1 - step * i)

    df3_plan = df_plan[(df_plan.y >= -0.4) & (df_plan.x <= 2) & (df_plan.x >= -2)]

    '''
    plt.plot(df2_plan.x, df2_plan.y, '.', label='Plan avant')
    plt.plot(df2_cylindre.x, df2_cylindre.y, '.', label='Cylindre avant')
    plt.plot(df3_cylindre.x, df3_cylindre.y, '.', label='Cylindre après')
    plt.plot(df3_plan.x, df3_plan.y, '.', label='Plan après')
    plt.legend()
    plt.title('Usure au niveau des mailles centrales')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    '''

    # Réécriture d'un nouveau fichier (qui écrasera le précédent). 

    # Réécriture de plan : 
    tab = df_plan.reset_index().to_numpy()
    tab_aux = tab[:, 0].copy()
    nodes = list(tab_aux.astype(int).astype(str))
    list_x = list(tab[:, 1].astype(str))
    list_y = list(tab[:, 2].astype(str))

    nouvelles_coordonnees_plan = []
    for i in range(30684 - 24288):
        nouvelles_coordonnees_plan.append(nodes[i] + ", " + list_x[i] + ", " + list_y[i] + "\n")

    # Réécriture du cylindre : 
    tab = df_cylindre.reset_index().to_numpy()
    tab_aux = tab[:, 0].copy()
    nodes = list(tab_aux.astype(int).astype(str))
    list_x = list(tab[:, 1].astype(str))
    list_y = list(tab[:, 2].astype(str))

    nouvelles_coordonnees_cylindre = []
    for i in range(9355-20):
        nouvelles_coordonnees_cylindre.append(nodes[i] + ", " + list_x[i] + ", " + list_y[i] + "\n")

    # Réécriture du fichier : 
    with open("projet_MIG.inp", "w") as f:
        for i in range(20):
            f.write(lines[i])
        for i in range(9355 - 20):
            f.write(nouvelles_coordonnees_cylindre[i])
        for i in range(9355, 24288):
            f.write(lines[i])
        for i in range(30684 - 24288):
            f.write(nouvelles_coordonnees_plan[i])
        for i in range(30684, len(lines)):
            f.write(lines[i])
    
    # On copie dans les nouveaux odb et inp parmi les fichiers extraits.


    os.chdir(wd)
    os.system(f"cp projet_MIG.inp projet_MIG_{ite_calcul}.inp")
    # os.system(f"copy projet_MIG.odb projet_MIG_{ite_calcul}.odb")
    os.system(f"mv projet_MIG_{ite_calcul}.inp " + wd + f"/Extracted_files/Extracted_files{ite_calcul}")
    # os.system(f"move projet_MIG_{ite_calcul}.odb /Users/flaneerUser/Desktop/MIG_simulation/Extracted_files/Extracted_files{ite_calcul}")

    if ite_calcul == iteration_calcul:
        os.system(f"cp projet_MIG.odb projet_MIG_{ite_calcul}.odb")
        os.system(f"mv projet_MIG_{ite_calcul}.odb " + wd + f"/Extracted_files/Extracted_files{ite_calcul}")


    # os.system(f"find dev /maxdepth -1 -name 'projet_MIG*' | xargs cp -t /Users/flaneerUser/Desktop/MIG_simulation/Extracted_files/Extracted_files{ite_calcul}")

    os.system("python courbes.py")
    print(f"ITERATION {ite_calcul} : DONE. ")


# os.system(f"find dev /maxdepth -1 -name 'projet_MIG*' -type f -delete")
# os.chdir("Fichiers avant modification")
# os.system("find dev /maxdepth -1 -name 'projet_MIG*' | xargs cp -t /Users/flaneerUser/Desktop/MIG_simulation")




print("CALCUL TERMINE.")

