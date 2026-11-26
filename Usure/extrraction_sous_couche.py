
from abaqus import *
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
from math import *
from numpy import *
import time
import visualization
import sys
from decimal import *
from os import chdir


def extraction_var_sous(nom_fichier_odb,largeur_plan,profondeur,tailleMaille,DossierDatFileInc):

    # Fonction d extraction des donnees qui nous interessent dans lODB abaqus
    odb = openOdb(path=str(nom_fichier_odb))		# on recupere l odb
    myAssembly = odb.rootAssembly					# on recupere l assembly
	
    print "\nDEBUT FCT EXCTRACTION\n"
	
    nbr_frames = len(odb.steps['Fretting_fatigue'].frames)		# donne le nombre de frame (commence a zero)
    nbr_elem_par_box_plan     = int(round((largeur_plan/tailleMaille)*(profondeur/tailleMaille)))
    nbr_noeuds_setContactPlan = int(round(largeur_plan /tailleMaille+1))
    
    print "Parametres d entree :\n"
    
    print "largeur_plan = %f\n" % largeur_plan
    print "profondeur   = %f\n" % profondeur
    print "tailleMaille = %f\n" % tailleMaille
    print "nbr_elem_par_box_plan     = %f\n" % nbr_elem_par_box_plan
    print "nbr_noeuds_setContactPlan = %f\n" % nbr_noeuds_setContactPlan
    

    dernierIncrAbsoluStepPrecedente = 0

    nbStep       = len(odb.steps.keys())
    nomStep      = odb.steps.keys()
    stepFret     = odb.steps['Fretting_fatigue']			# on recupere la step Fretting_fatigue
    
    StepKeys     = odb.steps.keys()
    stepInitiale = odb.steps[StepKeys[0]]
    stepNum      = zeros((nbStep,1),'i')
    nbIncrStep   = zeros((nbStep,1),'i')
    
    dernierIncrAbsoluStep = zeros((nbStep,1),'i')
    MatIncrStep = zeros((nbStep,3),'i')
    i=0 
    for step in nomStep:
        stepNum[i]=int(odb.steps[step].number)
        # On evalue le nombre de frames par step via le nombre de frames enregistrees
        nbIncrStep[i]                   = len(odb.steps[step].frames)
        dernierIncrAbsoluStep[i]        = dernierIncrAbsoluStepPrecedente + nbIncrStep[i]
        dernierIncrAbsoluStepPrecedente = dernierIncrAbsoluStep[i]
        i+=1


    ### Enregistrement de la matrice de suivi des increments de calcul ----------------------------------
    CompleteName = os.path.join(DossierDatFileInc, 'res_incr_step.dat')
    res          = open(CompleteName, 'w')
    nb_frames=0
    for j in range(nbStep):
        res.write('%i\t%i\t%i\n' % (stepNum[j][0],nbIncrStep[j][0],dernierIncrAbsoluStep[j][0]))
    res.close()
	
	
### Enregistrement de la matrice des points temporels de la step fretting  #################
    CompleteName = os.path.join(DossierDatFileInc, 'res_step_fretting.dat')
    res          = open(CompleteName,'w')
    for fram in stepFret.frames:
        nb_frames = nb_frames + 1
        Tps       = fram.frameValue
        res.write('%0.6f\n' % Tps)
    res.close()


    # basePlan=myAssembly.nodeSets['BASE_PLAN']
    boxPlan = myAssembly.elementSets['BOX_PLAN']

###################################################################################   
##   ____  _                                   _____           _ _               ##
##  / ___|(_) __ _ _ __ ___   __ _            | ____|_ __  ___(_) | ___  _ __    ##
##  \___ \| |/ _` | '_ ` _ \ / _` |   _____   |  _| | '_ \/ __| | |/ _ \| '_ \   ##
##   ___) | | (_| | | | | | | (_| |  |_____|  | |___| |_) \__ \ | | (_) | | | |  ##
##  |____/|_|\__, |_| |_| |_|\__,_|           |_____| .__/|___/_|_|\___/|_| |_|  ##
##           |___/                                  |_|                          ##
###################################################################################  
##     
##    EXTRACTION DONNEES DE CONTRAINES SIGMA POUR LE NSET BoxPan    ################################
    print "\nExtraction contraintes"
    # On va s'occuper ici d'enregistrer les donnees de contraintes aux noeuds de notre boite.
    # Les noeuds sont ranges dans un certain ordre : parcours colonne par colonne avec abscisses 
    # croissants. Cad on part de la colonne la plus a gauche de la boite et on recupere les donnees
    # en descendant dans cette colonne, puis dans celle qui se trouve a droite.. etc.

    # On stocke les valeurs des contraintes S11,S22,S33 et S12 dans des variables Sii.
    # On sinteresse seulement aux instants extremes du fretting afin de minimiser
    # le temps de calcul.  --> max correspond a t=0.25 et min a t=0.75
    # On stocke sous le format integer et en kPa afin de limiter le volume du fichier.
    print "nbr_elem_par_box_plan= ",nbr_elem_par_box_plan,"\n"

	# Partie en option : enregistrement des champs de contrainte pour tous les increments
    
    S11=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    S22=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    S33=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    S12=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    print "nbr_elem_par_box_plan= ",nbr_elem_par_box_plan,"\n"
    ind_colonne=0
    for incr in stepFret.frames :
        fieldValues_S = incr.fieldOutputs['S'].getSubset(region=boxPlan).values
        for dat in range(int(nbr_elem_par_box_plan)) : 
            S11[dat][ind_colonne]=1000*fieldValues_S[dat].data[0]
            S22[dat][ind_colonne]=1000*fieldValues_S[dat].data[1]
            S33[dat][ind_colonne]=1000*fieldValues_S[dat].data[2]
            S12[dat][ind_colonne]=1000*fieldValues_S[dat].data[3]
        ind_colonne+=1
    nameS11 = os.path.join(DossierDatFileInc, 'S11.dat')
    nameS22 = os.path.join(DossierDatFileInc, 'S22.dat')
    nameS33 = os.path.join(DossierDatFileInc, 'S33.dat')
    nameS12 = os.path.join(DossierDatFileInc, 'S12.dat')
    s11  =  open(nameS11, 'w')
    s22  =  open(nameS22, 'w')
    s33  =  open(nameS33, 'w')
    s12  =  open(nameS12, 'w')
    
    for jj in range(int(nb_frames)) :
        for ii in range(int(nbr_elem_par_box_plan)) : 
            s11.write('%0.2f\t' % (S11[ii][jj]))
            s22.write('%0.2f\t' % (S22[ii][jj]))
            s33.write('%0.2f\t' % (S33[ii][jj]))
            s12.write('%0.2f\t' % (S12[ii][jj]))
        s11.write('\n')
        s22.write('\n')
        s33.write('\n')
        s12.write('\n')
    s11.close()
    s22.close()
    s33.close()
    s12.close()
    ###################################################################################
    
    
# PARTIE EXTRACTION DES COMPOSANTES DE DEFORMATION ELASTIQUE E11 E22 ET E12
    print "\nExtraction deformation"
    Eps11=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    Eps22=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    Eps33=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    Eps12=zeros((nbr_elem_par_box_plan,nb_frames))*1.
    print "nbr_elem_par_box_plan= ",nbr_elem_par_box_plan,"\n"
    ind_colonne=0
    for incr in stepFret.frames :
        fieldValues_Eps = incr.fieldOutputs['E'].getSubset(region=boxPlan).values
        for dat in range(int(nbr_elem_par_box_plan)) : 
            Eps11[dat][ind_colonne] = 1000000 * fieldValues_Eps[dat].data[0]
            Eps22[dat][ind_colonne] = 1000000 * fieldValues_Eps[dat].data[1]
            Eps33[dat][ind_colonne] = 1000000 * fieldValues_Eps[dat].data[2]
            Eps12[dat][ind_colonne] = 1000000 * fieldValues_Eps[dat].data[3]
        ind_colonne+=1
    nameEps11 = os.path.join(DossierDatFileInc, 'Eps11.dat')
    nameEps22 = os.path.join(DossierDatFileInc, 'Eps22.dat')
    nameEps33 = os.path.join(DossierDatFileInc, 'Eps33.dat')
    nameEps12 = os.path.join(DossierDatFileInc, 'Eps12.dat')
    
    
    eps11=open(nameEps11, 'w')
    eps22=open(nameEps22, 'w')
    eps33=open(nameEps33, 'w')
    eps12=open(nameEps12, 'w')
    
    for jj in range(int(nb_frames)) :
        for ii in range(int(nbr_elem_par_box_plan)) : 
            eps11.write('%0.2f\t' % (Eps11[ii][jj]))
            eps22.write('%0.2f\t' % (Eps22[ii][jj]))
            eps33.write('%0.2f\t' % (Eps33[ii][jj]))
            eps12.write('%0.2f\t' % (Eps12[ii][jj]))
        eps11.write('\n')
        eps22.write('\n')
        eps33.write('\n')
        eps12.write('\n')
    eps11.close()
    eps22.close()
    eps33.close()
    eps12.close()
