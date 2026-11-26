

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

def extraction_var_surf(nom_fichier_odb,largeur_cyl,largeur_plan,profondeur,tailleMaille,DossierDatFileInc):

    # Fonction d'extraction des donnees qui nous interessent dans l ODB abaqus
    odb = openOdb(path=str(nom_fichier_odb))		# on recupere l odb
    myAssembly = odb.rootAssembly					# on recupere l'assembly
	
    print "\nDEBUT FCT EXCTRACTION\n"
	
    nbr_frames = len(odb.steps['Fretting_fatigue'].frames)		# donne le nombre de frame (commence a zero)
    nbr_elem_par_box_plan     = int(round((largeur_plan/tailleMaille)*(profondeur/tailleMaille)))
    nbr_noeuds_setContactPlan = int(round(largeur_plan /tailleMaille+1))
    nbr_noeuds_setContactCyl  = int(round(largeur_cyl  /tailleMaille+1)) 
    
    print "Parametres d entree :\n"
    
    print "largeur_plan = %f\n" % largeur_plan
    print "largeur_cyl  = %f\n" % largeur_cyl
    print "profondeur   = %f\n" % profondeur
    print "tailleMaille = %f\n" % tailleMaille
    print "nbr_elem_par_box_plan     = %f\n" % nbr_elem_par_box_plan
    print "nbr_noeuds_setContactPlan = %f\n" % nbr_noeuds_setContactPlan
    print "nbr_noeuds_setContactCyl  = %f\n" % nbr_noeuds_setContactCyl
    
    

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


## RECUPERATION DES NSETS UTILES POUR L'EXTRACTION  ################################
    setContactPlan = myAssembly.nodeSets['NOEUDS_CONTACT_PLAN']     
    setContactCyl  = myAssembly.nodeSets['NOEUDS_CONTACT_CYLINDRE']
    RP1            = myAssembly.instances['ASSEMBLY'].nodes[0]
    # basePlan=myAssembly.nodeSets['BASE_PLAN']
    boxPlan = myAssembly.elementSets['BOX_PLAN']

	
## EXTRACTION DONNEES DE DEPLACEMENT U1 ET U2 AU POINT DE REFERENCE RP1   ############
    CompleteName = os.path.join(DossierDatFileInc, 'res_delta.dat')
    res          = open(CompleteName, 'w')			# on ouvre un fichier en ecriture
    print "Extraction res_delta"
    print "Nombre d'increments dans la step fretting =",nbIncrStep[3]
    print "frames de la step",stepFret.frames
    incr = 0
    for incr in range(nbIncrStep[3][0]):
        Uvaluesdisp = stepFret.frames[incr].fieldOutputs['U'].getSubset(region=RP1).values
        delta       = Uvaluesdisp[0].data[0]
        res.write('%0.8f\n' % delta)
    res.close()


## EXTRACTION DONNEES DE FORCE DE REACTION Q (RF1 dans RF) AU POINT DE REFERENCE RP1 ##
    CompleteName = os.path.join(DossierDatFileInc, 'res_Q.dat')
    Q    = open(CompleteName, 'w')
    incr = 0
    print "Extraction res_Q"
    print "\n\n" 
    for incr in range(nbIncrStep[3][0]):
        RFValRP = stepFret.frames[incr].fieldOutputs['RF'].getSubset(region=RP1).values
        Q.write('%0.8f\n' % (RFValRP[0].data[0]))
    Q.close()
    
#####################################
##     ____                    _   ##
##    / ___|___   ___  _ __ __| |  ##
##   | |   / _ \ / _ \| '__/ _` |  ##
##   | |__| (_) | (_) | | | (_| |  ##
##    \____\___/ \___/|_|  \__,_|  ##
##                                 ##
#####################################
##
### COORDONNEES : COORD (X et Y) POUR LE NSET NCP en debut de fretting ##########
    NameCoordCyl  = os.path.join(DossierDatFileInc, 'res_coord_surf_cyl.dat')
    res_coord_cyl = open(NameCoordCyl, 'w')
    print "\nExtraction res_coord_surf_cyl"
    CoordValContact = stepFret.frames[0].fieldOutputs['COORD'].getSubset(region=setContactCyl).values
    coord = zeros((nbr_noeuds_setContactCyl,2))*1.
    i=0
    for noeud_i in CoordValContact:
        coord[i][0] = noeud_i.data[0]
        coord[i][1] = noeud_i.data[1]
        i+=1

    for j in range(int(nbr_noeuds_setContactCyl)):
        res_coord_cyl.write('%0.8f\t%0.8f\n' % (coord[j][0],coord[j][1]))
    res_coord_cyl.close()

    
    
    
    
    NameCoordPlan  = os.path.join(DossierDatFileInc, 'res_coord_surf_plan.dat')	
    res_coord_plan = open(NameCoordPlan, 'w')
    
    print "\nExtraction res_coord_surf_plan"
    CoordValContact = stepFret.frames[0].fieldOutputs['COORD'].getSubset(region=setContactPlan).values
    coord = zeros((nbr_noeuds_setContactPlan,2))*1.
    i = 0
    for noeud_j in CoordValContact:
        coord[i][0] = noeud_j.data[0]
        coord[i][1] = noeud_j.data[1]
        i += 1 

    for j in range(int(nbr_noeuds_setContactPlan)):
        res_coord_plan.write('%0.8f\t%0.8f\n' % (coord[j][0],coord[j][1]))
    res_coord_plan.close()
  
	
	
###############################################
##   _____                   _               ##
##  |  __ \                 (_)              ##
##  | |__) | __ ___  ___ ___ _  ___  _ __    ##
##  |  ___/ '__/ _ \/ __/ __| |/ _ \| '_ \   ##
##  | |   | | |  __/\__ \__ \ | (_) | | | |  ##
##  |_|   |_|  \___||___/___/_|\___/|_| |_|  ##
##                                           ##
###############################################         
## CPRESS POUR LES NSET  NCC  ET  NCP  ##########################

    NamePressCyl = os.path.join(DossierDatFileInc, 'res_pression_cyl.dat')
    resPressCyl  = open(NamePressCyl, 'w')
    print "\nExtraction res_pression_cyl"
    CPRESSValContact = stepFret.frames[10].fieldOutputs['CPRESS'].getSubset(region=setContactCyl).values
    i = 0
    for noeud_i in CPRESSValContact:
        resPressCyl.write('%15.8f\n' % (noeud_i.data))
        i += 1
    resPressCyl.close()
	
    
    
    NamePressPlan = os.path.join(DossierDatFileInc, 'res_pression_plan.dat')		
    resPressPlan  = open(NamePressPlan, 'w')
    print "\nExtraction res_pression_plan"
    CPRESSValContact = stepFret.frames[10].fieldOutputs['CPRESS'].getSubset(region=setContactPlan).values
    j=0
    for noeud_j in CPRESSValContact:
        resPressPlan.write('%15.8f\n' % (noeud_j.data))
        j+=1
    resPressPlan.close()



##########################################################################################
##    ____     _                        ____     _ _         _   _                      ##
##   / ___|___| |__   ___  __ _ _ __   / ___|___| (_)_ __   | | | |   ___  _ __   ___   ##
##  | |   / __| '_ \ / _ \/ _` | '__| | |   / __| | | '_ \  | | | |  / _ \| '_ \ / _ \  ##
##  | |___\__ \ | | |  __/ (_| | |    | |___\__ \ | | |_) | | |_| | | (_) | | | |  __/  ##
##   \____|___/_| |_|\___|\__,_|_|     \____|___/_|_| .__/   \___/   \___/|_| |_|\___|  ##
##                                                  |_|                                 ##
##########################################################################################
##
## CSHEAR1 pour les NSET NCP et NCC                                   ##
## CSLIP1 pour les NSET NCP et NCC                                    ##
## DEFORMATION : U (U1=X et U2=Y) pour les NSET NCP et NCC (fretting) ##

    NameShearCyl = os.path.join(DossierDatFileInc, 'res_cisaillement_cyl.dat')
    NameSlipCyl  = os.path.join(DossierDatFileInc, 'res_cdisp_cyl.dat')
    NameDeplaCyl = os.path.join(DossierDatFileInc, 'res_deplac_surf_cyl.dat')	
        
    res_cisail_cyl = open(NameShearCyl, 'w')
    res_gliss_cyl  = open(NameSlipCyl , 'w')
    res_deplac_cyl = open(NameDeplaCyl, 'w')
    
    print "\nExtraction sur frames : Cyl"
    ii = 0
    jj = 0
    # Stockage des donnees recuperees dans une matrice une ligne correspond a un noeud, et une colonne correspond a un increment
    M_CSHEAR_Cyl = zeros((nbr_noeuds_setContactCyl,nbr_frames))*1.
    M_CSLIP_Cyl  = zeros((nbr_noeuds_setContactCyl,nbr_frames))*1.
    M_U_Cyl      = zeros((nbr_noeuds_setContactCyl,nbr_frames))*1.
    
    for incr in stepFret.frames : 
        CSHEAR_NCC_incr_jj = incr.fieldOutputs['CSHEAR1'].getSubset(region=setContactCyl).values
        CSLIP_NCC_incr_jj  = incr.fieldOutputs['CSLIP1' ].getSubset(region=setContactCyl).values
        U_NCC_incr_jj      = incr.fieldOutputs['U'      ].getSubset(region=setContactCyl).values

        for noeud_ii in CSHEAR_NCC_incr_jj :
            M_CSHEAR_Cyl[ii][jj] = noeud_ii.data
            ii += 1
        ii = 0
        for noeud_ii in CSLIP_NCC_incr_jj :
            M_CSLIP_Cyl[ii][jj]  = noeud_ii.data
            ii += 1
        ii = 0
        for noeud_ii in U_NCC_incr_jj :
            M_U_Cyl[ii][jj] = noeud_ii.data[0]
            ii += 1
        ii  = 0
        jj += 1
	# Enregistrement des donnees de la matrice dans le fichier .dat
    for ii in range(int(nbr_noeuds_setContactCyl)) :
        for jj in range(int(nbr_frames)) :
            res_cisail_cyl.write('%15.8f   ' % M_CSHEAR_Cyl[ii][jj])
            res_gliss_cyl.write( '%15.8f   ' % M_CSLIP_Cyl[ ii][jj])
            res_deplac_cyl.write('%15.8f   ' % M_U_Cyl[     ii][jj])
        res_cisail_cyl.write('\n')
        res_gliss_cyl.write( '\n')
        res_deplac_cyl.write('\n')

    res_cisail_cyl.close()
    res_gliss_cyl.close()
    res_deplac_cyl.close()





    NameShearPlan = os.path.join(DossierDatFileInc, 'res_cisaillement_plan.dat')
    NameSlipPlan  = os.path.join(DossierDatFileInc, 'res_cdisp_plan.dat')
    NameDeplaPlan = os.path.join(DossierDatFileInc, 'res_deplac_surf_plan.dat')		
        
    res_cisail_plan = open(NameShearPlan, 'w')
    res_gliss_plan  = open(NameSlipPlan , 'w')
    res_deplac_plan = open(NameDeplaPlan, 'w')
    print "\nExtraction sur frames : Plan "
    ii = 0
    jj = 0
    # Stockage des donnees recuperees dans une matrice : une ligne correspond a un noeud, et une colonne correspond a un increment
    M_CSHEAR_Plan = zeros((nbr_noeuds_setContactPlan,nbr_frames))*1.
    M_CSLIP_Plan  = zeros((nbr_noeuds_setContactPlan,nbr_frames))*1.
    M_U_Plan      = zeros((nbr_noeuds_setContactPlan,nbr_frames))*1.
    for incr in stepFret.frames : 
        CSHEAR_NCP_incr_jj = incr.fieldOutputs['CSHEAR1'].getSubset(region=setContactPlan).values
        CSLIP_NCP_incr_jj  = incr.fieldOutputs['CSLIP1' ].getSubset(region=setContactPlan).values
        U_NCP_incr_jj      = incr.fieldOutputs['U'      ].getSubset(region=setContactPlan).values
        for noeud_ii in CSHEAR_NCP_incr_jj :
            M_CSHEAR_Plan[ii][jj] = noeud_ii.data
            ii += 1
        ii = 0 
        for noeud_ii in CSLIP_NCP_incr_jj :
            M_CSLIP_Plan[ ii][jj] = noeud_ii.data
            ii += 1
        ii = 0
        for noeud_ii in U_NCP_incr_jj :
            M_U_Plan[     ii][jj] = noeud_ii.data[0]
            ii += 1
        ii  = 0
        jj += 1
	# Enregistrement des donnees de la matrice dans le fichier .dat
    for ii in range(int(nbr_noeuds_setContactPlan)) :
        for jj in range(int(nbr_frames)) :
            res_cisail_plan.write('%15.8f   ' % M_CSHEAR_Plan[ii][jj])
            res_gliss_plan.write( '%15.8f   ' % M_CSLIP_Plan[ ii][jj])
            res_deplac_plan.write('%15.8f   ' % M_U_Plan[     ii][jj])
        res_cisail_plan.write('\n')
        res_gliss_plan.write('\n')
        res_deplac_plan.write('\n')

    res_cisail_plan.close()
    res_gliss_plan.close()
    res_deplac_plan.close()	


