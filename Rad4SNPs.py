#!/usr/bin/python

# PRE-COND : 
#1 - extraction de la base de donnees stacks ==> fichier fasta non vide
#2- repertoire -o out accessible en ecriture
#3- bwa installe
#4- samtools installe

# EXECUTION
# a partir du fichier fasta on lance toutes les commandes BWA et pileup
# et tous les filtres bwa samtools et pileup
# extraction finale de la dbb


# -i input : Fasta en sortie  RAD4SNPSSqlextract
# -i1 input : fichier de Read1
# -i2 input : fichier de read1

 
# POST-COND : multi fichiers fasta de read 2 au nom du catalog_id
#G.Lassalle&C.DELORD
#june 2017
###############################################################################
import argparse
import os
import sys
import MySQLdb
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='InputDir', help='Working Directory') #Liste de repertoires ?
parser.add_argument('-d', action='store', dest='database', help='Database Stacks') #Fichier texte avec liste de databases ?
parser.add_argument('-c', action='store', dest='CodeSp', help='code espece') #Fichier texte avec liste de codes Sp ?
parser.add_argument('-i1', action='store', dest='R11', help='first R1 file')
parser.add_argument('-i2', action='store', dest='R12', help='Second R1 file')
parser.add_argument('-i3', action='store', dest='R21', help='Third R1 file')
parser.add_argument('-i4', action='store', dest='R22', help='Fourth R1 file')
#parser.add_argument('-i5', action='store', dest='R21', help='first R2 file')
#parser.add_argument('-i6', action='store', dest='R22', help='Second R2 file')
#parser.add_argument('-i7', action='store', dest='R211', help='Third R2 file')
#parser.add_argument('-i8', action='store', dest='R221', help='Fourth R2 file')


parser.add_argument('--version', action='version', version='%(prog)s 0.1')

results = parser.parse_args()
print 'input directory		=', results.InputDir



##############################################################################
# On teste les arguments
##############################################################################    
#argument 1
#on verifie le infile + Initialisation des OUTPUT
if results.InputDir:
    if os.path.isdir(results.InputDir):
        print " repertoire de travail valide."
        
    else :
        print "Repertoire Input non trouve Fin du programme."
        sys.exit()
    
else :
    print "-i non renseigne => fin du programme"
    sys.exit()
    
#argument 2 et 3
#on verifie le infile + Initialisation des OUTPUT
if results.R11:
    if os.path.isfile(results.R11):
        print "premier Fichier R1 renseigne ...!!!" 
    else :
        print "Ceci n'est pas un fichier => Fin du programme."
        sys.exit()
        
else :
    print "fichier de Reads 1 non renseigne => fin du programme"
    sys.exit()
    

##############################################################################
#on ne teste pas le second fichier de R1
if results.database:
    db = MySQLdb.connect(host="127.0.0.1",      # your host, usually localhost
        user="root",                                # your username
        passwd="rif-35",                            # your password
        db=results.database)                        # name of the data base
    cur1= db.cursor()                               #initialisation de la connexion
else:
    print "Pas de base de donnees => fin du programme"
    sys.exit()  
###############################################################################    
    
if results.CodeSp:    
    CodeEspece=str(results.CodeSp)
    if CodeEspece[:1]!="_":
        CodeEspece=str(results.CodeSp)+str("_")
else:
    CodeEspece="std_"
        
#argument 4
"""        
#on verifie le OUTPUT
if results.OutDir:
    if os.path.isdir(results.OutDir):
        print "Repertoire de sortie ok ...!!!" 
        WorkDir=os.path.abspath(results.OutDir)
    else :
        print " Repertoire non valide => Fin du programme."
        sys.exit()        
else :
    print "-o non renseigne => fin du programme"
    sys.exit()
"""    

########################## Variables environnement #############################
WorkDir=os.path.abspath(results.InputDir)
FastaCatalog=str(WorkDir)+"/"+str(results.CodeSp)+"Catalog.fasta"

     
###############################################################################
#prog principal
###############################################################################
     
####
#on teste les programmes BWA et samtools     
if os.path.isfile("/usr/bin/bwa"):
    print "BWA ...ok" 
    
else :    
    print " BWA non trouve => Fin du programme."
    sys.exit()      
###samtools
if os.path.isfile("/usr/bin/samtools"):
    print "samtools ...ok" 
    
else :    
    print " samtools non trouve => Fin du programme."
    sys.exit()      

#####################################################
#Working Dir Writable

filepath = results.InputDir+'/file.txt'

try:
    filehandle = open( filepath, 'w' )
except IOError:
    sys.exit( 'Repertoire de travail non accessible ' + filepath )





##############################################################################
#Commandes !!!!!!!!!

###############################################################################
#fusion des fichiers de reads1 si ils sont plusieurs
#print "Repertoire de travail : "+str(WorkDir) 

print os.getcwd()

commandeExtractFasta="./Rad4SnpsSQL2Fasta.py -o "+str(FastaCatalog)+" -d "+str(results.database)+" -c "+str(CodeEspece) 

print "Extraction du fichier fasta" 
print commandeExtractFasta
os.system(commandeExtractFasta)


##############################Fusion de reads 1################################
if results.R11:
    if results.R12:
        commandeFusion="cat "+str(results.R11)+" "+str(results.R12)+" > "+str(WorkDir)+"/allR1.fq.gz"
        
    else :
        commandeFusion="cp "+str(results.R11)+" "+str(WorkDir)+"/allR1.fq.gz"
#############################fin de fusion       
        
##############################Fusion de reads 2################################
if results.R21:
    if results.R22:
        commandeFusion1="cat "+str(results.R21)+" "+str(results.R22)+" > "+str(WorkDir)+"/allR2.fq.gz"
        
    else :
        commandeFusion1="cp "+str(results.R21)+" "+str(WorkDir)+"/allR2.fq.gz"
#############################fin de fusion   
    
#################################### BWA Etape 1 ############################## 
#creation de l'index BWA
commande1="bwa index "+str(FastaCatalog)

commande2="bwa mem -a -M "+str(FastaCatalog)+" "+str(WorkDir)+"/allR1.fq.gz > "+str(WorkDir)+"/PremierAlign.sam"

commande3="samtools view -Sb "+str(WorkDir)+"/PremierAlign.sam | samtools sort - "+str(WorkDir)+"/PremierAlign1Sorted "

commande4="samtools view -F4 "+str(WorkDir)+"/PremierAlign1Sorted.bam > "+str(WorkDir)+"/PremierAlign1Sorted-F4.sam"


#commandex="samtools mpileup -vcf "+str(HelixAlign1Sorted-F4R1Filtered.fa) +" "+str(HelixAlignR1FilteredSorted.bam) +" > "+str(Candidats2.pileup)


print "Execution de "+str(commandeFusion)
os.system(commandeFusion)

print "Execution de "+str(commandeFusion1)
os.system(commandeFusion1)

print "Execution de "+str(commande1)
os.system(commande1)

print "Execution de "+str(commande2)
os.system(commande2)

print "Execution de "+str(commande3)
os.system(commande3)

print "Execution de "+str(commande4)
os.system(commande4)

print "     ************************************************************************"
print "Execution du script python Rad4SnpsR1SamFilterV2.py avec les parametres par defauts"
print "     ************************************************************************"

print os.getcwd()

commande5="./Rad4SnpsR1SamFilterV2.py -i "+str(WorkDir)+"/PremierAlign1Sorted-F4.sam"
os.system(commande5)


######################### FIN ETAPE ###########################################

Candidatfasta1=str(WorkDir)+"/PremierAlign1Sorted-F4R1Filtered.fa"

if os.path.isfile(Candidatfasta1):
    print "Fichier candidat Fasta cree ....passage a l etape 2"
else :
    sys.exit( '****ERROR**** fichier fasta non cree verifiez les etapes')


#### ON RECOMMENCE ##################
#creation de l'index BWA
commande21="bwa index "+str(Candidatfasta1)

commande22="bwa mem -a -M "+str(Candidatfasta1)+" "+str(WorkDir)+"/allR1.fq.gz > "+str(WorkDir)+"/SecondAlign.sam"

commande23="samtools view -Sb "+str(WorkDir)+"/SecondAlign.sam | samtools sort - "+str(WorkDir)+"/SecondAlign1Sorted"

#commande24="samtools view -F 4 "+str(WorkDir)+"/SecondAlign1Sorted.bam > "+str(WorkDir)+"/SecondAlign1Sorted-F4.sam"
#commande24 : deprecated

commande25="samtools index "+str(WorkDir)+"/SecondAlign1Sorted.bam"

commande25bis="samtools faidx "+str(Candidatfasta1)

commande26="samtools mpileup -d 1000 -O --ff 4 -f "+str(Candidatfasta1) +" "+ str(WorkDir)+"/SecondAlign1Sorted.bam"+" > "+str(WorkDir)+"/Candidats2.pileup"

print "Execution de "+str(commande21)
os.system(commande21)

print "Execution de "+str(commande22)
os.system(commande22)

print "Execution de "+str(commande23)
os.system(commande23)

#print "Execution de "+str(commande24)
#os.system(commande24)
print "Execution de "+str(commande25)
os.system(commande25)
os.system(commande25bis)

print "Execution de "+str(commande26)
os.system(commande26)


##################################################################################
#Execution RAd4Snps pileup filter 
##################################################################################
print os.getcwd()

commande27="./Rad4SnpsR2PileupFilter.py -i "+str(WorkDir)+"/Candidats2.pileup"
print "Execution de "+str(commande27)
os.system(commande27)

commande28="./Rad4SnpsFinalSqlExtractV2.py -i"+str(WorkDir)+"/Candidats2NoMulti.txt -d "+str(results.database)+" -c "+str(CodeEspece)+" > "+str(WorkDir)+"/CandidatFin.fasta"
print "Execution de "+str(commande28)
os.system(commande28) 

commande28bis="sed -i '1d' "+str(WorkDir)+"/CandidatFin.fasta"
os.system(commande28bis) 

#### ON RECOMMENCE AVEC LES READS2#################
#creation de l'index BWA

CandidatFin=str(WorkDir)+"/CandidatFin.fasta"

commande29="bwa index "+str(CandidatFin)

commande30="bwa mem -a -M "+str(CandidatFin)+" "+str(WorkDir)+"/allR2.fq.gz > "+str(WorkDir)+"/ThirdAlign.sam"

commande31="samtools view -Sb "+str(WorkDir)+"/ThirdAlign.sam | samtools sort - "+str(WorkDir)+"/ThirdAlign2Sorted"

commande32="samtools index "+str(WorkDir)+"/ThirdAlign2Sorted.bam"

commande32bis="samtools faidx "+str(CandidatFin)

commande33="samtools mpileup -d 1000 -O --ff 4 -f "+str(CandidatFin)+" "+str(WorkDir)+"/ThirdAlign2Sorted.bam"+" > "+str(WorkDir)+"/Candidats3.pileup"

print "Execution de "+str(commande29)
os.system(commande29)

print "Execution de "+str(commande30)
os.system(commande30)

print "Execution de "+str(commande31)
os.system(commande31)

print "Execution de "+str(commande32)
os.system(commande32)
os.system(commande32bis)

print "Execution de "+str(commande33)
os.system(commande33)

##################################################################################
#Execution RAd4Snps pileup filter de nouveau
##################################################################################
print os.getcwd()

commande34="./Rad4SnpsR2PileupFilter.py -i "+str(WorkDir)+"/Candidats3.pileup"
print "Execution de "+str(commande34)
os.system(commande34)

commande35="./Rad4SnpsFinalSqlExtractV2.py -i"+str(WorkDir)+"/Candidats3NoMulti.txt -d "+str(results.database)+" -c "+str(CodeEspece)+" > "+str(WorkDir)+"/SNPs_out.fasta"
print "Execution de "+str(commande35)
os.system(commande35) 