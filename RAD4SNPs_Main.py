#!/usr/bin/python



########################################################################## RAD4SNPs:############################################################################## 
# A set of Python scripts to select and validate independent SNPs markers from a list of read files #
##################################################################################################################################################################

# MAIN PROGRAM


# Authors: G.LASSALLE (gilles.lassalle@inra.fr) & C.DELORD (chrystelle.delord@inra.fr)
# Last update: AUGUST 2017


#################### PRE-CONDITIONS 
#- [-i] Working directory where to store results of the pipeline for the focal species X
#- [-d] exact name of MySQL database where denovo_map.pl Stacks data are available for the focal species X
#- [-i1] single-end reads (reads 1) for focal species X duplicate 1
#- [-i2] single-end reads (reads 1) for focal species X duplicate 2
#- [-i3] paired-end reads (reads 2) for focal species X duplicate 1
#- [-i4] paired-end reads (reads 2) for focal species X duplicate 2

#- BWA and SAMtools available
#- Connexion to the Stacks MySQL database available: databases of Stacks 'denovo_map output' for each species.

 

###############################################################################
import argparse
import os
import sys
import MySQLdb
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='InputDir', help='Working Directory') 
parser.add_argument('-d', action='store', dest='database', help='Stacks database')
parser.add_argument('-c', action='store', dest='CodeSp', help='ID of the species')
parser.add_argument('-i1', action='store', dest='R11', help='First R1 file')
parser.add_argument('-i2', action='store', dest='R12', help='Second R1 file')
parser.add_argument('-i3', action='store', dest='R21', help='First R2 file')
parser.add_argument('-i4', action='store', dest='R22', help='Second R2 file')



parser.add_argument('--version', action='version', version='%(prog)s 0.1')

results = parser.parse_args()
print 'input directory		=', results.InputDir



##############################################################################
# Arguments testing
##############################################################################  

if results.InputDir:
    if os.path.isdir(results.InputDir):
        print "Working directory is valid."
        
    else :
        print "Caution: working directory is invalid, please ckeck [-i]."
        sys.exit()
    
else :
    print "Please insert path for working directory [-i]. End of program."
    sys.exit()
    
##############################################################################
if results.database:
    db = MySQLdb.connect(host="",                   # your host, usually localhost
        user="",                                    # your username
        passwd="",                                  # your password
        db=results.database)                        # name of the database
    cur1= db.cursor()                               # connexion
    print "Currently working on MySQL database:      "+str(results.database)
    
else:
    print "Incorrect ID for database: database not found, please check [-d]"
    sys.exit()  
############################################################################### 


    
#
if results.R11:
    if os.path.isfile(results.R11):
        print "First file of single-end reads: found." 
    else :
        print "Path to single-end reads data is not a file: please check out [-i1]."
        sys.exit()
        
else :
    print "Please insert path to single-end read files [-i1]. End of program."
    sys.exit()


#    
if results.R12:
    if os.path.isfile(results.R12):
        print "Second file of single-end reads: found." 
    else :
        print "Path to single-end reads data is not a file: please check out [-i2]."
        sys.exit()
        
else :
    print "Please insert path to single-end read files [-2]. End of program."
    sys.exit() 
    
       
#   
if results.R21:
    if os.path.isfile(results.R21):
        print "First file of paired-end reads: found." 
    else :
        print "Path to paired-end reads data is not a file: please check out [-i3]."
        sys.exit()
        
else :
    print "Please insert path to paired-end read files [-i3]. End of program."
    sys.exit()
    
       
#   
if results.R22:
    if os.path.isfile(results.R22):
        print "Second file of paired-end reads: found." 
    else :
        print "Path to paired-end reads data is not a file: please check out [-i4]."
        sys.exit()
        
else :
    print "Please insert path to paired-end read files [-i4]. End of program."
    sys.exit()    
    
###############################################################################       
if results.CodeSp:    
    CodeEspece=str(results.CodeSp)
    if CodeEspece[:1]!="_":
        CodeEspece=str(results.CodeSp)+str("_")
else:
    CodeEspece="std_"
###############################################################################         


WorkDir=os.path.abspath(results.InputDir) # Current working directory
FastaCatalog=str(WorkDir)+"/"+str(results.CodeSp)+"Catalog.fasta" # Formatting name of candidates fasta file -output of MySQL filtering

     
###############################################################################
# Main program
###############################################################################
     
if os.path.isfile("/usr/bin/bwa"):
    print "BWA program is found." 
    
else :    
    print "Cannot find BWA: please check out pipeline requirements."
    sys.exit()      
###samtools
if os.path.isfile("/usr/bin/samtools"):
    print "SAMtools program is found." 
    
else :    
    print "Cannot find SAMtools: please check out pipeline requirements."
    sys.exit()      
    

#####################################################
# Working directory writable

filepath = results.InputDir+'/file.txt'

try:
    filehandle = open( filepath, 'w' )
except IOError:
    sys.exit( 'Working directory is not accessible' + filepath )

###############################################################################
# Pipeline commands:
###############################################################################



#################################### FIRST FILTERING ##########################
print os.getcwd()

commandeExtractFasta="./RAD4SNPs_SQL2Fasta.py -o "+str(FastaCatalog)+" -d "+str(results.database)+" -c "+str(CodeEspece) 

print "Extraction du fichier fasta" 
print commandeExtractFasta
os.system(commandeExtractFasta)


############################## Fusion of single-end reads #####################
if results.R11:
    if results.R12:
        commandFusionR1="cat "+str(results.R11)+" "+str(results.R12)+" > "+str(WorkDir)+"/allR1.fq.gz"
        
    else :
        commandFusionR1="cp "+str(results.R11)+" "+str(WorkDir)+"/allR1.fq.gz"
#############################fin de fusion       
        
############################## Fusion of paired-end reads #####################
if results.R21:
    if results.R22:
        commandFusionR2="cat "+str(results.R21)+" "+str(results.R22)+" > "+str(WorkDir)+"/allR2.fq.gz"
        
    else :
        commandFusionR2="cp "+str(results.R21)+" "+str(WorkDir)+"/allR2.fq.gz"
    


#################################### SECOND FILTERING (1) #####################

command1="bwa index "+str(FastaCatalog) # Indexing

command2="bwa mem -a -M "+str(FastaCatalog)+" "+str(WorkDir)+"/allR1.fq.gz > "+str(WorkDir)+"/PremierAlign.sam" # SE reads alignment

command3="samtools view -Sb "+str(WorkDir)+"/PremierAlign.sam | samtools sort - "+str(WorkDir)+"/PremierAlign1Sorted" # Conversion to bam file

command4="samtools view -F4 "+str(WorkDir)+"/PremierAlign1Sorted.bam > "+str(WorkDir)+"/PremierAlign1Sorted-F4.sam" # Elimination of unmapped SE reads


print "SE reads merging: "+str(commandFusionR1)
os.system(commandFusionR1)

print "PE reads merging: "+str(commandFusionR2)
os.system(commandFusionR2)

print "BWA indexing: "+str(command1)
os.system(command1)

print "Alignment: "+str(command2)
os.system(command2)

print "Conversion to bam file: "+str(command3)
os.system(command3)

print "Elimination of unmapped SE reads: "+str(command4)
os.system(command4)


print "     ************************************************************************"
print "                  Second filtering (1) with default parameters               "
print "     ************************************************************************"

print os.getcwd()

commande5="./RAD4SNPs_SamFilter.py -i "+str(WorkDir)+"/PremierAlign1Sorted-F4.sam"
os.system(commande5)


Candidatfasta1=str(WorkDir)+"/PremierAlign1Sorted-F4R1Filtered.fa" # Obtention of incomplete SE-validated fasta file

if os.path.isfile(Candidatfasta1):
    print "SE-validated fasta file about to be completed. Re-aligning to complete second filtering."
else :
    sys.exit( '****ERROR**** A problem occurred. Please check out alignment outputs.')


#################################### SECOND FILTERING (2) #####################


command21="bwa index "+str(Candidatfasta1)

command22="bwa mem -a -M "+str(Candidatfasta1)+" "+str(WorkDir)+"/allR1.fq.gz > "+str(WorkDir)+"/SecondAlign.sam"

command23="samtools view -Sb "+str(WorkDir)+"/SecondAlign.sam | samtools sort - "+str(WorkDir)+"/SecondAlign1Sorted"

command25="samtools index "+str(WorkDir)+"/SecondAlign1Sorted.bam"

command25bis="samtools faidx "+str(Candidatfasta1)

command26="samtools mpileup -d 1000 -O --ff 4 -f "+str(Candidatfasta1) +" "+ str(WorkDir)+"/SecondAlign1Sorted.bam"+" > "+str(WorkDir)+"/CandidatsR1.pileup"

print "BWA indexing: "+str(command21)
os.system(command21)

print "Alignment: "+str(command22)
os.system(command22)

print "Conversion to bam file: "+str(command23)
os.system(command23)

print "Indexing of bam file: "+str(command25)
os.system(command25)

print "Indexing for pileup file: "+str(command25bis)
os.system(command25bis)

print "Construction of SE pileup file:  "+str(command26)
os.system(command26)


print "     ************************************************************************"
print "                  Second filtering (2) with default parameters               "
print "     ************************************************************************"

print os.getcwd()

command27="./RAD4SNPs_PileupFilter.py -i "+str(WorkDir)+"/CandidatsR1.pileup"
print "End of second filtering: elimination of flanking variants: "+str(command27)
os.system(command27)

command28="./RAD4SNPs_FinalSQLExtract.py -i"+str(WorkDir)+"/CandidatsR1NoMulti.txt -d "+str(results.database)+" -c "+str(CodeEspece)+" > "+str(WorkDir)+"/CandidatFin.fasta"
print "Complete SE-validated fasta file: "+str(command28) 
os.system(command28)

command28bis="sed -i '1d' "+str(WorkDir)+"/CandidatFin.fasta"
os.system(command28bis) 


#################################### THIRD FILTERING ########################## 


CandidatFin=str(WorkDir)+"/CandidatFin.fasta"

if os.path.isfile(CandidatFin):
    print "SE-validated fasta file is completed. Re-aligning to perform third filtering."
else :
    sys.exit( '****ERROR**** A problem occurred. Please check out alignment and/or pileup outputs.')
    

command29="bwa index "+str(CandidatFin)

command30="bwa mem -a -M "+str(CandidatFin)+" "+str(WorkDir)+"/allR2.fq.gz > "+str(WorkDir)+"/ThirdAlign.sam"

command31="samtools view -Sb "+str(WorkDir)+"/ThirdAlign.sam | samtools sort - "+str(WorkDir)+"/ThirdAlign2Sorted"

command32="samtools index "+str(WorkDir)+"/ThirdAlign2Sorted.bam"

command32bis="samtools faidx "+str(CandidatFin)

command33="samtools mpileup -d 1000 -O --ff 4 -f "+str(CandidatFin)+" "+str(WorkDir)+"/ThirdAlign2Sorted.bam"+" > "+str(WorkDir)+"/Candidats3.pileup"


print "BWA indexing: "+str(command29)
os.system(command29)

print "Alignment: "+str(command30)
os.system(command30)

print "Conversion to bam file: "+str(command31)
os.system(command31)

print "Indexing of bam file: "+str(command32)
os.system(command32)

print "Indexing for pileup file: "+str(command32bis)
os.system(command32bis)

print "Construction of PE pileup file:  "+str(command33)
os.system(command33)


print "     ************************************************************************"
print "                    Third filtering with default parameters                  "
print "     ************************************************************************"

print os.getcwd()

command34="./RAD4SNPs_PileupFilter.py -i "+str(WorkDir)+"/Candidats3.pileup"
print "End of third filtering: elimination of flanking variants: "+str(command34)
os.system(command34)

command35="./RAD4SNPs_FinalSQLExtract.py -i"+str(WorkDir)+"/CandidatsR2NoMulti.txt -d "+str(results.database)+" -c "+str(CodeEspece)+" > "+str(WorkDir)+"/SNPs_out.fasta"
print "Complete PE-validated fasta file: "+str(command35)
os.system(command35) 

# End.
