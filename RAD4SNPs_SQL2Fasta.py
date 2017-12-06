#!/usr/bin/python



########################################################################## RAD4SNPs:############################################################################## 
# A set of Python scripts to select and validate independent SNPs markers from a list of read files #
##################################################################################################################################################################

# SUB-PROGRAM FOR FIRST FILTERING FROM MYSQL DATABASE


# Authors: G.LASSALLE (gilles.lassalle@inra.fr) & C.DELORD (chrystelle.delord@inra.fr)
# Last update: AUGUST 2017


#################### PRE-CONDITIONS 
#- [-i] File 
#- Connexion to the Stacks MySQL database available: databases of Stacks 'denovo_map output' for each species.




import argparse
import os
import csv
import sys
import MySQLdb

parser = argparse.ArgumentParser()
# General arguments
parser.add_argument('-o', action='store', dest='OutputFile', help='Fichier Fasta de sortie') # Name of the candidate sequences file
parser.add_argument('-d', action='store', dest='Database', help='Stacks database')
parser.add_argument('-c', action='store', dest='CodeSp', help='ID of the species')
# SQL request arguments
parser.add_argument('-dMin', action='store', dest='depthAlleleMin', help='')
parser.add_argument('-dMax', action='store', dest='depthAlleleMax', help='')
parser.add_argument('-p5', action='store', dest='pos5prim', help='')
parser.add_argument('-p3', action='store', dest='pos3prim', help='')
parser.add_argument('-q', action='store', dest='numSamplesHomologous', help='')
parser.add_argument('-p', action='store', dest='bothSamplesPolymorphic', help='Default "TRUE")

parser.add_argument('--version', action='version', version='%(prog)s 0.1')
results = parser.parse_args()


if results.depthAlleleMin:
	dMin=int(results.depthAlleleMin)
else:    
	dMin=5

if results.depthAlleleMax:
	dMax=int(results.depthAlleleMax)
else:    
	dMax=100
 
if results.pos5prim:
	p5=int(results.pos5prim)
else:    
	p5=30 
 
if results.pos3prim:
	p3=int(results.pos3prim)
else:    
	p3=110
 
if results.nbMaxAlleles:
	x=int(results.nbMaxAlleles)
else:    
	x=2
 
if results.numSamplesHomologous:
	q=int(results.numSamplesHomologous)
else:    
	q=2
 
if results.bothSamplesPolymorphic:
	p=str("\'aa:1(50.0%);ab:1(50.0%);\'")
else:    
	p=str("\'ab:2(100.0%);\'")


#########################################
# Server connexion
########################################
DBHost="127.0.0.1"# localhost 


##############################################################################
# Output file checking (.fasta file in text format for the listing of candidate SNPs)
##############################################################################

if results.OutputFile:
    try:
        FastaFileOut=open(results.OutputFile,"w")
        print "Creation of candidate SNPs fasta file: done."
    except:
        print "Creation of candidate SNPs fasta file has failed."
        print("Unexpected error:", sys.exc_info()[0])   
        sys.exit()           
else :
    print "Please check out inserted path for Output file [-o] of candidate sequences from MySQL. End of program."
    sys.exit()

    
##############################################################################
# Database connexion checking
##############################################################################

print "Initialisation of Stacks MySQL database" # normalement deja fait avec le prog principal! a valider avec Gilles.
if results.Database:
    db = MySQLdb.connect(host="",          # your host, usually localhost
        user="",                           # your username
        passwd="",                         # your password
        db=results.Database)               # name of the data base
    cur1= db.cursor()                      # connexion
else:
    print "Database not found, please check out inserted ID for the database [-d]. End of program."
    sys.exit()  
    
    
if results.CodeSp:    
    CodeEspece=str(results.CodeSp)    
else:
    CodeEspece=""
print  "The species ID that will be used is: "+str(CodeEspece)    



##############################################################################
# MySQL request
##############################################################################

print "Retrieving candidate markers from Stacks MySQL database..."

Request="select distinct ma.catalog_id,col,rank_1,rank_2,seq from matches ma join catalog_snps cs on ma.catalog_id=cs.tag_id join catalog_tags ct on ma.catalog_id=ct.tag_id join catalog_index ci on cs.tag_id=ci.tag_id where ma.catalog_id not in (select catalog_id FROM matches where depth < "+str(dMin)+") and ma.catalog_id not in (select catalog_id FROM matches where depth > "+str(dMax)+") and cs.col > "+str(p5)+" and cs.col < "+str(p3)+" and ci.snps = 1 and ci.alleles = 2 and ci.parents = "+str(q)+" and ci.ratio = "+str(p)+" order by ma.depth;"  
print Request


###############################################################################
# Preparation to write candidate markers information in the output .fasta file
###############################################################################

def AjoutFasta(FastaFile,IdSeq,Seq,code):
    IdSeq=">"+str(CodeEspece)+str(IdSeq)
    FastaFile.write("%s\n" % (IdSeq))
    FastaFile.write("%s\n" %	(Seq))   
    

###############################################################################
# Running the sub-program for first filtering
###############################################################################


cur1.execute(Request) 
# Print all the first cell of all the rows

for row in cur1.fetchall():
    IdSeq=str(row[0])+"|"+str(row[1])+"|["+str(row[2])+"/"+str(row[3])+"]"
    Seq=str(row[4])
    print "idseq="+str(IdSeq)
    print "seq="+str(Seq)
    
    AjoutFasta(FastaFileOut,IdSeq,Seq,CodeEspece) # Writing candidate markers information in the output .fasta file
db.close()
