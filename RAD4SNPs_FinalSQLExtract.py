#!/usr/bin/python



########################################################################## RAD4SNPs:############################################################################## 
# A set of Python scripts to select and validate independent SNPs markers from a list of read files #
##################################################################################################################################################################

# SUB-PROGRAM FOR EXTRACTING SUCCESSFUL CANDIDATE SEQUENCES (i.e. that passed the pileup filtering) FROM STACKS DATABASE
# En realite, il serait plus logique de simplement recuperer ces bonnes sequences dans le fichier fasta ayant servi de base a la construction du pileup.


# Authors: G.LASSALLE (gilles.lassalle@inra.fr) & C.DELORD (chrystelle.delord@inra.fr)
# Last update: AUGUST 2017


#################### PRE-CONDITIONS 
#- [-i] File 

###############################################################################
import argparse
import os
import csv
import sys
import MySQLdb
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='InputFile', help='.txt of candidate sequences that passed the .pileup filter')
parser.add_argument('-d', action='store', dest='Database', help='Stacks database')
parser.add_argument('-c', action='store', dest='CodeSp', help='ID of the species')

parser.add_argument('--version', action='version', version='%(prog)s 0.1')
results = parser.parse_args()


print 'inputfile		=', results.InputFile


##############################################################################
# Input file checking
##############################################################################   

if results.InputFile:
    if os.path.isfile(results.InputFile):
        
        ReadCsvFile=csv.reader(open(results.InputFile,"r"),delimiter='|')
        Base, BaseExt = os.path.splitext(results.InputFile)
        if BaseExt != ".txt" :
            print "The intput file extension is not .txt. End of program." 
            sys.exit() # On teste le format du fichier via son extension, serait mieux avec PySam.
    else :
        print "Input file is not found. End of program."
        sys.exit()
        
else :
    print "Please check out inserted path for the .txt input file of good sequences [-i]. End of program."
    sys.exit()
    
###############################################################################
# INITIATE Mysql connect
###############################################################################
db = MySQLdb.connect(host="",             # your host, usually localhost
                     user="",             # your username
                     passwd="",           # your password
                     db=results.Database) # name of the data base

# you must create a Cursor object. It will let you execute all the queries you need.
cur1= db.cursor()

###############################################################################
# PROCEDURES (a completer)
###############################################################################

def AjoutFasta(FastaFile,IdSeq,Seq):
   
    FastaFile.write("%s\n" % (IdSeq))
    FastaFile.write("%s\n" %	(Seq))

def AjoutSequenom():
    a=5




###############################################################################
# Prog Principal
###############################################################################

for refligne in ReadCsvFile:
    
    idSeq=refligne[0]
    idSeq1=idSeq[1:].split('_')    
    idSeq=idSeq1[1]
    
    requete2= "select distinct ma.catalog_id,col,rank_1,rank_2,seq from matches ma join catalog_snps cs on ma.catalog_id=cs.tag_id join catalog_tags ct on ma.catalog_id=ct.tag_id where ma.catalog_id=\""+idSeq+"\";"
    cur1.execute(requete2)

    for row in cur1.fetchall():
        IdSeq=">"+str(results.CodeSp)+str(row[0])+"|"+str(row[1])+"|["+str(row[2])+"/"+str(row[3])+"]"
        Seq=str(row[4])
        # Generating a fasta file of candidate sequences that passed the .pileup filter.
        print IdSeq          
        print Seq

        
        # PreSeq=Seq[0:int(row[1])]
        # PostSeq=Seq[(int(row[1])+1):]              
        
        # Generating the list of successful candidates sequences in a shape of a sequence with the SNP indicated between brackets [].       
        # MassARRAYId=str(results.CodeSp)+str(row[0])+"\t"+PreSeq+"["+str(row[2])+"/"+str(row[3])+"]"+PostSeq
        # print MassARRAYId
        
       
        
db.close()






