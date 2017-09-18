#!/usr/bin/python



########################################################################## RAD4SNPs:############################################################################## 
# A set of Python scripts to select and validate independent SNPs markers from a list of read files #
##################################################################################################################################################################

# SUB-PROGRAM FOR SECOND FILTERING PART 2, AND THIRD FILTERING : Filtering on a .pileup file for SE reads then for PE reads.


# Authors: G.LASSALLE (gilles.lassalle@inra.fr) & C.DELORD (chrystelle.delord@inra.fr)
# Last update: AUGUST 2017


#################### PRE-CONDITIONS 
#- [-i] File 


###############################################################################
import argparse
import os
import csv
import sys
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='InputFile', help=' sorted SAM file')
parser.add_argument('-nv', action='store', dest='depthFlankMax', help='Couverture maximale autorisee variant flanquant (def=2)')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
results = parser.parse_args()


print 'inputfile		=', results.InputFile


##############################################################################
# Input file checking
##############################################################################  

if results.InputFile:
    if os.path.isfile(results.InputFile):
        print "Inserted input is a file: OK." 
        ReadCsvFile=csv.reader(open(results.InputFile,"r"),delimiter='\t')
        Base, BaseExt = os.path.splitext(results.InputFile)
    # tester acces ecriture        
        GoodFileOut=str(Base)+"NoMulti.txt"
        FileIn0=open(GoodFileOut,"w")
        BadFileOut=str(Base)+"WithMulti.txt"
        FileIn1=open(BadFileOut,"w")

        if BaseExt == ".pileup" :
            print ".pileup extension: OK."
        else:
            print "The intput file extension is not .pileup. End of program." 
            sys.exit()
    else :
        print "Input file is not found. End of program."
        sys.exit()
        
else :
    print "Please check out inserted path for the .pileup input file [-i]. End of program."
    sys.exit()


if results.depthFlankMax :
    depthFlankMax=int(results.depthFlankMax)
else:
    depthFlankMax=2
    

##############################################################################
# Sub-programs for analyse of input .pileup file
##############################################################################

def testPosition(lignelue):
   
    idSeqCurrent=str(lignelue[0])  # Identifiant of the candidate sequence.
    PosSnp= (idSeqCurrent.split('|'))[1] 
    Poslue=int(lignelue[1]) 
    Variant=str(lignelue[4])
    CompteurVariant=0
    if int(Poslue)<>int(PosSnp)+1: # In Stacks, bp positions start at zero. A SNP given at the 8th position according to Stacks is actually the 9th base pair.
        for Nuc in Variant:
            if Nuc.lower() in listNuc:
                CompteurVariant+=1
                
    if int(CompteurVariant) >= int(depthFlankMax) :
        print str(idSeqCurrent)+" keep ="+str(keep)          
        return 0
        # 0 = Discard of candidate sequences where a flanking variant is found with depth coverage higher than depthFlankMax.
    
###############################################################################    
    
def AjoutLineTxt(Fichier, IdSeq):
   
   Fichier.write("%s%s\n" % (">",IdSeq)) # Sorting candidate sequences between those that passed second filtering (2) or third filtering / those that did not.
    
    
        
###############################################################################
# Intialisation of variables        
###############################################################################          

# Variables for a given read aligned against the candidate sequence
Variant=''       # sequence du read 

# Variables for alignments on a candidate sequence
idSeqRef=0      # ID of the analyzed candidate
keep=1          # The alignment/current candidate sequence is kept by default   

# Variables for the whole sub-program
SeqTotal=0      # Counting number of total candidate sequences  
SeqGood=0       # Counting number of SE validated or PE-validated candidate sequences
SeqWrong=0      # Counting number of discarded candidate sequences


###############################################################################

print "Beginning to read the .pileup file..."        

# Initialisation:

idSeqRef=''
listNuc=('a','t','c','g')       
CompteurPileup=0

for refligne in ReadCsvFile:    
    
    idSeqCurrent=str(refligne[0])                     
    if idSeqRef=='':
        # At the beginning of the .pileup file screening, idSeqRef is empty.
        idSeqRef=idSeqCurrent 
        if testPosition(refligne)==0:keep=0
        
                
      
    else:
        if idSeqCurrent==idSeqRef:
            if testPosition(refligne)==0:keep=0
        
        
        else :  # Switching to the next candidate sequence alignment:
            

            if keep==1:
                AjoutLineTxt(FileIn0, idSeqRef) # Listing the IDs of candidate sequences that passed the filter
                         
            if keep==0:                            
                AjoutLineTxt(FileIn1, idSeqRef) # Listing the IDs of candidate sequences that did not passed the filter (i.e. discarded sequences)
              
                
          
# Re-Initialisation:
          
            idSeqRef=idSeqCurrent 
            keep=1
            CompteurPileup+=1
         
            if testPosition(refligne)==0:keep=0

                                                  
#########################################################################                
        # .pileup file screening is complete: analyzing the last candidate.    
########################################################################  
                                                  

if keep==1:
    AjoutLineTxt(FileIn0, idSeqRef)                        
else :
    idSeqRef=idSeqCurrent                 
    AjoutLineTxt(FileIn1, idSeqRef)
    
CompteurPileup+=1
print "Number of analyzed candidates = "+str(CompteurPileup)


# End.
