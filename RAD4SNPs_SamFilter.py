#!/usr/bin/python



########################################################################## RAD4SNPs:############################################################################## 
# A set of Python scripts to select and validate independent SNPs markers from a list of read files #
##################################################################################################################################################################

# SUB-PROGRAM FOR SECOND FILTERING PART 1 : Filtering on a .sam file.


# Authors: G.LASSALLE (gilles.lassalle@inra.fr) & C.DELORD (chrystelle.delord@inra.fr)
# Last update: AUGUST 2017


#################### PRE-CONDITIONS 
#- [-i] A sorted .sam file of aligned single-end reads
#- Autres arguments a valider definitivement avec Gilles




import argparse
import os
import csv
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='InputFile', help='Sorted .sam file from the first alignment of single-end reads')
parser.add_argument('-x', action='store', dest='nbReadSam', help='Maximum number of outlier reads in an alignment against a candidate marker (def=3)')
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
        if BaseExt == ".sam" :
            print ".sam extension: OK."
            # Output files initialisation
            FastaFileName=str(Base)+"R1Filtered.fa"
            FastaFileOut=open(FastaFileName,"w")      
            
        else:
            print "The intput file extension is not .sam. End of program." 
            sys.exit() 
    else :
        print "Input file is not found. End of program."
        sys.exit()
        
else :
    print "Please check out inserted path for the .sam input file [-i]. End of program."
    sys.exit()


if results.nbReadSam :
    nbRead=int(results.nbReadSam)
else:
    nbRead=3


##############################################################################
#                               Constant values
##############################################################################

MaxMismatch=2
# MaxRatio=0.45

##############################################################################
# Sub-programs for analyse of input .sam file
##############################################################################

def testRead(ligneSamLue,MaxMismatch):
    # Input : a line of the .sam file
    # OutPut : The state of a read (outlier or not) + The state of the candidate marker   
    # A candidate with more than >nbReads outliers is discarded.
    # An "outlier" read is a read that does not align end-to-end to the candidate sequence, or display >2 mismatches with the candidate sequence.
 
    ReadState=1 # The state of a read is considered as not-outlier by default.
    NbreNMMax=0 
    idSeqCurrent=str(ligneSamLue[2]) 
    Flag=str(ligneSamLue[1])
    PosSnp=(idSeqCurrent.split('|'))[1]          
    Seq=str(ligneSamLue[9])     
    CIGAR=str(ligneSamLue[5])
    NM=(ligneSamLue[11].split(":"))[2]    
    if int(Flag)==0: 
        if str(CIGAR)=="145M":
            SNP=Seq[int(PosSnp)]
            if int(NM)> int(MaxMismatch):
                ReadState=0
            if int(NM)==int(MaxMismatch):
                NbreNMMax=1
        else:
            SNP=''
            ReadState=0 # Aligned reads flagged as secondary (i.e. with Flag=="246") are automatically considered as outliers.
    else:
        SNP=''
        ReadState=0        
    return idSeqCurrent,ReadState,SNP,NbreNMMax        

##############################################################################

def ajoutAllele(Dico,allele): 
    # Input : aligned reads against a candidate sequence
    # OutPut : a dictionary containing the number of reads per allelic state at the SNP candidate marker
    
    if str(allele)<>'':
        if str(allele) in Dico.keys():
            Dico[allele]=int(Dico[allele])+1
        else:
            Dico[allele]=1      
    return Dico       

##############################################################################

def RatioAllele(listeTuple): # unused in the default settings
    SommeEffectif=0
    for allele in SortTabSnp:
        SommeEffectif= SommeEffectif+allele[1]    

    ratio=(float(listeTuple[1][1])/float(SommeEffectif))
    return ratio,SommeEffectif
    
##############################################################################    

def AjoutLineFasta(Fichier, ListSamLine): # Writing the partially SE pre-validated fasta file.
    # Fichier is linked to a writable file.
    for SamLine in ListSamLine:       
        CIGAR=str(SamLine[5])
        NM=(SamLine[11].split(":"))[2]
        if str(CIGAR)=="145M":
            if int(NM)==0:
                idSeqRef=str(SamLine[2])     
                Seq=str(SamLine[9])
                Fichier.write("%s%s\n" % (">",idSeqRef))
                Fichier.write("%s\n" %	(Seq))
                break
        
    
    
    
###############################################################################        

def AjoutLineSam(FichierOut, ListSamLine): # Writing a new .sam file containing only candidate sequences (and associated aligned reads) that passed the second filtering (1).
    for SamLine in ListSamLine:
        NewLine=""
        for elts in SamLine:
            NewLine=NewLine+str(elts)+"\t"
        FichierOut.write("%s\n" % (NewLine))
    
    
###############################################################################
# Intialisation of variables        
###############################################################################       

# Variables for alignments on a candidate sequence
TabSnp={}         # Dictionary of alleles of a candidate SNP
Alignement=list() # Buffer table of the alignment correspong to idSeqCurrent (i.e, the current candidate sequence)
NbreBadRead=0     # Counting the number of outlier reads in an alignment
keep=1            # The alignment/current candidate sequence is kept by default  
MaxNM=0           # Counting reads with >MaxMismatch in an alignment

# Variables for the whole sub-program
SeqTotal=0        # Counting number of total candidate sequences 
SeqGood=0         # Counting number of SE pre-validated candidate sequences
SeqWrong=0        # Counting number of discarded candidate sequences


###############################################################################


print "Beginning to read the .sam file..."        

# Initialisation:

idSeqRef='' # Current read identifiant.
        
for refligne in ReadCsvFile:
    idSeqCurrent=str(refligne[2]) 
    # At the beginning of the .sam file screening, idSeqRef is empty.
    if idSeqRef=='': 
       idSeqRef=idSeqCurrent
       ResultatsRead=testRead(refligne,MaxMismatch)
       # Counting outliers read in an alignment.
       if int(ResultatsRead[1])==0:NbreBadRead+=1
       if int(ResultatsRead[3])==1:MaxNM+=1     
       # Writing a new .sam file with retained alignments.       
       Alignement.append(refligne)    
       # Filling table of alleles                  
       TabSnp=ajoutAllele(TabSnp,ResultatsRead[2])
        
      
       
       
       
   # Pursuing the .sam file screening...   
    else:
        if idSeqRef==idSeqCurrent: # While being in the same alignment (i.e., still screening reads of a given candidate sequence)
         
           ResultatsRead=testRead(refligne,MaxMismatch) # Checking the state of the read (outlier or not)
           # Incrementing the counter of outlier reads inside the alignment
           if int(ResultatsRead[1])==0:NbreBadRead+=1
           if int(ResultatsRead[3])==1:MaxNM+=1                
           # on ajoute la ligne SAM dans le futur fichier SAM Bon ou Mauvais        
           Alignement.append(refligne)    
           #on rempli le tableau d alleles
           TabSnp=ajoutAllele(TabSnp,ResultatsRead[2])   
                 
           #print ResultatsRead 
           
         
        else: # idSeqRef is not emplty but different from idSeqCurrent, which means we switched to another candidate sequence alignment.
        # The former candidate sequence alignment is then analyzed as a whole:        
            
            # Sorting the allele table as a tuple list        
            SortTabSnp=sorted(TabSnp.iteritems(), key=lambda item: -item[1])    
            # Discard of the candidate if too many outlier reads:
            if int(NbreBadRead)>int(nbRead):keep=0             
            # Discard of the candidate if supplementary alleles (e.g. tri- quadriallelic marker) is present with depth coverage >1.
            if len(SortTabSnp)==4:keep=0
            if len(SortTabSnp)==3: 
                if int(SortTabSnp[int(len(SortTabSnp))-1][1])>1:
                    keep=0
                else:
                    CalcRatio=RatioAllele(SortTabSnp)
            
            # SNP information:     
            if len(SortTabSnp)==2: 
                CalcRatio=RatioAllele(SortTabSnp) # Obtaining the ratio between depth coverage of each two alleles for information.
                
            # It is possible to discard a candidate if the ratio between depth coverage of the two main allele exceed a certain value. This option has not been retained in Delord et al. due to low depth coverage that may be unrepresentative of actual allele frequencies in the discovery pool.
            # With a higher depth coverage (e.g. >50) it is possible to discard an alignment where the ratio approaches 0.50 as it might result from paralogy.
            # as follows:
            # if float(CalcRatio[0])>float(MaxRatio):keep=0
         
            # Checking the retained candidates sequences:          
            if keep==1: 
                SeqGood+=1
                print idSeqRef 
                print SortTabSnp
                print "Ratio = "+str(CalcRatio[0])
                print "Depth coverage of each allele = "+str(CalcRatio[1])
                print "Size of the new .sam file ="+str(len(Alignement))
                print ""    
                
                # Generating the SE pre-validated fasta file:                
                AjoutLineFasta(FastaFileOut,Alignement)

##################################################################################           
           
# Re-Initialisation:
           
            del Alignement[:]
            TabSnp={}
            NbreBadRead=0
            keep=1
            MaxNM=0
            # on affecte les nouvelles valeurs
            idSeqRef=idSeqCurrent              
            ResultatsRead=testRead(refligne,MaxMismatch)
            if int(ResultatsRead[1])==0:NbreBadRead+=1            
            if int(ResultatsRead[3])==1:MaxNM+=1                
            Alignement.append(refligne)             
            TabSnp=ajoutAllele(TabSnp,ResultatsRead[2])
           
            
print "Number of candidate sequences that passed Second filtering (1) = "+str(SeqGood)    
    
 
#########################################################################                
        # .sam file screening is complete: analyzing the last candidate.    
########################################################################               

if int(NbreBadRead)>int(nbRead):keep=0

# Sorting the allele table as a tuple list       
SortTabSnp=sorted(TabSnp.iteritems(), key=lambda item: -item[1])    

if len(SortTabSnp)==4:keep=0
if len(SortTabSnp)==3: 
    if int(SortTabSnp[int(len(SortTabSnp))-1][1])>1:
        keep=0
    else:
        CalcRatio=RatioAllele(SortTabSnp)

# SNP information:   
if len(SortTabSnp)==2: 
    CalcRatio=RatioAllele(SortTabSnp)     
          
# if float(CalcRatio[0])>float(MaxRatio):keep=0 
if keep==1: 
    SeqGood+=1
    print idSeqRef 
    print SortTabSnp
    print "ratio = "+str(CalcRatio[0])
    print "Depth coverage of each allele = "+str(CalcRatio[1])
    print "Size of the new .sam file ="+str(len(Alignement))
    print ""    
              
    # Completing the SE pre-validated fasta file: 
    AjoutLineFasta(FastaFileOut,Alignement)
 
FastaFileOut.close
 
## End.   
