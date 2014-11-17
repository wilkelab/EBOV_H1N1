# SJS. quick script to check out aa/nuc alignment variation in 2014 outbreaks.
from Bio import AlignIO
import numpy as np
import sys

codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT", "NNN"]
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codon_dict   = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F", "NNN":"X"}
nuc = ['A', 'C', 'G', 'T']

for gene in ['np', 'l', 'gp', 'vp24_', 'vp30_', 'vp35_', 'vp40_']:
    print "\n\n========== " + gene + " =========="
    aln = AlignIO.read("EBOV_2014/" + gene + "2014_unique.aln", "fasta")
    ids = []
    numseq = len(aln)
    print numseq
    continue
    alnlen = len(aln[0])
    aamat = np.empty( [numseq, alnlen/3], dtype = 'int8' )
    nucmat = np.empty( [numseq, alnlen], dtype = 'int8')
    i=0
    for entry in aln:
        id = str(entry.id)
        seq = str(entry.seq)
        for j in range(0,alnlen,3):
            try:
                amino = codon_dict[ seq[j:j+3] ]
                aamat[i][j/3] = amino_acids.index(amino)
            except:
                aamat[i][j/3] = 100
        for j in range(alnlen):
            try:
                nucmat[i][j] = nuc.index(seq[j])
            except:
                nucmat[i][j] = 4
        ids.append( str(entry.id) )      
        i += 1
    
    
    informative = []
    for col in aamat.T:
        colset = list(set(col))
        if len(colset) == 1 or len(colset)==2 and 100 in colset:
            continue
        else:
            informative.append(col)


    print "amino informative:"
    x = 0
    for col in np.array(informative).T:
        colstr = ''
        for row in col:
            try:
                colstr += amino_acids[row]
            except:
                colstr += 'X'
        print ids[x], colstr
        x += 1

  
    informative = []
    for col in nucmat.T:
        colset = list(set(col))
        if len(colset) == 1 or len(colset)==2 and 4 in colset:
            continue
        else:
            informative.append(col)  

    print "\n\nnuc informative:"
    x=0
    for col in np.array(informative).T:
        colstr = ''
        for row in col:
            try:
                colstr += nuc[row]
            except:
                colstr += 'N'
        print ids[x], colstr
        x += 1
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


