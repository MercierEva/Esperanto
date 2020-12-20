#-*- coding: utf-8 -*-​
import sys
from Bio import SeqIO
import shutil

#create dictionary with IUPAC codes
IUPAC = {}
IUPAC['A'] = "A"
IUPAC['C'] = "C"
IUPAC['G'] = "G"
IUPAC['T'] = "T"
IUPAC['M'] = "AC"
IUPAC['R'] = "AG"
IUPAC['W'] = "AT"
IUPAC['S'] = "CG"
IUPAC['Y'] = "CT"
IUPAC['K'] = "GT"
IUPAC['V'] = "ACG"
IUPAC['H'] = "ACT"
IUPAC['D'] = "AGT"
IUPAC['B'] = "CGT"
IUPAC['X'] = "GATC"
IUPAC['N'] = "GATC"

def MatchLetter(a, b):
    global IUPAC
    try:
        sa = IUPAC[a.upper()]
    except:
        return False
    try:
        sb = IUPAC[b.upper()]
    except:
        return False
    for ca in sa:
        if ca in sb:
            return True
    return False

def MatchPrefix(Seq, Primer):
    L = len(Seq)
    n = len(Primer)
    if L < n:
        n = L
    Diffs = 0
    for i in range(0, n):
        if not MatchLetter(Seq[i], Primer[i]):
            Diffs += 1
    return Diffs


with open(sys.argv[1]) as f:  
    for line in f:
        line = line.strip()
        if not line:           
            continue
        if line.startswith(">"):
            nb_seq_centroid = line.split('with_')[1:][0]
            if int(nb_seq_centroid) > 20 :
                with open(sys.argv[2]+"03_vsearch/PASS/fasta_filter_" + sys.argv[3] + ".fasta", 'a') as f :      
                    f.write(line+'\n')
        else:
            if int(nb_seq_centroid) > 20 :  
                with open(sys.argv[2]+"03_vsearch/PASS/fasta_filter_" + sys.argv[3] + ".fasta", 'a') as f :      
                    f.write(line+'\n')


if len(sys.argv) == 4 :

    primer = sys.argv[5]
    handle = open(sys.argv[2] + "03_vsearch/PASS/fasta_filter_" + sys.argv[3] + ".fasta", 'rU')
    SeqRecords = SeqIO.parse(handle, 'fasta')
    SeqRecords_Forward = []
    for rec in SeqRecords:
        Seq = str(rec.seq)
        Diffs = MatchPrefix(Seq, primer) #count diffs between Seq and primer
        if Diffs < 2 :  #if Diffs less than some threshold, then just print record
            SeqRecords_Forward.append(rec.reverse_complement())
        else : 
            SeqRecords_Forward.append(rec)
    SeqIO.write(SeqRecords_Forward, sys.argv[4], 'fasta')

else: 
    try :
        shutil.move(sys.argv[2] + "03_vsearch/PASS/fasta_filter_" + sys.argv[3] + ".fasta", sys.argv[4])
    except FileNotFoundError :
        with open(sys.argv[4], mode='a'): pass


