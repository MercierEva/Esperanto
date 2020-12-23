import sys
import math

def count_seq_cluster(sample, folder):
    path = folder + "02_vsearch/PASS/consensus_" + sample +".fasta" 
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count=line.split('seqs=')[1:][0]
                break
    return float(count)

def calcul_quality_final(quality, depth):     
    a=float()
    b=float()
    if quality == 17 :
        a=16.72
        b=11.285
    elif quality == 16 :
        a=15.736
        b=10.621
    elif quality == 15 :
        a=14.753
        b=9.957
    elif quality == 14 :
        a=13.769
        b=9.293
    elif quality == 13 :
        a=12.786
        b=8.63
    elif quality == 12 :
        a=11.802
        b=7.966
    elif quality == 11 : 
        a=10.819
        b=7.302
    elif quality == 10 :
        a=9.835
        b=6.638
    elif quality == 9 :
        a=8.8518
        b=5.9743
    elif quality == 8 :
        a=7.8682
        b=5.3105
    elif quality == 7 :
        a=6.8847
        b=4.6467
    else :
        a=5.9012
        b=3.9829
    
    ln_depth = math.log(depth)
    calcul=a*(ln_depth)+b
    return calcul 


def calcul(qual):
    prob=10**(-int(qual)/10)
    print(prob)
    return prob

depth=count_seq_cluster(sys.argv[2], sys.argv[3])
qual_final=calcul_quality_final(sys.argv[1], depth)
calcul(qual_final)