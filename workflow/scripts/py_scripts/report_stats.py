from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import sys
import gzip
import yaml
import subprocess
import math

class Report_Stat:

    def __init__(self):
        self.d = {'Sample':[], "Quality_Consensus_Final": [], "QThreshold_1":[], "Depth_1":[], "IdentityPercent_2":[], "Depth_2":[], "Mapping depth":[], "Covered length":[],  "Length":[], "Sequence fasta":[]}

    def complete_array(self, sample, qual_fin, qual, countfastq, nb_id, countcluster, map_depth, cov, length, fasta):
        self.d["Sample"].append(sample)
        self.d["Quality_Consensus_Final"].append(round(qual_fin))
        self.d["QThreshold_1"].append(int(qual))
        self.d["Depth_1"].append(countfastq)
        self.d["IdentityPercent_2"].append(nb_id)
        self.d["Depth_2"].append(int(countcluster))
        self.d["Mapping depth"].append(map_depth)
        self.d["Covered length"].append(cov)
        self.d["Length"].append(length)
        self.d["Sequence fasta"].append(fasta)
        return self.d

    def return_qual(self, sample, folder):
        with open(folder +"07_stats/Temp/quality_" + sample +".temp", "r") as qual:
            q = qual.read()
        return q


    def count_seq_fasta(self, sample, folder):
        path = folder + "01_nanofilt/PASS/" + sample +"_filt.fasta" 
        count=0
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count+=1
        return count

    def read_nb_id(self, quality):
        quality=int(quality)
        if quality == 17 :
	        nb_id=0.95
        elif quality == 16 :
	        nb_id=0.94
        elif quality == 15 :
	        nb_id=0.93
        elif quality == 14 :
	        nb_id=0.92
        elif quality == 13 :
	        nb_id=0.89
        elif quality == 12 :
            nb_id=0.87
        elif quality == 11 : 
            nb_id=0.83
        elif quality == 10 :
            nb_id=0.79
        elif quality == 9 :
            nb_id=0.74
        elif quality == 8 :
            nb_id=0.68
        elif quality == 7 :
            nb_id=0.59
        elif quality == 6 :
            nb_id=0.49
        
        return nb_id
    
    def count_seq_cluster(self, sample, folder):
        path = folder + "02_vsearch/PASS/consensus_" + sample +".fasta" 
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count=line.split('seqs=')[1:][0]
                    break
        return count
        
    def calcul_depth(self, sample, folder):
        cmd = "samtools mpileup " + folder +"04_multialignment/"+ sample +"_MSA.sorted.bam | awk -v X=5 '$4>=X' | wc -l"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        out, err = proc.communicate()
        depth = float(out.decode('utf-8'))
        return depth       
    
    def length_seq(self, sample, folder):
        total=0
        input_file=folder+"03_porechop/"+sample+"_ont.fasta"
        handle = open(input_file, 'rU')       
        SeqRecords = SeqIO.parse(handle, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)
            total += len(Seq)
        return total
    
    def calcul_coverage(self, depth, length):
        coverage = depth/length*100
        return round(float(coverage), 2)

    def write_array(self, d, output):
        df = pd.DataFrame(data=d, columns=['Sample', "Quality_Consensus_Final", "QThreshold_1", "Depth_1", "IdentityPercent_2", "Depth_2", "Mapping depth", "Covered length", "Length", "Sequence fasta"])
        df.to_csv(output, sep='\t', encoding='utf-8', index=False)

    def open_fasta(self, sample, folder):
        sequence = ''
        fastafile=folder+"03_porechop/"+sample+"_ont.fasta"
        with open(fastafile) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    continue
                sequence += line  
        return sequence   

    def calcul_quality_final(self, quality, depth):
        a=float()
        b=float()
        quality=int(quality)
        if quality == 17 :
	        a, b = 16.72, 11.285
        elif quality == 16 :
	        a, b = 15.736, 10.621
        elif quality == 15 :
	        a, b = 14.753, 9.957
        elif quality == 14 :
	        a, b = 13.769, 9.293
        elif quality == 13 :
	        a, b = 12.786, 8.63
        elif quality == 12 :
            a, b = 11.802, 7.966
        elif quality == 11 : 
            a, b = 10.819, 7.302
        elif quality == 10 :
            a, b = 9.835, 6.638
        elif quality == 9 :
            a, b = 8.8518, 5.9743
        elif quality == 8 :
            a, b = 7.8682, 5.3105
        elif quality == 7 :
            a, b = 6.8847, 4.6467
        else :
            a, b = 5.9012, 3.9829
        ln_depth = math.log(float(depth))
        calcul=a*ln_depth+b
        return calcul 
        
if __name__ == "__main__":
    new_report = Report_Stat()
    quality = new_report.return_qual(sys.argv[1], sys.argv[2])
    countdepth1 = new_report.count_seq_fasta(sys.argv[1], sys.argv[2])
    nb_id = new_report.read_nb_id(quality)
    countdepth2 = new_report.count_seq_cluster(sys.argv[1], sys.argv[2])
    map_depth = new_report.calcul_depth(sys.argv[1], sys.argv[2])
    length = new_report.length_seq(sys.argv[1], sys.argv[2])
    cov = new_report.calcul_coverage(map_depth, length)
    sequence=new_report.open_fasta(sys.argv[1], sys.argv[2])
    qual_fin=new_report.calcul_quality_final(quality, countdepth2)
    array = new_report.complete_array(sys.argv[1], qual_fin, quality, countdepth1, nb_id, countdepth2, map_depth, cov, length, sequence)
    new_report.write_array(array, sys.argv[3])
        

