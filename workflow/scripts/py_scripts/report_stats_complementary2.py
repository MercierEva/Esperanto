import pandas as pd
import sys
import subprocess
from collections import OrderedDict
from Bio import SeqIO

class Report_Stat:

    def __init__(self):
        self.d = {'sample':[],"IdentityPercent_2":[], "Depth_2":[], 'mapping depth':[], 'covered length':[]}

    def complete_array(self, sample, nb_id, count, map_depth, cov):
        self.d['sample'].append(sample)
        self.d["IdentityPercent_2"].append(nb_id)
        self.d["Depth_2"].append(count)
        self.d['mapping depth'].append(map_depth)
        self.d['covered length'].append(cov)
        return OrderedDict(self.d)
    
    def read_nb_id(self, folder, sample):
        with open(folder + "06_stats/Temp/perc_cons_"+sample+".temp", "r") as f : 
            nb_id = f.read()
        return round(float(nb_id), 2)

    def count_seq(self, sample, folder):
        path = folder + "03_vsearch/PASS/consensus_" + sample +".fasta" 
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count=line.split('seqs=')[1:][0]
                    break
        return count

    def write_array(self, d1, d2, output):
        df1 = pd.read_csv(d1, sep='\t') 
        df2 = pd.DataFrame(data=d2)  
        df = pd.merge(df1, df2, on='sample', how='right')    
        df.to_csv(output, sep='\t', encoding='utf-8', index=False)

    def length_seq(self, input_file):
        total=0
        handle = open(input_file, 'rU')
        SeqRecords = SeqIO.parse(handle, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)
            total += len(Seq)
        return total
    
    def calcul_depth(self, bamfile):
        cmd = "samtools mpileup " + bamfile + " | awk -v X=5 '$4>=X' | wc -l"
        depth = subprocess.check_output(cmd, shell=True)
        return float(depth)       
    
    def calcul_coverage(self, depth, length):
        coverage = depth/length*100
        return round(float(coverage), 2)
        
if __name__ == "__main__":
    new_report = Report_Stat()
    count = new_report.count_seq(sys.argv[1], sys.argv[7])
    nb_id = new_report.read_nb_id(sys.argv[7], sys.argv[1])
    length = new_report.length_seq(sys.argv[2] + '/' + sys.argv[4])
    map_depth = new_report.calcul_depth(sys.argv[2] + '/' + sys.argv[3])
    cov = new_report.calcul_coverage(map_depth, length)
    array_final = new_report.complete_array(sys.argv[1], nb_id, count, map_depth, cov)
    new_report.write_array(sys.argv[5], array_final, sys.argv[6])
