import pandas as pd
import sys
import subprocess
from collections import OrderedDict
import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

class Report_Stat:

    def __init__(self):
        self.d = {'sample':[], "QThreshold_1":[], "Depth_1":[], "IdentityPercent_2":[], "Depth_2":[], 'mapping depth':[], 'covered length':[]}

    def complete_array(self, sample, qual, count_1, nb_id, count_2, map_depth, cov):
        self.d['sample'].append(sample)
        self.d["QThreshold_1"].append(qual)
        self.d["Depth_1"].append(count_1)
        self.d["IdentityPercent_2"].append(nb_id)
        self.d["Depth_2"].append(count_2)
        self.d['mapping depth'].append(map_depth)
        self.d['covered length'].append(cov)
        return OrderedDict(self.d)
        

    def return_qual(self, sample, folder):
        with open(folder +"/06_stats/Temp/quality_" + sample +".temp", "r") as qual:
            q = qual.read()
        return q

    def read_nb_id(self, folder, sample):
        with open(folder + "06_stats/Temp/perc_cons_"+sample+".temp", "r") as f : 
            nb_id = f.read()
        return round(float(nb_id), 2)
    
    def count_seq_fastq(self, input_file):
        count=0
        with gzip.open(input_file, 'rt') as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                count += 1
        return count

    def count_seq(self, input_file):
        count = 0 
        with open(input_file) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                count += 1
        return count
    
    def length_seq(self, input_file):
        length = 0 
        with open(input_file) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                length = len(seq)
        return length

    def write_array(self, d, output):
        df = pd.DataFrame(data=d)  
        df.to_csv(output, sep='\t', encoding='utf-8', index=False)
    
    def calcul_l_depth(self, bamfile):
        cmd = "samtools mpileup " + bamfile + " | awk -v X=1 '$4>=X' | wc -l"
        depth = subprocess.check_output(cmd, shell=True)
        return float(depth)

    def calcul_depth(self, bamfile):
        try :
            cmd = "samtools depth -a "+ bamfile +" | awk '{sum+=$3} END { print sum/NR}'"
            d = subprocess.check_output(cmd, shell=True)
            return float(d)  
        except ZeroDivisionError:
            d = 0
            return float(d)
  
    def calcul_coverage(self, depth, length):
        try:
            coverage = depth/length*100
            return round(float(coverage), 2)
        except ZeroDivisionError:
            coverage = "NA"
    
    def append_on_init(self, folder, sample, output):
        path=folder + "06_stats/Temp/" + sample + "*.tsv" 
        list_of_array_by_sample = sorted(glob.glob(path))
        list_of_dataframes = []
        for filename in list_of_array_by_sample:
            list_of_dataframes.append(pd.read_csv(filename, sep='\t'))
        merged_df = pd.concat(list_of_dataframes, ignore_index=True)
        merged_df.to_csv(output, sep='\t', index=False, encoding='utf-8')
        
if __name__ == "__main__":
    new_report = Report_Stat()
    quality =new_report.return_qual(sys.argv[1], sys.argv[9])
    depth=new_report.count_seq_fastq(sys.argv[9]+"01_nanofilt/PASS/"+sys.argv[1]+"_filt.fastq.gz")
    count = new_report.count_seq(sys.argv[2])
    nb_id = new_report.read_nb_id(sys.argv[9], sys.argv[1])
    length = new_report.length_seq(sys.argv[3] + '/' + sys.argv[5])
    map_depth = new_report.calcul_l_depth(sys.argv[3] + '/' + sys.argv[4])
    cov = new_report.calcul_coverage(map_depth, length)
    d = new_report.calcul_depth(sys.argv[3] + '/' + sys.argv[4])
    array_final = new_report.complete_array(sys.argv[1], quality, depth, nb_id, count, d, cov)
    path_cluster =  sys.argv[9] +"06_stats/Temp/" + sys.argv[1] + "_cluster_" + sys.argv[7]+ ".tsv"
    array_prog = new_report.write_array(array_final, path_cluster)
    if sys.argv[7]==sys.argv[8]:
        print("condition_finale")
        df = new_report.append_on_init(sys.argv[9], sys.argv[1], sys.argv[6])
    
        
