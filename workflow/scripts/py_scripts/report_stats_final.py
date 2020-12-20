import pandas as pd
import sys
from collections import OrderedDict
from Bio import SeqIO

class Report_Stat:

    def __init__(self):
        self.d = {'sample':[], 'length':[], 'sequence fasta':[]}

    def complete_array(self, sample, length, fasta):
        self.d['sample'].append(sample)
        self.d['length'].append(length)
        self.d['sequence fasta'].append(fasta)
        return self.d
    
    def open_fasta(self, fastafile):
        sequence = ''
        with open(fastafile) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    continue
                sequence += line  
        return sequence

    def write_array(self, d1, d2, output):
        df1 = pd.read_csv(d1, sep='\t') 
        df2 = pd.DataFrame(data=OrderedDict(d2))  
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
        
if __name__ == "__main__":
    new_report = Report_Stat()
    fasta=new_report.open_fasta(sys.argv[2])
    length=new_report.length_seq(sys.argv[2])
    array_final = new_report.complete_array(sys.argv[1], length, fasta)
    new_report.write_array(sys.argv[3], array_final, sys.argv[4])
        
