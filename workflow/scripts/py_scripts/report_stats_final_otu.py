import pandas as pd
import sys
from collections import OrderedDict
from Bio import SeqIO
import os
sys.path.append(os.getcwd() + "/b3fbe6da678d3841a976")
from b3fbe6da678d3841a976 import FastaMLtoSL

class Report_Stat:

    def __init__(self, input_fastas_files):
        self.seq = []
        self.lengths = []
        self.inFile = input_fastas_files
        self.SL_fastas = FastaMLtoSL.FastaMLtoSL(self.inFile)
        self.SL_fastas.returnSLFasta()

    def separate_files(self):    
        handle = open(self.inFile+".out", 'rU')
        SeqRecords = SeqIO.parse(handle, 'fasta')
        number=1
        for rec in SeqRecords:
            Seq = str(rec.seq)
            while str(rec.id).endswith(str(number)) == False :
                self.seq.append("")
                number+=1
            self.seq.append(Seq)
            number+=1
        return self.seq
   
    def write_array(self, d1, d2, d3, output):
        df1 = pd.read_csv(d1, sep='\t') 
        df1_sorted = df1.sort_values('Depth_2',ascending=False)
        ds2 = pd.Series(d2, name='lengths')  
        ds3 = pd.Series(d3, name='sequencies')
        dfbis = pd.concat([df1, ds2], axis=1)  
        df = pd.concat([dfbis, ds3], axis=1)   
        df.to_csv(output, sep='\t', encoding='utf-8', index=False)
    
    def length_seq(self):
        handle = open(self.inFile + ".out", 'rU')
        SeqRecords = SeqIO.parse(handle, 'fasta')
        number=1
        for rec in SeqRecords:
            Seq = str(rec.seq)
            while str(rec.id).endswith(str(number)) == False :
                self.lengths.append("")
                number+=1
            self.lengths.append(len(Seq))
            number+=1
        return self.lengths
   
        
if __name__ == "__main__":
    new_report = Report_Stat(sys.argv[1])
    fastas=new_report.separate_files()
    lengths=new_report.length_seq()
    new_report.write_array(sys.argv[2], lengths, fastas, sys.argv[3])
        
