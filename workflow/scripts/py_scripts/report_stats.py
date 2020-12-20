from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import sys
import gzip
from collections import OrderedDict
import yaml

class Report_Stat:

    def __init__(self):
        self.d = {'sample':[], "QThreshold_1":[], "Depth_1":[]}

    def complete_array(self, sample, qual, count):
        self.d['sample'].append(sample)
        self.d["QThreshold_1"].append(qual)
        self.d["Depth_1"].append(count)
        return self.d

    def write_array(self, d, output):
        df = pd.DataFrame(data=OrderedDict(d))
        df.to_csv(output, sep='\t', encoding='utf-8', index=False)

    def count_seq_fastq(self, input_file):
        count=0
        with gzip.open(input_file, 'rt') as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                count += 1
        return count

    def return_qual(self, sample, folder):
        with open(folder +"06_stats/Temp/quality_" + sample +".temp", "r") as qual:
            q = qual.read()
        return q
        
if __name__ == "__main__":
    new_report = Report_Stat()
    quality = new_report.return_qual(sys.argv[1], sys.argv[4])
    count = new_report.count_seq_fastq(sys.argv[3])
    array = new_report.complete_array(sys.argv[1], quality, count)
    new_report.write_array(array, sys.argv[2])
        

