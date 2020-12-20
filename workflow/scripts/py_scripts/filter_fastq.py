from Bio import SeqIO
import os
import math
import gzip
import sys

name = sys.argv[1] 
min =float(sys.argv[2])
max =float(sys.argv[3])
qs = float(sys.argv[4])


with gzip.open(sys.argv[5], 'wt') as output_handle:
    with gzip.open(name, "rt") as handle:
        for rec in SeqIO.parse(handle, "fastq") :
            rec.letter_annotations["phred_quality"]
            probs = []
            for q in rec.letter_annotations["phred_quality"]:
                e = float(math.pow(10.0,-1*(float(q)/10.0)))
                probs.append(e)
            av_prob = float(sum(probs))/float(len(rec.letter_annotations["phred_quality"]))
            av_q = float(-10.0*(math.log10(float(av_prob))))
		
            if len(rec.seq) >=min and len(rec.seq) <=max and av_q >= qs:
		# Add this record to our list
                SeqIO.write(rec, output_handle, "fastq")
			
