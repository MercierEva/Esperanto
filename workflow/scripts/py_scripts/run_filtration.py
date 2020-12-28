import subprocess
import sys
import gzip
import os
from math_calcul import Mathematica_1
import shutil

input_file=sys.argv[1]
min_length=sys.argv[2] 
max_length=sys.argv[3]
rd=sys.argv[4]
quality_cons=sys.argv[5]
sample=sys.argv[6]
output=sys.argv[7]
folder=sys.argv[8]
file_interm=str(folder) + "01_nanofilt/" + str(sample) + "_filt.fastq.gz"


try:
    os.mkdir(folder+"07_stats")
    os.mkdir(folder+"07_stats/Temp/")
except FileExistsError:
    pass

for q in range(17, 6, -1) : 
    print(q)
    unzipping=subprocess.run(["gunzip", "-c", str(input_file)] , check=True, capture_output=True)
    filtering = subprocess.run(["NanoFilt", "-l", str(min_length),"--maxlength",str(max_length),"-q",str(q),"--readtype", str(rd)], input=unzipping.stdout, capture_output=True)
    subprocess.run(["gzip"], input=filtering.stdout, stdout=gzip.open(file_interm, 'wb'))
    init_cov=Mathematica_1(q)
    cov=init_cov.calcul_depth_min(quality_cons)
    try:
        with gzip.open(file_interm, 'rb') as f:
            file_content = f.read(1)
            if len(file_content) > 0 :
                reads = []
                for line in f:
                    if line.decode('utf-8').startswith("@"):
                        reads.append(line)
                totalReads = len(reads)
            if totalReads > cov :
                break 
    except:
        pass

shutil.move(file_interm, output)

with open(folder+"07_stats/Temp/quality_"+sample+".temp", "w") as quality_file :
    quality_file.write(str(q))

