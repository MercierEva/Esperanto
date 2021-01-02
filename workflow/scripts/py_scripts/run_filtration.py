import subprocess
import sys
import gzip
import os
import shutil

input_file=sys.argv[1]
min_length=sys.argv[2] 
max_length=sys.argv[3]
rd=sys.argv[4]
sample=sys.argv[5]
output=sys.argv[6]
folder=sys.argv[7]
file_interm=str(folder) + "01_nanofilt/" + str(sample) + "_filt.fastq.gz"

def count_seq(input_file):
    proc = subprocess.Popen("zcat " + input_file + " | wc -l ", shell=True, stdout=subprocess.PIPE) 
    out, err = proc.communicate()
    count= int(out.decode('utf-8'))//4
    return count 

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
    if count_seq(file_interm) > 75 : 
        break
    else:
        continue

shutil.move(file_interm, output)

with open(folder+"07_stats/Temp/quality_"+sample+".temp", "w") as quality_file :
    quality_file.write(str(q))

