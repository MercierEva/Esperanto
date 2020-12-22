#!bin/bash

input=$1
min_length=$2 
max_length=$3 
rd=$4
cov=$5
sample=$6
output=$7
folder=$8
file_interm="${folder}01_nanofilt/${sample}_filt.fastq.gz"

if [ ! -d ${folder}06_stats/ ]; then
    mkdir ${folder}06_stats
fi

if [ ! -d ${folder}06_stats/Temp/ ]; then
    mkdir ${folder}06_stats/Temp
fi


for (( q = 17; q > 8 ; q-- )) ; do
    gunzip -c $input | NanoFilt -l $min_length --maxlength $max_length -q $q --readtype $rd | gzip > ${file_interm}
    #python scripts/py_scripts/filter_fastq.py $input $min_length $max_length $q $file_interm
    if gzip -t $file_interm && [[ $(gunzip -c $file_interm | head -c1 | wc -c) != 0 ]] ; then
        count=$(( $( zcat $file_interm | wc -l) / 4 ))
        if [[ "$count" -gt "$cov"  ]] ;then
            mv $file_interm $output
            break 
	    fi
    fi
    if [[ "$q" -eq 9 ]]; then
        mv $file_interm $output
           break 
    fi	    
done 


echo "$q" >  ${folder}06_stats/Temp/quality_${sample}.temp

