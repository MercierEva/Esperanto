input=$1
min_length=$2 
max_length=$3 
rd=$4
quality_cons=$5
sample=$6
output=$7
folder=$8
file_interm="${folder}01_nanofilt/${sample}_filt.fastq.gz"

if [ ! -d ${folder}07_stats/ ]; then
    mkdir ${folder}07_stats
fi

if [ ! -d ${folder}07_stats/Temp/ ]; then
    mkdir ${folder}07_stats/Temp
fi


for (( q = 17; q > 5 ; q-- )) ; do
    gunzip -c $input | NanoFilt -l $min_length --maxlength $max_length -q $q --readtype $rd | gzip > ${file_interm}  
    cov=$(python scripts/py_scripts/math_calcul.py ${q} ${quality_cons})
    if gzip -t $file_interm && [[ $(gunzip -c $file_interm | head -c1 | wc -c) != 0 ]] ; then
        count=$(( $( zcat $file_interm | wc -l) / 4 ))
        if [[ "$count" -gt "$cov"  ]] ;then
            mv $file_interm $output
            break 
	    fi
    fi
    if [[ "$q" -eq 6 ]]; then
        mv $file_interm $output
        break 
    fi	
done 

echo "$q" >  ${folder}07_stats/Temp/quality_${sample}.temp

