input=$1
threads=$2
sample=$3
cluster_dir=$4
output=$5
folder=$6

qual=$(cat ${folder}07_stats/Temp/quality_${sample}.temp)

if [ $qual = '17' ];then
	nb_id=0.95
elif [ $qual = '16' ];then
	nb_id=0.94
elif [ $qual = '15' ];then
	nb_id=0.93
elif [ $qual = '14' ];then
	nb_id=0.92
elif [ $qual = '13' ];then
	nb_id=0.89
elif [ $qual = '12' ];then
	nb_id=0.87
elif [ $qual = '11' ];then 
	nb_id=0.83
elif [ $qual = '10' ];then
	nb_id=0.79
elif [ $qual = '9' ];then
	nb_id=0.74
elif [ $qual = '8' ];then
	nb_id=0.68
elif [ $qual = '7' ];then
	nb_id=0.59
elif [ $qual = '6' ];then
	nb_id=0.49
fi

vsearch --cluster_size ${input} --id ${nb_id} --iddef 0 --strand both --clusterout_sort --threads ${threads} --consout ${output} --clusters ${cluster_dir} --fasta_width 0
