input=$1
threads=$2
sample=$3
cluster_dir=$4
cov=$5
output=$6
folder=$7

qual=$(cat ${folder}06_stats/Temp/quality_${sample}.temp)

if (( qual == 17 )); then
	nb_id=0.96
elif (( qual==16 ));then
	nb_id=0.95
elif (( qual==15 )); then
	nb_id=0.937
elif (( qual==14 )); then
	nb_id=0.92
elif (( qual==13 ));  then
	nb_id=0.90
elif (( qual==12 )); then
	nb_id=0.875
elif (( qual==11 )); then 
	nb_id=0.841
elif (( qual==10 )); then
	nb_id=0.80
elif (( qual==9 )); then
	nb_id=0.748
elif (( qual==8 )); then
	nb_id=0.683
elif (( qual==7)); then
	nb_id=0.60
elif ((qual<=6)); then
	nb_id=0.498
fi

vsearch --cluster_size ${input} --id ${nb_id} --strand both --clusterout_sort --threads ${threads} --consout ${folder}03_vsearch/consensus_${sample}.fasta --clusters ${cluster_dir} --fasta_width 0

echo $nb_id > ${folder}06_stats/Temp/perc_cons_${sample}.temp

mv ${folder}03_vsearch/consensus_${sample}.fasta ${output}
