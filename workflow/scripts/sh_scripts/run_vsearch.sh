input=$1
threads=$2
sample=$3
cluster_dir=$4
cov=$(echo "$5 * 0.75" | bc )
output=$6
folder=$7

qual=$(cat ${folder}06_stats/Temp/quality_${sample}.temp)


if (( qual == 17 ));then
	nb_id=0.96
elif (( qual==16 ));then
	nb_id=0.95
elif (( qual==15 ));then
	nb_id=0.937
elif (( qual==14 ));then
	nb_id=0.92
elif (( qual==13 ));then
	nb_id=0.90
elif (( qual==12 ));then
	nb_id=0.875
elif (( qual==11 ));then 
	nb_id=0.841
elif (( qual==10 ));then
	nb_id=0.80
elif (( qual==9 ));then
	nb_id=0.748
elif (( qual==8 ));then
	nb_id=0.683
elif (( qual==7 ));then
	nb_id=0.60
elif (( qual<=6 ));then
	nb_id=0.498
fi

vsearch --cluster_size ${input} --id ${nb_id} --iddef 2 --strand both --clusterout_sort --threads ${threads} --consout ${folder}03_vsearch/consensus_${sample}.fasta --clusters ${cluster_dir} --fasta_width 0
cluster_bigger=$(ls -S ${cluster_dir}* | head -1)
count=$(grep -c '>' $cluster_bigger )
while [[ -e $cluster_bigger ]] && (( $(echo "$count < $cov" | bc -l) )) ; do 
    if (( $(echo "$nb_id < 0.5" | bc -l) )) ; then
        mv ${folder}03_vsearch/consensus_${sample}.fasta ${output}
    else
	nb_id=$( echo "$nb_id - 0.02" | bc -l )
	vsearch --cluster_size ${input} --id ${nb_id} --iddef 2 --strand both --clusterout_sort --threads ${threads} --consout ${folder}03_vsearch/consensus_${sample}.fasta --clusters ${cluster_dir} --fasta_width 0  
	cluster_bigger=$(ls -S ${cluster_dir}* | head -1)
  	count=$(grep -c '>' $cluster_bigger )
    fi	  
done
echo $nb_id >  ${folder}06_stats/Temp/perc_cons_${sample}.temp

if (( $(echo "$count >= $cov" | bc -l) )) ; then
    mv ${folder}03_vsearch/consensus_${sample}.fasta ${output}
fi
