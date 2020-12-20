cluster_dir=$1
input=$2
threads=$3
model=$4
output=$5
sample=$6
report_out=$7
folder=$8


while read line ; do
  if [ ${line:0:1} == ">" ] ; then
    ((k++))
    outfile=$(echo ${line:1:${#line}})
    echo $line > "${folder}"02_vsearch/"${k}"_"${outfile}".fa
  else
    echo $line >> "${folder}"02_vsearch/"${k}"_"${outfile}".fa
  fi
done < $input

len=$(( $(ls -1v ${folder}02_vsearch/*${sample}*.fa | wc -l) ))
for ((k=1; k <= $len ; k++))
do
    ech=$(ls -1v ${folder}02_vsearch/*${sample}*.fa | sed -n ${k}p)
    echo $ech
    d=$(ls -S ${cluster_dir}* | egrep 'cluster_([0-9]{1,5})$' | head -$len | sed -n ${k}p)
    echo $d
    medaka_consensus -i $ech -d $d -t $threads -o "${folder}"03_medaka/"${sample}"_medaka_cluster_"${k}" -m $model
    python scripts/py_scripts/report_stats_complementary2_otu.py ${sample} $d "${folder}"03_medaka/"${sample}"_medaka_cluster_"${k}" calls_to_draft.bam consensus.fasta ${report_out} ${k} ${len} ${folder}
    cat "${folder}"03_medaka/"${sample}"_medaka_cluster_"${k}"/consensus.fasta | awk '/^>/{split($0,a," "); print a[1]"_medaka_'${k}'" ; next}{print $0}' >> ${output} 
done
