input=$1
output_SNP=$2
output_INDEL=$3
output_CNS=$4
sample=$5
folder=$6

qual=$(cat ${folder}07_stats/Temp/quality_${sample}.temp)

varfreq=$(python workflow/scripts/py_scripts/math_calcul2.py $qual $sample $folder) 


varscan pileup2snp $input --p-value 0.01 --min-var-freq $varfreq  > $output_SNP  
varscan pileup2indel $input --p-value 0.01 --min-var-freq $varfreq > $output_INDEL
varscan pileup2cns $input --p-value 0.01 --min-var-freq $varfreq > $output_CNS
