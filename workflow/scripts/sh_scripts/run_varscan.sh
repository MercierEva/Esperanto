input=$1
output_SNP=$2
output_INDEL=$3
output_CNS=$4
sample=$5
folder=$6

qual=$(cat ${folder}06_stats/Temp/quality_${sample}.temp)

if (( qual == 17 )); then
	varfreq=0.02
elif (( qual==16 )); then
	varfreq=0.025
elif (( qual==15 )); then
	varfreq=0.03
elif (( qual==14 )); then
	varfreq=0.04
elif (( qual==13 )); then
	varfreq=0.05
elif (( qual==12 )); then
	varfreq=0.06
elif (( qual==11 )); then 
	varfreq=0.08
elif (( qual==10 )); then
	varfreq=0.1
fi


varscan pileup2snp ${input} --p-value 0.01 --min-var-freq ${varfreq}  > ${output_SNP}  
varscan pileup2indel ${input} --p-value 0.01 --min-var-freq ${varfreq} > ${output_INDEL}
varscan pileup2cns ${input} --p-value 0.01 --min-var-freq ${varfreq} > ${output_CNS}