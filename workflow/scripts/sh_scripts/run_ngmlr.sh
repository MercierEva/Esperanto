threads=$1
ref=$2
cluster_dir=$3
output=$4

reads=$(ls -S ${cluster_dir}* | head -1)

ngmlr -x ont -t ${threads} -r ${ref} -q ${reads} -o ${output}