input=$2
threads=$3
output=$4
model=$5
cluster_bigger=$(ls -S ${1}* | head -1 || true)
medaka_consensus -i $cluster_bigger -d $input -t $threads -o $output -m $model || true

