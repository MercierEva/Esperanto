count=$(zcat $1 | awk '(NR%4==2)' | wc -l )

extract=$(( ${3}*98/${count} ))
filtlong -p $extract $1 | gzip > $2

