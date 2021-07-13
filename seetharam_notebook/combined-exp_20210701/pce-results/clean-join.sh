#!/bin/bash
f=$1
sed 's/,/\t/g' $f |sed 's/"//g' |sed 's/^Log10PValue/CellType\tLog10PValue/g' > .temp1
cat *.txt |sed 's/,/\t/g' |sed 's/"//g' |grep -v "tegenes" | sed '1 i CellType\tGenes' > .temp2
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2;next}{ print $0, a[$1]}' .temp2 .temp1 > ../${f%.*}_with-genes.tsv
