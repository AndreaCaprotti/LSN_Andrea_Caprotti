#!bin/bash
  
mu_start=$1
mu_end=$2
s_start=$3
s_end=$4
tot=$5

for ((i=0;i<${tot};i++));
do
	for ((j=0;j<${tot};j++));
	do
		./08_1.exe ${mu_start} ${i} ${mu_end} ${s_start} ${j} ${s_end} ${tot}
	done
done	
