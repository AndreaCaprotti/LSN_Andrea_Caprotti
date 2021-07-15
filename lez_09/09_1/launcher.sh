#!bin/bash

tot=$1

for ((i=0;i<${tot};i++));
do		    
	for ((j=0;j<${tot};j++));
	do
		./09_1.exe ${i} 10 ${j} ${tot}
	done
done	
