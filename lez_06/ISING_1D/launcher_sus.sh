#!bin/bash

bias=$1
therm=$2
field=$3
meas=$4

for ((i=1;i<${bias};i++));
do
	./Monte_Carlo_ISING_1D.exe  0 1 0 0 ${therm} ${i} ${meas}
	for ((j=0;j<${therm};j++));
	do
		./Monte_Carlo_ISING_1D.exe 0 1 1 0 ${therm} ${i} ${meas}
	done
	./Monte_Carlo_ISING_1D.exe 0 1 1 1 ${therm} ${i} ${meas}
done
