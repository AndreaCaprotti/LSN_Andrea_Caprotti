#!bin/bash

simulation_no=$1
therm=$2
field=$3
meas=$4

for ((i=0;i<${simulation_no};i++));
do
	./Monte_Carlo_ISING_1D.exe  ${i} ${simulation_no} 0 0 ${field} ${meas}
	for ((j=0;j<${therm};j++));
	do
		./Monte_Carlo_ISING_1D.exe  ${i} ${simulation_no} 1 0 ${field} ${meas}
	done
	./Monte_Carlo_ISING_1D.exe ${i} ${simulation_no} 1 1 ${field} ${meas}
done
