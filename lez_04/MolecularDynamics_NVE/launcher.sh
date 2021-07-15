#!bin/bash

therm=$1

cp config.fcc config.0
cp input.solid input.dat
./MolDyn_NVE.exe 0 0 0
    for ((j=0;j<${therm};j++));
    do
		./MolDyn_NVE.exe 1 0 0
	done
./MolDyn_NVE.exe 1 1 0
cp config.fcc config.0
cp input.liquid input.dat
./MolDyn_NVE.exe 0 0 1
for ((j=0;j<${therm};j++));
do
    ./MolDyn_NVE.exe 1 0 1
done
./MolDyn_NVE.exe 1 1 1
cp config.final config.0
cp input.gas input.dat
./MolDyn_NVE.exe 0 0 2
for ((j=0;j<${therm};j++));
do
    ./MolDyn_NVE.exe 1 0 2
done
./MolDyn_NVE.exe 1 1 2
