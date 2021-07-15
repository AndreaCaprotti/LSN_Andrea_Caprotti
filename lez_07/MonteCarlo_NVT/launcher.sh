#!bin/bash

therm=$1

cp config.fcc config.0
cp input.solid input.dat
./Monte_Carlo_NVT.exe solid 1000 1 0 0
    for ((j=0;j<${therm};j++));
    do
		./Monte_Carlo_NVT.exe solid 1000 1 1 0
	done
./Monte_Carlo_NVT.exe solid 500000 1 1 1
cp config.fcc config.0
cp input.liquid input.dat
./Monte_Carlo_NVT.exe liquid 1000 1 0 0
for ((j=0;j<${therm};j++));
do
    ./Monte_Carlo_NVT.exe liquid 1000 1 1 0
done
./Monte_Carlo_NVT.exe liquid 500000 1 1 1
cp config.final config.0
cp input.gas input.dat
./Monte_Carlo_NVT.exe gas 1000 1 0 0
for ((j=0;j<${therm};j++));
do
    ./Monte_Carlo_NVT.exe gas 1000 1 1 0
done
./Monte_Carlo_NVT.exe gas 500000 1 1 1
