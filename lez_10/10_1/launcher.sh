#!bin/bash

temp=$1

for((i=0;i<${temp};i++));
do
	./10_1.exe ${i}
done;
