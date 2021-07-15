CC = g++
CFLAGS = -Wall -O3 --std=c++11

LIB_PATH:= ../../libraries/
INCLUDES:= $(LIB_PATH)random.h Monte_Carlo_NVT.h

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe

c_gas :
	rm output.*.gas

c_liquid :
	rm output.*.liquid
	
c_solid :
	rm output.*.solid
