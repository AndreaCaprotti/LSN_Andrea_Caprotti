CC = g++
CFLAGS = -Wall -O3 --std=c++11
DEBUG= -g -gdwarf-2

LIB_PATH:= ../../libraries/
EXECUTABLES:= 10_1
INCLUDES:= $(LIB_PATH)random.h $(LIB_PATH)random_walk.h $(LIB_PATH)block_stat.h $(LIB_PATH)chromosome.h

$(EXECUTABLES).exe : $(EXECUTABLES).o random.o
	$(CC) $(DEBUG) random.o $(EXECUTABLES).o -o $(EXECUTABLES).exe

$(EXECUTABLES).o : $(EXECUTABLES).cpp
	$(CC) $(DEBUG) -c $(EXECUTABLES).cpp -o $(EXECUTABLES).o $(CFLAGS)

random.o : $(LIB_PATH)random.cpp $(LIB_PATH)random.h
	$(CC) -c $(LIB_PATH)random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o $(EXECUTABLES).exe seed.out

clean_file:
	rm *.txt *.dat
