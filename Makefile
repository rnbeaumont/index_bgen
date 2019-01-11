# makefile for poe program
CC=g++
BOOST=/mnt/Data10/rnb203/programs/boost_1_59_0/boost_1_59_0/lib/
EX_FLAGS=-std=c++11 -lz -lboost_iostreams -O3 -Wall
OB_FLAGS= -c
SOURCES=main.cpp general_functions.cpp index.cpp
EXECUTABLE=index_bgen

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(SOURCES) $(FLAGS) -o $@ $(EX_FLAGS)

clean:
	rm $(OBJECTS) $(EXECUTABLE)
