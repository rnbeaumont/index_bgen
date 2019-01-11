# makefile for index_bgen
CC=g++
EX_FLAGS=-std=c++11 -lz -lboost_iostreams -O3 -Wall
SOURCES=main.cpp general_functions.cpp index.cpp
EXECUTABLE=index_bgen

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(SOURCES) $(FLAGS) -o $@ $(EX_FLAGS)

clean:
	rm $(OBJECTS) $(EXECUTABLE)
