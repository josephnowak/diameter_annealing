CCX = g++
CFLAGS= -Wall

all: main

main:
	$(CCX) -o RUN test_imports.cpp diameter_annealing/*.cpp $(CFLAGS)

