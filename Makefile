CCX = g++
CFLAGS= -Wall

all: main

main:
	$(CCX) -O3 -std=c++17 -o RUN benchmarking.cpp diameter_annealing/*.cpp graph_generators/*.cpp $(CFLAGS)

