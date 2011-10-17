all: held_karp
held_karp: held_karp.cpp held_karp.h utils.h
	g++	held_karp.cpp -o held_karp -Wall -O3 -lgurobi45 -lemon
clean:
	rm -rf held_karp
