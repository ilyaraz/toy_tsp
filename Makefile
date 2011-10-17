all: hk
hk: hk.cpp held_karp.h utils.h
	g++	hk.cpp -o hk -Wall -O3 -lgurobi45 -lemon
clean:
	rm -rf hk
