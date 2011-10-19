all: hk hk_viz
hk_viz: mongoose.o hk_viz.cpp held_karp.h utils.h heuristics.h
	g++ hk_viz.cpp mongoose.o -o hk_viz -Wall -O3 -pthread -ldl -lgurobi45 -lemon
mongoose.o: mongoose.h mongoose.c
	gcc -c mongoose.c -o mongoose.o -Wall -O3
hk: hk.cpp held_karp.h utils.h heuristics.h
	g++	hk.cpp -o hk -Wall -O3 -lgurobi45 -lemon
clean:
	rm -rf hk hk_viz mongoose.o
