tsp: main.cpp held_karp.h utils.h
	g++	main.cpp -o tsp -Wall -O3 -lgurobi45 -lemon
clean:
	rm -rf tsp
