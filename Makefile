all:
	g++	main.cpp -o tsp -Wall -O3 -lgurobi45 -lemon
clean:
	rm -rf tsp
