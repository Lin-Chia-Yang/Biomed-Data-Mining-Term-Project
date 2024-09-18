all:motif.cpp
	g++ -fopenmp -O3 motif.cpp -o motif.out
clean:
	rm -f motif.out
