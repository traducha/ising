ising: main.c tools.c
	gcc -O3 -Wall -I/usr/local/include/igraph -L/usr/local/lib -o main main.c tools.c -ligraph -lm
