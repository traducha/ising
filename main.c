#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <igraph.h>
#include "graph_tools.h"


const int T = 10;  // time steps in main loop
const int N = 10;  // number of nodes
const double M = 20;  // number of edges


void print_array_int(int* array, int n)
/*
 * Prints array
 */
{
    int i;
    for(i = 0; i < n; i++)
    {
        printf("%d\n", array[i]);
    }
}


int sum_array_int(int array[], int n)
/*
 * Sums array's elements
 */
{
    int i, res = 0;
    for(i = 0; i < n; i++)
    {
        res += array[i];
    }
    return res;
}


void save_to_csv_int(int* mag, int n)
/*
 * Saves array to .csv file in a manner (index;array[index])
 */
{
    FILE *file;
    file = fopen("magnetization.csv", "w");
    int i;
    for(i = 0; i < n; i++)
    {
        fprintf(file, "%d;%d\n", i, mag[i]);
    }
    fclose(file);
}


igraph_t initial_graph(int spins[N])
/*
 * Generates initial ER graph for simulation.
 */
{
    igraph_t graph;
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, N, M,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
                            
    int i;
    double r;
    srand(time(NULL));
    
    for (i = 0; i < N; i++)
    {
        r = rand();
        if (r / RAND_MAX < 0.5)
        {
            spins[i] = 1;
        }
        else
        {
            spins[i] = -1;
        }
    }
    return graph;
}


void main_loop(igraph_t graph, int spins[N], int* mag, int times)
/*
 * 
 */
{
    int i, v_index, e_index;
    double r;
    srand(time(NULL));
    
    for (i = 0; i < N; i++)
    {
        // select vertex and edge
        r = rand();
        v_index = (int)((N * r / RAND_MAX) + 0.5);
        r = rand();
        e_index = (int)((M * r / RAND_MAX) + 0.5);
        
        // TODO: do metropolis
        spins[v_index] = -spins[v_index];        
        
        // compute magnetization
        mag[i] = sum_array_int(spins, N);
    }
}


int main(void)
{
    int spins[N], mag[T];
    igraph_t graph = initial_graph(spins);
    
    main_loop(graph, spins, mag, T);
    print_array_int(mag, T);
    
    save_to_csv_int(mag, T);
    
    igraph_destroy(&graph);
    return 0;
}
