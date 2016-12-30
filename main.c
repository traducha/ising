#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <igraph.h>
#include "graph_tools.h"


const int T = 1000000;  // time steps in main loop
const int N = 10000;  // number of nodes
const int M = 20000;  // number of edges <k>=2*M/N
const double J = 1.0;  // J in hamiltonian
const double h = 0.0;  // h in hamiltonian
const double B = 1.0;  // 1/kbT inverse of the temperature


void print_array_int(int *array, int n)
/*
 * Prints array of size n
 */
{
    printf("[ ");
    long int i;
    for(i = 0; i < n; i++)
    {
        if (i < n - 1)
        {
            printf("%d, ", *(array + i));  // the same as array[i]
            continue;
        }
        printf("%d ", *(array + i));
    }
    printf("]\n");
    return;
}


void print_vector_igraph(igraph_vector_t *v)
/*
 * Prints igraph vector
 */
{
    printf("[ ");
    long int i, n = igraph_vector_size(v);
    for (i = 0; i < n; i++) {
        if (i < n - 1)
        {
            printf("%li, ", (long int) VECTOR(*v)[i]);
            continue;
        }
        printf("%li ", (long int) VECTOR(*v)[i]);
    }
    printf("]\n");
    return;
}


int sum_array_int(int *array, int n)
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


void save_to_csv_int(char *fname, int *array, int n)
/*
 * Saves int array to .csv file in a manner <index;array[index]>
 */
{
    FILE *file;
    file = fopen(fname, "w");
    if (file == NULL)
    {
        perror ("Cannot open a file");
        return;
    }
    
    fprintf(file, "index;");
    int i;
    for(i = 0; i < n; i++)
    {
        if (i < n - 1)
        {
            fprintf(file, "%d;", i);
            continue;
        }
        fprintf(file, "%d\n", i);
    }
    
    fprintf(file, "value;");
    for(i = 0; i < n; i++)
    {
        if (i < n - 1)
        {
            fprintf(file, "%d;", array[i]);
            continue;
        }
        fprintf(file, "%d", array[i]);
    }
    fclose(file);
    return;
}


igraph_t initial_graph(int *spins, int nodes, int edges)
/*
 * Generates initial ER graph for simulation.
 */
{
    igraph_t graph;
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, nodes, edges,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    
    srand(time(NULL));
    double r;
    int i;
    for (i = 0; i < nodes; i++)
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


double energy_change(int spin_change, igraph_vector_t *neigs, int *spins)
/*
 * Calculates the change in hamiltonian after given spin_change,
 * spin_change must be equal <new_state - old_state>
 */
{
    double delta = 0.0;
    int i, neig_index;
    for (i = 0; i < igraph_vector_size(neigs); i++)
    {
        neig_index = VECTOR(*neigs)[i];
        delta += spins[neig_index];
    }
    delta = - (J * delta + h) * spin_change;
    return delta;
}


int new_vertex(int g_size, int current_v)
/*
 * Randomly selects new vertex (different than the current one),
 * with equal probability for every node
 */
{
    double r = rand();
    int new_v = (( (g_size - 1) * r / RAND_MAX) + 0.5);
    if (new_v == current_v)
    {
        new_v = new_vertex(g_size, current_v);
    }
    return new_v;
}


void main_loop(igraph_t *graph, int *spins, int nodes, int edges, int *mag, int times)
/*
 * Main loop of the algorithm, currently with switching only one end on a node
 */
{
    srand(time(NULL));
    igraph_vector_t neigs, edges_vector;
    igraph_vector_init(&neigs, 1);
    igraph_vector_init(&edges_vector, 1);
    igraph_integer_t  e_from, e_to, v_to_new;
    igraph_es_t es;
    double r, delta, delta_v;
    int v_index, e_index, v_to_detach, v_to_keep;
    
    int i;
    for (i = 0; i < times; i++)
    {
        // select vertex and edge
        r = rand();
        v_index = (int)(((nodes - 1) * r / RAND_MAX) + 0.5);
        r = rand();
        e_index = (int)(((edges - 1) * r / RAND_MAX) + 0.5);

        // spin switching
        igraph_neighbors(graph, &neigs, v_index, IGRAPH_ALL);
        delta = energy_change(-2 * spins[v_index], &neigs, spins);  // -2 * spin is the optional energy change
        if (delta <= 0.0)
        {
            spins[v_index] = -spins[v_index];
        }
        else
        {
            r = rand();
            if (r / RAND_MAX < exp(- B * delta))
            {
                spins[v_index] = -spins[v_index];
            }
            
        }
        
        // edge rewiring
        igraph_edge(graph, e_index, &e_from, &e_to);
        r = rand();
        if (r / RAND_MAX < 0.5)
        {
            v_to_detach = e_from;
            v_to_keep = e_to;
        }
        else
        {
            v_to_detach = e_to;
            v_to_keep = e_from;
        }
        v_to_new = new_vertex(nodes, v_to_detach);
        
        delta_v = - J * spins[v_to_keep] * (spins[v_to_new] - spins[v_to_detach]);
        if (delta_v <= 0.0)
        {
            //delete edge e_from<->e_to
            igraph_es_pairs_small(&es, IGRAPH_UNDIRECTED, e_from, e_to, -1);
            igraph_delete_edges(graph, es);
            igraph_es_destroy(&es);
            
            //create edge v_to_keep<->v_to_new
            igraph_vector_resize(&edges_vector, 2);
            VECTOR(edges_vector)[0] = v_to_keep;
            VECTOR(edges_vector)[1] = v_to_new;
            igraph_add_edges(graph, &edges_vector, 0);
            //printf("(%d, %d) -> (%d, %d) \n\n", e_from, e_to, v_to_keep, v_to_new);
        }
        else
        {
            r = rand();
            if (r / RAND_MAX < exp(- B * delta_v))
            {
                //delete edge e_from<->e_to
                igraph_es_pairs_small(&es, IGRAPH_UNDIRECTED, e_from, e_to, -1);
                igraph_delete_edges(graph, es);
                igraph_es_destroy(&es);
                
                //create edge v_to_keep<->v_to_new
                igraph_vector_resize(&edges_vector, 2);
                VECTOR(edges_vector)[0] = v_to_keep;
                VECTOR(edges_vector)[1] = v_to_new;
                igraph_add_edges(graph, &edges_vector, 0);
                //printf("(%d, %d) -> (%d, %d) \n\n", e_from, e_to, v_to_keep, v_to_new);
            }
        }
        
//        printf("(%d %li) -> (%li ) \n\n", e_index, e_from, e_to);
//        igraph_get_edgelist(graph, &edges_vector, 0);
//        print_vector_igraph(&edges_vector);
        
        // compute magnetization
        mag[i] = sum_array_int(spins, nodes);
    }
    
    igraph_vector_destroy(&neigs);
    igraph_vector_destroy(&edges_vector);
    return;
}


int main(void)
{
    clock_t t_start, t_end;
    int *spins, *mag;
    spins = calloc(N, sizeof(int));
    mag = calloc(T, sizeof(int));
    igraph_t graph = initial_graph(spins, N, M);
    
    t_start = clock();
    main_loop(&graph, spins, N, M, mag, T);
    t_end = clock();
    
    save_to_csv_int("magnetization.csv", mag, T);
    printf("Main loop finished in %f s.\n",
        (double)(t_end - t_start) / CLOCKS_PER_SEC);
    
    free(mag);
    free(spins);
    igraph_destroy(&graph);
    return 0;
}
