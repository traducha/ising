#include <igraph.h>

#ifndef TOOLS_H
#define TOOLS_H


void print_array_int(int *array, int n);
/*
 * Prints array of size n
 */


void print_array_double(double *array, int n);
/*
 * Prints array of size n
 */
 

void print_vector_igraph(igraph_vector_t *v);
/*
 * Prints igraph vector
 */


int index_of_value_int(int *array, int n, int value);
/*
 * Returns index of the first value in the array.
 * If value is not present in the array, returns -1.
 */


int index_of_value_vector(igraph_vector_t *v, int n, int value);
/*
 * Returns index of the first value in the vactor.
 * If value is not present in the vactor, returns -1.
 */


int sum_array_int(int *array, int n);
/*
 * Sums array's elements
 */


void save_to_csv_int(char *file_name, int *array, int n);
/*
 * Saves int array to .csv file in a manner:
 * index;0;1;2;3;4...
 * value;<value1>;<value2>;<value3>;<value4>...
 */
 
 
 void save_to_csv_double(char *file_name, double *array, int n);
/*
 * Saves double array to .csv file in a manner:
 * index;0;1;2;3;4...
 * value;<value1>;<value2>;<value3>;<value4>...
 */
 
 
void save_to_csv_xystd_double(char *file_name, double *x, double *y, double *std, int n);
/*
 * Saves double array1 and array2 to .csv file in a manner:
 * x;<value1>;<value2>;<value3>;<value4>...
 * y;<value1>;<value2>;<value3>;<value4>...
 * std;<value1>;<value2>;<value3>;<value4>...
 */
 
 
void initial_graph(igraph_t *graph, int *spins, int nodes, int edges);
/*
 * Generates initial ER graph for simulation with random spins.
 */


double energy_change(int spin_change, igraph_vector_t *neigs, int *spins, double J, double h);
/*
 * Calculates the change in hamiltonian after given spin_change,
 * spin_change must be equal <new_state - old_state>
 */


int new_vertex(int g_size, int exclude1, int exclude2, igraph_vector_t *exclude);
/*
 * Randomly selects new vertex
 * different than exclude1, exclude2 and not present in exclude array,
 * with equal probability for every node
 */


int * draw_new_edge(igraph_t *graph, int g_size, int exclude1, int exclude2);
/*
 * Randomly selects new edge (different than the current one),
 * with equal probability for every node,
 * without excluded nodes, without loops, without repeting existing edge.
 */
 
 
 double avg_double(double *array, int len);
/*
 * Calculates the average value for a given double array
 */


double std_double(double *array, int len, double avg);
/*
 * Calculates the standard deviation for a given double array
 */


double avg_int(int *array, int len);
/*
 * Calculates the average value for a given int array
 */


double std_int(int *array, int len, double avg);
/*
 * Calculates the standard deviation for a given int array
 */


double compute_energy(int g_size, igraph_vector_t *esge_list, int *spins, double J, double h);
/*
 * Calculates the total energy of the graph with given spins and connections,
 * esge_list must contain pairs of connected nodes
 */
 
 
int count_incompatible_links(igraph_vector_t *esge_list, int *spins);
/*
 * Counts the number of links connecting nodes with different spins.
 */


#endif  /* TOOLS_H */