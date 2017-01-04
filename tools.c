#include <stdlib.h>
#include <igraph.h>
#include "tools.h"


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


void print_array_double(double *array, int n)
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
            printf("%f, ", *(array + i));  // the same as array[i]
            continue;
        }
        printf("%f ", *(array + i));
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


int index_of_value_int(int *array, int n, int value)
/*
 * Returns index of the first value in the array.
 * If value is not present in the array, returns -1.
 */
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (array[i] == value)
        {
            return i;
        }
    }
    return -1;
}


int index_of_value_vector(igraph_vector_t *v, int n, int value)
/*
 * Returns index of the first value in the vactor.
 * If value is not present in the vactor, returns -1.
 */
{
    int i;
    for (i = 0; i < n; i++)
    {
        if ((long int) VECTOR(*v)[i] == value)
        {
            return i;
        }
    }
    return -1;
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


void save_to_csv_int(char *file_name, int *array, int n)
/*
 * Saves int array to .csv file in a manner:
 * index;0;1;2;3;4...
 * value;<value1>;<value2>;<value3>;<value4>...
 */
{    
    FILE *file;
    file = fopen(file_name, "w");
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


void save_to_csv_double(char *file_name, double *array, int n)
/*
 * Saves double array to .csv file in a manner:
 * index;0;1;2;3;4...
 * value;<value1>;<value2>;<value3>;<value4>...
 */
{    
    FILE *file;
    file = fopen(file_name, "w");
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
            fprintf(file, "%f;", array[i]);
            continue;
        }
        fprintf(file, "%f", array[i]);
    }
    fclose(file);
    return;
}


void save_to_csv_xystd_double(char *file_name, double *x, double *y, double *std, int n)
/*
 * Saves double array1 and array2 to .csv file in a manner:
 * x;<value1>;<value2>;<value3>;<value4>...
 * y;<value1>;<value2>;<value3>;<value4>...
 * std;<value1>;<value2>;<value3>;<value4>...
 */
{    
    FILE *file;
    file = fopen(file_name, "w");
    if (file == NULL)
    {
        perror ("Cannot open a file");
        return;
    }
    
    fprintf(file, "x;");
    int i;
    for(i = 0; i < n; i++)
    {
        if (i < n - 1)
        {
            fprintf(file, "%f;", x[i]);
            continue;
        }
        fprintf(file, "%f\n", x[i]);
    }
    
    fprintf(file, "y;");
    for(i = 0; i < n; i++)
    {
        if (i < n - 1)
        {
            fprintf(file, "%f;", y[i]);
            continue;
        }
        fprintf(file, "%f\n", y[i]);
    }
    
    fprintf(file, "std;");
    for(i = 0; i < n; i++)
    {
        if (i < n - 1)
        {
            fprintf(file, "%f;", std[i]);
            continue;
        }
        fprintf(file, "%f", std[i]);
    }
    fclose(file);
    return;
}


void initial_graph(igraph_t *graph, int *spins, int nodes, int edges)
/*
 * Generates initial ER graph for simulation with random spins.
 */
{
    igraph_erdos_renyi_game(graph, IGRAPH_ERDOS_RENYI_GNM, nodes, edges,
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
    return;
}


double energy_change(int spin_change, igraph_vector_t *neigs, int *spins, double J, double h)
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


int new_vertex(int g_size, int exclude1, int exclude2, igraph_vector_t *exclude)
/*
 * Randomly selects new vertex
 * different than exclude1, exclude2 and not present in exclude array,
 * with equal probability for every node
 */
{
    double r = rand();
    int new_v = (( (g_size - 1) * r / RAND_MAX) + 0.5);
    if (new_v == exclude1 || new_v == exclude2 ||
                index_of_value_vector(exclude, igraph_vector_size(exclude), new_v) != -1)
    {
        new_v = new_vertex(g_size, exclude1, exclude2, exclude);
    }
    return new_v;
}


int * draw_new_edge(igraph_t *graph, int g_size, int exclude1, int exclude2)
/*
 * Randomly selects new edge (different than the current one),
 * with equal probability for every node,
 * without excluded nodes, without loops, without repeting existing edge.
 */
{
    double r = rand(), w = rand();
    static int new_v[2];
    new_v[0] = (( (g_size - 1) * r / RAND_MAX) + 0.5);
    new_v[1] = (( (g_size - 1) * w / RAND_MAX) + 0.5);
    if (new_v[0] == new_v[1] || new_v[0] == exclude1 || new_v[1] == exclude1 ||
                new_v[0] == exclude2 || new_v[1] == exclude2)
    {
        int *p;
        p = draw_new_edge(graph, g_size, exclude1, exclude2);
        new_v[0] = *p;  // the same as p[0]
        new_v[1] = *(p + 1);  // the same as p[1]
    }
    else
    {
        igraph_integer_t eid;
        igraph_get_eid(graph, &eid, (igraph_integer_t) new_v[0],
                                    (igraph_integer_t) new_v[1], IGRAPH_UNDIRECTED, 0);
        if (((int) eid) != -1)
        {
            int *p;
            p = draw_new_edge(graph, g_size, exclude1, exclude2);
            new_v[0] = *p;  // the same as p[0]
            new_v[1] = *(p + 1);  // the same as p[1]
        }
    }
    return new_v;
}


double avg_double(double *array, int len)
/*
 * Calculates the average value for a given double array
 */
{
    double average = 0.0;
    int i;
    for (i = 0; i < len; i++)
    {
        average += (double) array[i];
    }
    average = average / len;
    
    return average;
}


double std_double(double *array, int len, double avg)
/*
 * Calculates the standard deviation for a given double array
 */
{
    double std = 0.0;
    int i;  
    for (i = 0; i < len; i++)
    {
        std += pow((double) array[i] - avg, 2.0);
    }
    std = sqrt(std / len);
    
    return std;
}


double avg_int(int *array, int len)
/*
 * Calculates the average value for a given int array
 */
{
    double average = 0.0;
    int i;
    for (i = 0; i < len; i++)
    {
        average += (double) array[i];
    }
    average = average / len;
    
    return average;
}


double std_int(int *array, int len, double avg)
/*
 * Calculates the standard deviation for a given int array
 */
{
    double std = 0.0;
    int i;  
    for (i = 0; i < len; i++)
    {
        std += pow((double) array[i] - avg, 2.0);
    }
    std = sqrt(std / len);
    
    return std;
}


double compute_energy(int g_size, igraph_vector_t *esge_list, int *spins, double J, double h)
/*
 * Calculates the total energy of the graph with given spins and connections,
 * esge_list must contain pairs of connected nodes
 */
{
    double energy = 0.0;
    int i, from, to, n = igraph_vector_size(esge_list);
    for (i = 0; i < n; i += 2)
    {
        from = VECTOR(*esge_list)[i];
        to = VECTOR(*esge_list)[i+1];
        energy += - (J * spins[from] * spins[to]);
    }
    if (h != 0.0)
    {
        for (i = 0; i < g_size; i++)
        {
            energy += - h * spins[i];
        }
    }
    return energy;
}
