#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <igraph.h>
#include "tools.h"


const int T = 100000;  // thermalization time
const int SAMPLE_TIME = 200000;  // time window to average mag ang energy over (after thermalization)
const int N = 500;  // number of nodes
const int M = 1000;  // number of edges <k>=2*M/N, c=M/N
const double J = 1.0;  // J in hamiltonian
const double h = 0.0;  // h in hamiltonian
const double B = 1.0;  // 1/kbT inverse of the temperature
const double MIN_TEMP = 0.99;  // min temperature
const double MAX_TEMP = 10.0;  // max temperature
const int TEMP_STEPS = 100;  // number of values for temperature


void algorithm_one(igraph_t *graph, int *spins, int nodes, int edges,
                        int *mag, double *energy, int times, double beta)
/*
 * Metropolis algorithm of Ising model with <times> time steps and saving magnetization and energy.
 * Currently with switching only one end of a node.
 */
{
    srand(time(NULL));
    igraph_vector_t neigs, edges_vector;
    igraph_vector_init(&neigs, 1);
    igraph_vector_init(&edges_vector, 1);
    igraph_integer_t  e_from, e_to, v_to_new;
    igraph_es_t es;
    double r, delta;
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
        delta = energy_change(-2 * spins[v_index], &neigs, spins, J, h);  // -2 * spin is the possible energy change
        if (delta <= 0.0)
        {
            spins[v_index] = -spins[v_index];
        }
        else
        {
            r = rand();
            if (r / RAND_MAX < exp(- beta * delta))
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
        
        // draw new vertex to switch to
        igraph_neighbors(graph, &neigs, v_to_keep, IGRAPH_ALL);
        v_to_new = new_vertex(nodes, v_to_detach, v_to_keep, &neigs);
        
        delta = - J * spins[v_to_keep] * (spins[v_to_new] - spins[v_to_detach]);
        r = rand();
        if (delta <= 0.0 || r / RAND_MAX < exp(- beta * delta))
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
        }
        
//        printf("(%d %li) -> (%li ) \n\n", e_index, e_from, e_to);
//        igraph_get_edgelist(graph, &edges_vector, 0);
//        print_vector_igraph(&edges_vector);
        
        // compute magnetization and internal energy
        igraph_get_edgelist(graph, &edges_vector, 0);
        energy[i] = compute_energy(nodes, &edges_vector, spins, J, h);
        mag[i] = sum_array_int(spins, nodes);
    }
    
    igraph_vector_destroy(&neigs);
    igraph_vector_destroy(&edges_vector);
    return;
}


void algorithm_two(igraph_t *graph, int *spins, int nodes, int edges,
                        int *mag, double *energy, int times, double beta)
/*
 * Metropolis algorithm of Ising model with <times> time staps and saving magnetization and energy.
 * This version implements switching both ends of an edge.
 * About 0.5% faster than first version.
 */
{
    srand(time(NULL));
    igraph_vector_t neigs, edges_vector;
    igraph_vector_init(&neigs, 1);
    igraph_vector_init(&edges_vector, 1);
    igraph_integer_t e_from, e_to;
    igraph_es_t es;
    double r, delta;
    int v_index, e_index, *new_e;
    
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
        delta = energy_change(-2 * spins[v_index], &neigs, spins, J, h);  // -2 * spin is the possible energy change
        if (delta <= 0.0)
        {
            spins[v_index] = -spins[v_index];
        }
        else
        {
            r = rand();
            if (r / RAND_MAX < exp(- beta * delta))
            {
                spins[v_index] = -spins[v_index];
            }
            
        }
        
        // edge rewiring
        igraph_edge(graph, e_index, &e_from, &e_to);
        
        // draw new ends for the edge
        new_e = draw_new_edge(graph, nodes, (int) e_from, (int) e_to);
        
        delta = - J * (spins[new_e[0]] * spins[new_e[1]] - spins[e_from] * spins[e_to]);
        r = rand();
        if (delta <= 0.0 || r / RAND_MAX < exp(- beta * delta))
        {
            //delete edge e_from<->e_to
            igraph_es_pairs_small(&es, IGRAPH_UNDIRECTED, e_from, e_to, -1);
            igraph_delete_edges(graph, es);
            igraph_es_destroy(&es);
            
            //create edge new_e[0]<->new_e[1]
            igraph_vector_resize(&edges_vector, 2);
            VECTOR(edges_vector)[0] = new_e[0];
            VECTOR(edges_vector)[1] = new_e[1];
            igraph_add_edges(graph, &edges_vector, 0);
        }
        
        // compute magnetization and internal energy
        igraph_get_edgelist(graph, &edges_vector, 0);
        energy[i] = compute_energy(nodes, &edges_vector, spins, J, h);
        mag[i] = sum_array_int(spins, nodes);
    }
    
    igraph_vector_destroy(&neigs);
    igraph_vector_destroy(&edges_vector);
    return;
}


int compare_thermalization(void)
/*
 * Function computing magnetization and energy in time
 * to check and compare thermalization times.
 */
{
    clock_t t_start, t_end;
    int *spins, *mag;
    double *energy, average, std;
    char *file_name;
    igraph_t graph;
    
    int i, R = 2;
    for (i = 0; i < R; i++)
    {
        spins = calloc(N, sizeof(int));
        mag = calloc(T, sizeof(int));
        energy = malloc(T * sizeof(double));
        initial_graph(&graph, spins, N, M);
        t_start = clock();
        algorithm_one(&graph, spins, N, M, mag, energy, T, B);
        t_end = clock();
        
        file_name = (char*)malloc(30 * sizeof(char));
        sprintf(file_name, "res/mag1_%d.csv", i);
        save_to_csv_int(file_name, mag, T);
        free(file_name);
        file_name = (char*)malloc(30 * sizeof(char));
        sprintf(file_name, "res/energy1_%d.csv", i);
        save_to_csv_double(file_name, energy, T);
        
        printf("Main loop 1 finished in %f s.\n",
            (double)(t_end - t_start) / CLOCKS_PER_SEC);
            
        average = avg_int(mag, T);
        printf("mag: %f\n", average);
        std = std_int(mag, T, average);
        printf("std: %f\n", std);
        
        average = avg_double(energy, T);
        printf("energy: %f\n", average);
        std = std_double(energy, T, average);
        printf("std: %f\n", std);
            
        free(mag);
        free(spins);
        free(energy);
        free(file_name);
        igraph_destroy(&graph);
    }
    
    
    for (i = 0; i < R; i++)
    {
        spins = calloc(N, sizeof(int));
        mag = calloc(T, sizeof(int));
        energy = malloc(T * sizeof(double));
        initial_graph(&graph, spins, N, M);    
        t_start = clock();
        algorithm_two(&graph, spins, N, M, mag, energy, T, B);
        t_end = clock();

        file_name = (char*)malloc(30 * sizeof(char));
        sprintf(file_name, "res/mag2_%d.csv", i);
        save_to_csv_int(file_name, mag, T);
        free(file_name);
        file_name = (char*)malloc(30 * sizeof(char));
        sprintf(file_name, "res/energy2_%d.csv", i);
        save_to_csv_double(file_name, energy, T);
        
        printf("Main loop 2 finished in %f s.\n",
            (double)(t_end - t_start) / CLOCKS_PER_SEC);
            
        average = avg_int(mag, T);
        printf("mag: %f\n", average);
        std = std_int(mag, T, average);
        printf("std: %f\n", std);
        
        average = avg_double(energy, T);
        printf("energy: %f\n", average);
        std = std_double(energy, T, average);
        printf("std: %f\n", std);
        
        free(mag);
        free(spins);
        free(energy);
        free(file_name);
        igraph_destroy(&graph);
    }
    return 0;
}


int phase_diagram_one(void)
/*
 * Computes magnetization and energy (after thermalization)
 * for different temperatures (B_ARRAY - with inverse of temperature)
 */
{
    clock_t t_start, t_end;
    int *spins, *mag;
    double *energy, tmp, temp[TEMP_STEPS];
    igraph_t graph;
    
    double mag_avg[TEMP_STEPS], mag_std[TEMP_STEPS];  // results
    double energy_avg[TEMP_STEPS], energy_std[TEMP_STEPS];  // results
    
    int i;
    for (i = 0; i < TEMP_STEPS; i++)
    {
        t_start = clock();
        spins = calloc(N, sizeof(int));
        initial_graph(&graph, spins, N, M);
        
        temp[i] = MIN_TEMP + ((MAX_TEMP - MIN_TEMP) / TEMP_STEPS) * i;       
        
        // first thermalization
        mag = calloc(T, sizeof(int));
        energy = malloc(T * sizeof(double));
        algorithm_one(&graph, spins, N, M, mag, energy, T, 1.0/temp[i]);  // TODO additional arg to omit calculationg mag and energy
        free(mag);
        free(energy);
        
        // actual calculations
        mag = calloc(SAMPLE_TIME, sizeof(int));
        energy = malloc(SAMPLE_TIME * sizeof(double));
        algorithm_one(&graph, spins, N, M, mag, energy, SAMPLE_TIME, 1.0/temp[i]);

        // computing results
        tmp = avg_int(mag, SAMPLE_TIME);
        mag_avg[i] = abs(tmp);  // take absolute value of magnetization
        mag_std[i] = std_int(mag, SAMPLE_TIME, tmp);
        
        energy_avg[i] = avg_double(energy, SAMPLE_TIME);
        energy_std[i] = std_double(energy, SAMPLE_TIME, energy_avg[i]);
        
        free(mag);
        free(spins);
        free(energy);
        igraph_destroy(&graph);
        t_end = clock();
        printf("%d Computing energy and magnetization for T=%f finished in %f s.\n",
            i, temp[i], (double)(t_end - t_start) / CLOCKS_PER_SEC);
    }
    
    char file_name1[] = "res/energy_vs_B_1.csv";
    save_to_csv_xystd_double(file_name1, temp, energy_avg, energy_std, TEMP_STEPS);
    
    char file_name2[] = "res/mag_vs_B_1.csv";
    save_to_csv_xystd_double(file_name2, temp, mag_avg, mag_std, TEMP_STEPS);
    return 0;
}


int phase_diagram_two(void)
/*
 * Computes magnetization and energy (after thermalization)
 * for different temperatures (B_ARRAY - with inverse of temperature)
 */
{
    clock_t t_start, t_end;
    int *spins, *mag;
    double *energy, tmp, temp[TEMP_STEPS];
    igraph_t graph;
    
    double mag_avg[TEMP_STEPS], mag_std[TEMP_STEPS];  // results
    double energy_avg[TEMP_STEPS], energy_std[TEMP_STEPS];  // results
    
    int i;
    for (i = 0; i < TEMP_STEPS; i++)
    {
        t_start = clock();
        spins = calloc(N, sizeof(int));
        initial_graph(&graph, spins, N, M);
        
        temp[i] = MIN_TEMP + ((MAX_TEMP - MIN_TEMP) / TEMP_STEPS) * i;       
        
        // first thermalization
        mag = calloc(T, sizeof(int));
        energy = malloc(T * sizeof(double));
        algorithm_two(&graph, spins, N, M, mag, energy, T, 1.0/temp[i]);  // TODO additional arg to omit calculationg mag and energy
        free(mag);
        free(energy);
        
        // actual calculations
        mag = calloc(SAMPLE_TIME, sizeof(int));
        energy = malloc(SAMPLE_TIME * sizeof(double));
        algorithm_two(&graph, spins, N, M, mag, energy, SAMPLE_TIME, 1.0/temp[i]);

        // computing results
        tmp = avg_int(mag, SAMPLE_TIME);
        mag_avg[i] = abs(tmp);  // take absolute value of magnetization
        mag_std[i] = std_int(mag, SAMPLE_TIME, tmp);
        
        energy_avg[i] = avg_double(energy, SAMPLE_TIME);
        energy_std[i] = std_double(energy, SAMPLE_TIME, energy_avg[i]);
        
        free(mag);
        free(spins);
        free(energy);
        igraph_destroy(&graph);
        t_end = clock();
        printf("%d Computing energy and magnetization for T=%f finished in %f s.\n",
            i, temp[i], (double)(t_end - t_start) / CLOCKS_PER_SEC);
    }
    
    char file_name1[] = "res/energy_vs_B_2.csv";
    save_to_csv_xystd_double(file_name1, temp, energy_avg, energy_std, TEMP_STEPS);
    
    char file_name2[] = "res/mag_vs_B_2.csv";
    save_to_csv_xystd_double(file_name2, temp, mag_avg, mag_std, TEMP_STEPS);
    return 0;
}


int main(void)
{
    //compare_thermalization();
    phase_diagram_one();
    return 0;
}