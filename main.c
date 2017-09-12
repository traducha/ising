#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <igraph.h>
#include "tools.h"


const int T = 100000;  // thermalization time
const int SAMPLE_TIME = 80000;  // time window to average mag ang energy over (after thermalization)
const int N = 100;  // number of nodes
const int M = 300;  // number of edges <k>=2*M/N, c=M/N
const double GAMMA = 2.0;  // power of the degree in hamiltonian (1.0 is neutral): -{sum_over_i} k_i^gamma
const double ALPHA = 2.0;  // power of Jij (0.0 is neutral)
const double J = 1.0;  // J in hamiltonian
const double h = 0.0;  // h in hamiltonian
const double FI = 0.5; // probability of switching edge instead of spin
const double FI_MIN = 0.5;
const double FI_MAX = 0.5;
const int FI_STEPS = 1;
const double B = 1.0;  // 1/kbT inverse of the temperature
const double MIN_TEMP = 0.1;  // min temperature
const double MAX_TEMP = 150.0;  // max temperature 0.004
const int TEMP_STEPS = 100;  // number of values for temperature


void algorithm_one(igraph_t *graph, int *spins, int nodes, int edges,
                        int *mag, double *energy, int times, double beta)
/* OLD, just for comparsion
 * Metropolis algorithm of Ising model with <times> time steps and saving magnetization and energy.
 * Currently with switching only one end of a node.
 */
{
    igraph_vector_t neigs, edges_vector;
    igraph_vector_init(&neigs, 1);
    igraph_vector_init(&edges_vector, 1);
    igraph_integer_t e_from, e_to, v_to_new;
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
        energy[i] = compute_energy(graph, nodes, &edges_vector, spins, ALPHA, h, GAMMA);
        mag[i] = sum_array_int(spins, nodes);
    }
    
    igraph_vector_destroy(&neigs);
    igraph_vector_destroy(&edges_vector);
    return;
}


void algorithm_two(igraph_t *graph, int *spins, int nodes, int edges, int *mag,
                   double *energy, int times, double beta)
/*
 * Metropolis algorithm of Ising model with <times> time staps and saving magnetization and energy.
 * This version implements switching both ends of an edge.
 * About 0.5% faster than first version.
 */
{
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
            delta = energy_change_alpha_spin(graph, v_index, &neigs, -2 * spins[v_index], // -2 * spin is the possible energy change
                                             spins, ALPHA, 2.0 * edges / nodes);
            delta +=  h * 2 * spins[v_index];
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
        
        delta = energy_change_alpha_egde(graph, (int) e_from, (int) e_to, new_e[0], new_e[1],
                                             spins, ALPHA, 2.0 * edges / nodes);
        delta += energy_change_gamma(graph, (int) e_from, (int) e_to, new_e[0], new_e[1], GAMMA);
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
        energy[i] = compute_energy(graph, nodes, &edges_vector, spins, ALPHA, h, GAMMA);
        mag[i] = sum_array_int(spins, nodes);
    }
    
    igraph_vector_destroy(&neigs);
    igraph_vector_destroy(&edges_vector);
    return;
}


void algorithm_two_complex(igraph_t *graph, int *spins, int nodes, int edges, int *mag, int *mag_abs,
            int *incompatible, double *energy, int *clust_num, int *largest_clust,
            int *largest_degree, int times, double beta, double fi)
/*
 * Metropolis algorithm of Ising model with <times> time staps and saving
 * magnetization (and its abs), energy, largest component, number of components,
 * and number of incompatible links, largest degree...
 * This version implements switching both ends of an edge.
 * About 0.5% faster than first version.
 */
{
    igraph_vector_t neigs, edges_vector, clusters;
    igraph_vector_init(&neigs, 1);
    igraph_vector_init(&edges_vector, 1);
    igraph_vector_init(&clusters, 1);
    igraph_integer_t e_from, e_to;
    igraph_es_t es;
    double r, delta, proper_fi = fi;
    int v_index, e_index, tmp_clust_num, *new_e;
    
    if (proper_fi > 1.0) {proper_fi = FI;}
    
    int i;
    for (i = 0; i < times; i++)
    {
        // select vertex and edge
        r = rand();
        v_index = (int)(((nodes - 1) * r / RAND_MAX) + 0.5);
        r = rand();
        e_index = (int)(((edges - 1) * r / RAND_MAX) + 0.5);

        r = rand();
        if (r / RAND_MAX > proper_fi)
        {
            // spin switching
            igraph_neighbors(graph, &neigs, v_index, IGRAPH_ALL);
            delta = energy_change_alpha_spin(graph, v_index, &neigs, -2 * spins[v_index], // -2 * spin is the possible energy change
                                             spins, ALPHA, 2.0 * edges / nodes);
            delta +=  h * 2 * spins[v_index]; 
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
        }
        else
        {
            // edge rewiring
            igraph_edge(graph, e_index, &e_from, &e_to);
            
            // draw new ends for the edge
            new_e = draw_new_edge(graph, nodes, (int) e_from, (int) e_to);
            
            delta = energy_change_alpha_egde(graph, (int) e_from, (int) e_to, new_e[0], new_e[1],
                                             spins, ALPHA, 2.0 * edges / nodes);
            delta += energy_change_gamma(graph, (int) e_from, (int) e_to, new_e[0], new_e[1], GAMMA);
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
        }
        
        // compute magnetization, internal energy, number and size of clusters, incompatible links
        igraph_get_edgelist(graph, &edges_vector, 0);
        energy[i] = compute_energy(graph, nodes, &edges_vector, spins, ALPHA, h, GAMMA);
        mag[i] = sum_array_int(spins, nodes);
        mag_abs[i] = abs(mag[i]);
        
        igraph_clusters(graph, NULL, &clusters, &tmp_clust_num, IGRAPH_WEAK);
        clust_num[i] = tmp_clust_num;
        largest_clust[i] = (int) igraph_vector_max(&clusters);
        
        incompatible[i] = count_incompatible_links(&edges_vector, spins);
        
        igraph_degree(graph, &neigs, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
        largest_degree[i] = (int) igraph_vector_max(&neigs);
    }
    
    igraph_vector_destroy(&neigs);
    igraph_vector_destroy(&edges_vector);
    igraph_vector_destroy(&clusters);
    return;
}


void algorithm_two_thermalize(igraph_t *graph, int *spins, int nodes, int edges,
                              int times, double beta, double fi)
/*
 * Copy of algorithm_two() without computing and returning any quantities,
 * for the faster thermalization.
 */
{
    igraph_vector_t neigs, edges_vector;
    igraph_vector_init(&neigs, 1);
    igraph_vector_init(&edges_vector, 1);
    igraph_integer_t e_from, e_to;
    igraph_es_t es;
    double r, delta, proper_fi = fi;
    int v_index, e_index, *new_e;
    
    if (proper_fi > 1.0) {proper_fi = FI;}
    
    int i;
    for (i = 0; i < times; i++)
    {
        // select vertex and edge
        r = rand();
        v_index = (int)(((nodes - 1) * r / RAND_MAX) + 0.5);
        r = rand();
        e_index = (int)(((edges - 1) * r / RAND_MAX) + 0.5);

        r = rand();
        if (r / RAND_MAX > proper_fi)
        {
            // spin switching
            igraph_neighbors(graph, &neigs, v_index, IGRAPH_ALL);
            delta = energy_change_alpha_spin(graph, v_index, &neigs, -2 * spins[v_index], // -2 * spin is the possible energy change
                                             spins, ALPHA, 2.0 * edges / nodes);
            delta +=  h * 2 * spins[v_index];
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
        }
        else
        {
            // edge rewiring
            igraph_edge(graph, e_index, &e_from, &e_to);
            
            // draw new ends for the edge
            new_e = draw_new_edge(graph, nodes, (int) e_from, (int) e_to);
            
            delta = energy_change_alpha_egde(graph, (int) e_from, (int) e_to, new_e[0], new_e[1],
                                             spins, ALPHA, 2.0 * edges / nodes);
            delta += energy_change_gamma(graph, (int) e_from, (int) e_to, new_e[0], new_e[1], GAMMA);
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
        }
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


int phase_diagram_two(double fi)
/*
 * Computes magnetization, energy, size of the largest component,
 * and number of components (after thermalization)
 * for different temperatures
 */
{
    clock_t t_start, t_end, loop_start, loop_end;
    int *spins, *mag, *mag_abs, *incompatible, *largest_clust, *clust_num, *largest_degree;
    double *energy, tmp, temp[TEMP_STEPS];
    igraph_t graph;
    
    // results
    double mag_avg[TEMP_STEPS], mag_std[TEMP_STEPS];
    double mag_abs_avg[TEMP_STEPS], mag_abs_std[TEMP_STEPS];
    double energy_avg[TEMP_STEPS], energy_std[TEMP_STEPS];
    double incompatible_avg[TEMP_STEPS], incompatible_std[TEMP_STEPS];
    double largest_clust_avg[TEMP_STEPS], largest_clust_std[TEMP_STEPS];
    double clust_num_avg[TEMP_STEPS], clust_num_std[TEMP_STEPS];
    double largest_degree_avg[TEMP_STEPS], largest_degree_std[TEMP_STEPS];
    
    loop_start = clock();
    int i;
    for (i = 0; i < TEMP_STEPS; i++)
    {
        t_start = clock();
        spins = calloc(N, sizeof(int));
        initial_graph(&graph, spins, N, M);
        
        temp[i] = MIN_TEMP + ((MAX_TEMP - MIN_TEMP) / TEMP_STEPS) * i;       
        
        // first of all thermalization
        algorithm_two_thermalize(&graph, spins, N, M, T, 1.0/temp[i], fi);
        
        // actual calculations
        mag = calloc(SAMPLE_TIME, sizeof(int));
        mag_abs = calloc(SAMPLE_TIME, sizeof(int));
        incompatible = calloc(SAMPLE_TIME, sizeof(int));
        clust_num = calloc(SAMPLE_TIME, sizeof(int));
        largest_clust = calloc(SAMPLE_TIME, sizeof(int));
        largest_degree = calloc(SAMPLE_TIME, sizeof(int));
        energy = malloc(SAMPLE_TIME * sizeof(double));
        algorithm_two_complex(&graph, spins, N, M, mag, mag_abs, incompatible, energy,
                      clust_num, largest_clust, largest_degree, SAMPLE_TIME, 1.0/temp[i], fi);

        // computing results
        tmp = avg_int(mag, SAMPLE_TIME);
        mag_avg[i] = abs(tmp);  // take absolute value of magnetization
        mag_std[i] = std_int(mag, SAMPLE_TIME, tmp);
        
        mag_abs_avg[i] = avg_int(mag_abs, SAMPLE_TIME);
        mag_abs_std[i] = std_int(mag_abs, SAMPLE_TIME, mag_abs_avg[i]);
        
        energy_avg[i] = avg_double(energy, SAMPLE_TIME);
        energy_std[i] = std_double(energy, SAMPLE_TIME, energy_avg[i]);
        
        incompatible_avg[i] = avg_int(incompatible, SAMPLE_TIME);
        incompatible_std[i] = std_int(incompatible, SAMPLE_TIME, incompatible_avg[i]);
        
        largest_clust_avg[i] = avg_int(largest_clust, SAMPLE_TIME);
        largest_clust_std[i] = std_int(largest_clust, SAMPLE_TIME, largest_clust_avg[i]);
        
        clust_num_avg[i] = avg_int(clust_num, SAMPLE_TIME);
        clust_num_std[i] = std_int(clust_num, SAMPLE_TIME, clust_num_avg[i]);
        
        largest_degree_avg[i] = avg_int(largest_degree, SAMPLE_TIME);
        largest_degree_std[i] = std_int(largest_degree, SAMPLE_TIME, largest_degree_avg[i]);
                
        // cleaning
        free(mag);
        free(mag_abs);
        free(incompatible);
        free(clust_num);
        free(largest_clust);
        free(spins);
        free(energy);
        igraph_destroy(&graph);
        t_end = clock();
        printf("%d Computing energy and magnetization for T=%f finished in %f s.\n",
            i, temp[i], (double)(t_end - t_start) / CLOCKS_PER_SEC);
    }
    loop_end = clock();
    
    // saving results
    char *file_name;
    file_name = malloc(100 * sizeof(char));
    
    sprintf(file_name, "res/energy_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, energy_avg, energy_std, TEMP_STEPS);
    
    sprintf(file_name, "res/mag_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, mag_avg, mag_std, TEMP_STEPS);
    
    sprintf(file_name, "res/mag_abs_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, mag_abs_avg, mag_abs_std, TEMP_STEPS);
    
    sprintf(file_name, "res/incompatible_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, incompatible_avg, incompatible_std, TEMP_STEPS);
    
    sprintf(file_name, "res/largest_clust_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, largest_clust_avg, largest_clust_std, TEMP_STEPS);
    
    sprintf(file_name, "res/clust_num_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, clust_num_avg, clust_num_std, TEMP_STEPS);
    
    sprintf(file_name, "res/largest_degree_vs_B_N%d_L%d_J%f_h%f_FI%f_GA%f_AL%f.csv", N, M, J, h, fi, GAMMA, ALPHA);
    save_to_csv_xystd_double(file_name, temp, largest_degree_avg, largest_degree_std, TEMP_STEPS);
    
    free(file_name);    
    printf("Computing energy and magnetization for %d times finished in %f s.\n",
            TEMP_STEPS, (double)(loop_end - loop_start) / CLOCKS_PER_SEC);
    return 0;
}


void phase_diagram_two_over_fi(void)
{
    double fi;
    
    int i;
    for (i = 0; i < FI_STEPS; i++)
    {
        fi = (double) FI_MIN + i * (FI_MAX - FI_MIN) / (FI_STEPS - 1);
        printf("%d/%d Starting computation for FI=%f\n", i, FI_STEPS, fi);
        phase_diagram_two(fi);
    }
}


void tests()
{
    char *test;
    test = malloc(100 * sizeof(char));
    sprintf(test, "%s %s %s %d", "We", "are", "in", 2012);
    printf("WYNIK: %s \n", test);
    sprintf(test, "%s %s %s %d", "We", "are", "in", 2013);
    printf("WYNIK: %s \n", test);
    free(test);
    return;
}


int main(void)
{
    srand(time(NULL));
    
    //compare_thermalization();
    phase_diagram_two(FI);
    //phase_diagram_two_over_fi();
    //tests();
    return 0;
}
