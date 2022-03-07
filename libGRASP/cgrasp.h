/*
 * cgrasp.h
 *
 *  Created on: Feb 25, 2010
 *      Author: Ricardo Martins de Abreu Silva
 */

#ifndef CGRASP_H_
#define CGRASP_H_

#define DEBUG_LEVEL0_ 0
#define DEBUG_LEVEL1_ 1
#define DEBUG_LEVEL2_ 2
#define DEBUG_LEVEL3_ 3
#define DEBUG_LEVEL4_ 4

#include <python2.7/Python.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ltdl.h>
#include <vector>
#include <string>

#include "simclist.h"
#include "mt19937ar.h"
#include "Funcao.h"

#define _infinity 1e+30
#define PI 3.14159265

#define ELITE_SIZE 10   // Tamaho máximo do conjunto elite


typedef struct  
{
        double *x;
        double cost;
} Solution;

typedef struct 
{
        Solution *solutions;   // lista de soluções da elite
        double  inf_threshold;  // limiar inferior da lista de custo (Pior custo da elite)
        int current_size;       // Tamanho atual da elite
        double *_template;       
} EliteSet;

typedef struct 
{
        int label;
        int cluster_size;
        
        double position;
        double *cluster;
} Centroid;

typedef enum  {
        SAME_TIME,
        AGAINST_CGRASP,
        AGAINST_OTHER_HEURISTICS,
        AGAINST_DCGRASP,
        WITH_HART_STOP_RULE
} Experiment;

enum ExperimentType { BY_ITERATION, BEFORE_AFTER, BY_CFOs, BY_GAP };

typedef struct _results_
{
        double  best_fo;
        double* best_solution;
        int number_of_cfos;

        double best_fo_before_mining;
        double* best_solution_before_mining;
        int number_of_cfos_before_mining;

} Results;

typedef struct _algorithmBehaviour_
{
        
        std::vector<double*> solutions_by_iteration;
        std::vector<double*> generated_solutions_by_iteration;

        std::vector<double> costs_by_iteration;
        std::vector<double> generated_costs_by_iteration;
        
        std::vector<double*> solutions_history;
        std::vector<double*> generated_solutions_history;

        std::vector<double> costs_history;
        std::vector<double> generated_costs_history;

        std::vector<double**> elite_by_iteration;
        std::vector<double*> elite_costs_by_iteration;
        
        std::vector<double**> elite_history;
        std::vector<double*> elite_costs_history;

        double*  elite_costs;
        double** elite_dataset;

        std::vector<double*> deviations_in_elite_by_dimension;
        std::vector<double> elite_avg_costs;
        std::vector<double> elite_best_costs;
        std::vector<double> elite_worst_costs;

} AlgorithmBehaviour;

double cgrasp( const char * funcCode,const char * mining_strategy, int isContinuousMining, double dmStartMoment, double patternPercentUsed, int eliteSize, int dimension, double *lower_bounds, double *upper_bounds, 
        Funcao *func, double hs, double he, double plo, int number_of_iterations, int max_functions_calls, unsigned long seed, bool &success, 
        double standardDeviation = 0);
        
double cgrasp_cec(int n, double ep, int seed, double h_s, double h_e, double ro, int max_points, int func_num, void (*func)(double *, double *, int, int, int));

#endif /*CGRASP_H_*/
