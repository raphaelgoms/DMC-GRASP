/*
 * cgrasp.c
 *
 *  Created on: Oct 25, 2010
 *      Author: Ricardo Martins de Abreu Silva
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "Util.h"

#include "cgrasp.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

//#include "module.h"

#if HAVE_LTDL_H
#include <ltdl.h>
#endif

#include <sys/resource.h>
#include <algorithm>
#include<numeric>
#include "xmeans.h"




using namespace std;

#include <random>
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(1); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0, 1);

int debugLevel = DEBUG_LEVEL0_;

int MAX_ITERATION;
int MAX_nCFOs;

//typedef double entrypoint ( double * arg);
//entrypoint * run;

float res; // retorno do cgrasp (melhor valor de f.o. encontrado)

double std_threshold;
int useDataMining = 0;
int g_elite_size;
int extraiu = 0;
int numDivisions = 0;
int CURRENTSIZE = 0;
double * deviations;

double currentStep;

int block = 0;

FILE *solutionsFile;
FILE *solutionsEliteFile;

FILE *foBeforePatternExtrationFile;
FILE *foAfterPatternExtrationFile;

FILE *resFile;

double *_template;

PyObject *pdict, *pName, *pModule, *pFunc, *pList, *pValue, *pArgs, *pvar, *pl, *pu;

typedef struct
{
	int i;		// i-th coordinate
	double z_i; // value for the i-th coordinate that minimizes the obj. function with the other n-1 coordinates of x held at their current values
	double g_i; // objective function value for x with the i-th coordinate set to z_i
} item;

size_t mymeter(const void *el)
{
	/* every element has the constant size of a item structure */
	return sizeof(item);
}

/*
 * compare items by coordinates
 *
 * this function compares two coordinates from two items:
 * <0: a greater than b
 * 0: a equivalent to b
 * >0: b greater than a
 */
int mycomparator(const void *a, const void *b)
{
	/* compare items */
	const item *A = (item *)a;
	const item *B = (item *)b;
	return (A->i < B->i) - (A->i > B->i);
}

/* return "match" when the contact ID matches the key */
int seeker(const void *el, const void *key)
{
	/* let's assume el and key being always != NULL */
	const item *cont = (item *)el;

	if (cont->i == *(int *)key)
		return 1;
	return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double urand()
{
	// return drand48();
	//    return (double(rand())/RAND_MAX);
	return genrand_real1() / RAND_MAX;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int urand_ind(int max_ind)
{
	//return rand()%max_ind;
	//return (int)(genrand_real1() * max_ind);

	//unsigned long r = genrand_int32();
	return floor (dis(gen)* max_ind);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int * dimSortedByDeviations;

void sortByDeviations(int size) 
{
	//printf("sort by deviation\n");
	int i, j;
	dimSortedByDeviations = (int*) malloc(size * sizeof(int));
	for (i = 0; i < size; i++) {
		dimSortedByDeviations[i]=i;
	}

	for (j = 1; j < size; j++)
	{
		int current = dimSortedByDeviations[j];
		double key = deviations[j];

		int i = j - 1;
		while (i > -1 && deviations[i] > key)
		{
			dimSortedByDeviations[i + 1] = dimSortedByDeviations[i];
			deviations[i + 1] = deviations[i];
			i--;
		}

		dimSortedByDeviations[i + 1] = dimSortedByDeviations[i];
		deviations[i + 1] = key;
	}
}

int selectByRank(int n)
{
	double r = genrand_real1();
	//printf("%f\n", r);

	int rank = n;
	int selected = 0;

	double full_sum = (n *(n+1))/2;
	double cum_prob = (double)rank/full_sum;
	while (r > cum_prob) {
		selected++;
		rank--;
		cum_prob += (double)rank/full_sum;
	}
	
	return selected;
}

int * getRandSubSet(int n, int k, int use_rank) {
	int i, *set, *aux;

	if (use_rank)
		sortByDeviations(n);
	
	set = (int*) malloc (k*sizeof(int));

	if (use_rank) {
		aux = dimSortedByDeviations;
	} else {
		aux = (int*) malloc (n*sizeof(int));
		for (i=0; i<n; i++) {
			aux[i] = i;
		}
	}

	int m = n-1, r;
	for (i=0; i<k; i++) {
		if (use_rank) {
			r = selectByRank(m);
		} else {
			r = urand_ind(m);
		}

		set[i] = aux[r]; 
		aux[r] = aux[m];
		m--;
	}

	return set;
}



void UnifRand(int n, double *l, double *u, double *x, int *fixed, double *__template_, double percent)
{
	int i;
	//double percent = 1.0;
	if (extraiu && useDataMining)
	{
		//printf("USING PATTERN!\n");
		if (percent == 1.0) {
			//printf("COMPLETE!\n");	
			for (i = 0; i < n; i++)
			{
				if(fixed[i])
					(*(x + i)) = __template_[i];
				else
					(*(x + i)) = (*(l + i)) + genrand_real1() * ((*(u + i)) - (*(l + i)));
			}
		} else {
			
			int sss = percent * n;
			int * sub = getRandSubSet(n, sss, 0);
			for (i = 0; i < n; i++)
			{
				(*(x + i)) = (*(l + i)) + genrand_real1() * ((*(u + i)) - (*(l + i)));
			}

			for (i = 0; i < sss; i++)
			{
				if(fixed[sub[i]])
					(*(x + sub[i])) = __template_[sub[i]];
			}
		}		
	}
	else
	{
		

		for (i = 0; i < n; i++)	{
			x[i] = l[i]  + dis(gen) * (u[i] - l[i]);
		}

		//Util::printX(x, n); cout << endl;
	}
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float timer()
{
	struct rusage r;

	getrusage(0, &r);
	return (float)(r.ru_utime.tv_sec + r.ru_utime.tv_usec / (float)1000000);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void print_solution(int n, double *x, FILE *outFile)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (i + 1 == n)
			fprintf(outFile, "%f", *(x + i));
		else
			fprintf(outFile, "%f, ", *(x + i));
	}
	fprintf(outFile, "\n");
}

void print_solution_in_screen(int n, double *x, FILE *outFile)
{
	int i;
	for (i = 0; i < n; i++)
	{
		printf("%f ", *(x + i));
	}
	printf("\nsolution: ");
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void copy_double_vector(double *dst, double *src, int sz)
{
	int i;
	for (i = 0; i < sz; i++)
	{
		*(dst + i) = *(src + i);
	}
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void define_precision(double *number, int precision)
{
	*(number) = floor(pow(10, precision) * (*(number))) / pow(10, precision);
}

void copy(double *s, double *other, int N)
{
	int i = 0;
	for (i = 0; i < N; i++)
	{
		s[i] = other[i];
	}
}

int compareTwoSolutions(Solution s1, Solution s2)
{
	return (s1.cost > s2.cost);
}

void sortSolutions(Solution *solutions, int n)
{
	int j;
	for (j = 1; j < n; j++)
	{
		Solution current = solutions[j];
		double key = current.cost;

		int i = j - 1;
		while (i > -1 && solutions[i].cost > key)
		{
			solutions[i + 1] = solutions[i];
			i--;
		}
		solutions[i + 1] = current;
	}
}

bool updateEliteSet(double *solution, double cost, EliteSet *elite, int N)
{
	int i, j;
	int position = -1;
	bool wasUpdated = false;

	if (CURRENTSIZE < g_elite_size)
	{
		Solution *s = &elite->solutions[CURRENTSIZE];
		s->x = (double *)malloc(sizeof(double) * N);
		copy_double_vector(s->x, solution, N);
		s->cost = cost;

		//printf("cost: %lf\n", cost);

		CURRENTSIZE = CURRENTSIZE + 1;
		wasUpdated = true;
	}
	else if (cost < elite->solutions[CURRENTSIZE - 1].cost)
	{
		Solution *s = &elite->solutions[CURRENTSIZE - 1];
		//s->x = (double*)malloc(sizeof(double)*N);
		copy_double_vector(s->x, solution, N);
		s->cost = cost;
		wasUpdated = true;
	}

	//printf("cost: %d\n", CURRENTSIZE);
	if (wasUpdated)
		sortSolutions(elite->solutions, CURRENTSIZE);

	return wasUpdated;
}

double calculateAverageCostElite(EliteSet *elite, int elite_size)
{
	double acc = 0.0;
	for (int j = 0; j < elite_size; j++) {
		acc += elite->solutions[j].cost;
	}

	return acc / elite_size;
}

double getWorstCostElite(EliteSet *elite, int elitesize)
{
	double worst = -INFINITY;
	for (int j = 0; j < elitesize; j++) {
		if (worst < elite->solutions[j].cost)
			worst = elite->solutions[j].cost;
	}

	return worst;
}

double getBestCostElite(EliteSet *elite, int elitesize)
{
	double best = INFINITY;
	for (int j = 0; j < elitesize; j++) {
		if (best > elite->solutions[j].cost)
			best = elite->solutions[j].cost;
	}

	return best;
}


double* calculateAverageInElite(EliteSet *elite, int elite_size, int number_of_variables)
{
	double* averages = new double[number_of_variables];

	for (int i = 0; i < number_of_variables; i++) {
		double acc = 0.0;
		for (int j = 0; j < elite_size; j++) {
			acc += elite->solutions[j].x[i];
		}

		averages[i] = acc / elite_size;
	}

	return averages;
}

double getNormalizedValue(double value, double min, double max)
{
	return (value - min) / (max - min);
}

double* calculateDeviationsInElite(EliteSet *elite, int elite_size,
 int number_of_variables, double *min, double *max, double* averages = NULL)
{
	if (!averages) {
		averages = calculateAverageInElite(elite, elite_size, number_of_variables);
	}
	
	double *standard_deviantions = new double[number_of_variables];
	for (int i = 0; i < number_of_variables; i++)	{
		double acc = 0.0;
		double normalized_mean = getNormalizedValue(averages[i], min[i], max[i]);
		for (int j = 0; j < elite_size; j++)
		{ 
			Solution s = elite->solutions[j];
			double normalized_value = getNormalizedValue(s.x[i], min[i], max[i]);

			acc += (normalized_mean - normalized_mean) * (normalized_mean - normalized_mean);
		}

		standard_deviantions[i] = sqrt(acc / elite_size);
	}

	return standard_deviantions;
}

// retorna conjuntos de posições fixas
void FindPatterns(EliteSet *elite, int problem_dimension, int *toFix, double *l, double *u, double h)
{
	double StdDevThresold = std_threshold;
	int elite_size = CURRENTSIZE;
	double* pattern = calculateAverageInElite(elite, elite_size, problem_dimension);
	deviations = calculateDeviationsInElite(elite, elite_size, problem_dimension, l, u,  pattern);

	// construir o conjunto de posições de fixas
	for (int i = 0; i < problem_dimension; i++) {
		elite->_template[i] = pattern[i];

		if (StdDevThresold > 0) {
			if (StdDevThresold > deviations[i]) {
				toFix[i] = 1;
			}
			else {
				toFix[i] = 0;
			}
		}
		else {
			toFix[i] = 1;
		}
	}

	delete pattern;
}

void FindPatternsWithXmeans(EliteSet *elite, int problem_dimension, int *toFix, double *l, double *u, double h, bool clustering_before)
{
	Xmeans xmeans;
	vector<vector<double>> points;
	vector<double> lower_bounds;
	vector<double> upper_bounds;
	int elite_size = CURRENTSIZE;
	double StdDevThresold = std_threshold;

	for (int i = 0; i < elite_size; i++)
	{
		vector<double> point;
		for (int j = 0; j < problem_dimension; j++)
		{
			point.push_back(elite->solutions[i].x[j]);
		}
		points.push_back(point);
		
	}

    for (int i = 0; i < problem_dimension; i++)
	{
		lower_bounds.push_back(l[i]);
		upper_bounds.push_back(u[i]);
	}
	
	map<int, double> pattern;
	if (clustering_before)
		pattern = xmeans.extractPatternAfterClustering(points, lower_bounds, upper_bounds);
	else
		pattern = xmeans.extractPatternXMeans(points, lower_bounds, upper_bounds);

	for (int i = 0; i < problem_dimension; i++) {
		toFix[i] = 0;
	}

	double* average = calculateAverageInElite(elite, elite_size, problem_dimension);
	deviations = calculateDeviationsInElite(elite, elite_size, problem_dimension, l, u, average);
	std::map<int, double>::iterator it = pattern.begin();
    while (it != pattern.end())
    {
		elite->_template[it->first] = it->second;

		if (StdDevThresold == 0 || deviations[it->first] < StdDevThresold) {
			// selecionar somente variaveis menos dispersas
			toFix[it->first] = 1; 
		}

		it++;
	}
}

void ChooseRandomPattern(EliteSet *elite, int problem_dimension, int *toFix, double *l, double *u, double h, bool clustering_before)
{
	Xmeans xmeans;
	vector<vector<double>> points;
	vector<double> lower_bounds;
	vector<double> upper_bounds;
	int elite_size = CURRENTSIZE;
	double StdDevThresold = std_threshold;

	for (int i = 0; i < elite_size; i++)
	{
		vector<double> point;
		for (int j = 0; j < problem_dimension; j++)
		{
			point.push_back(elite->solutions[i].x[j]);
		}
		points.push_back(point);
		
	}

    for (int i = 0; i < problem_dimension; i++)
	{
		lower_bounds.push_back(l[i]);
		upper_bounds.push_back(u[i]);
	}

	vector<map<int, double>> patterns = xmeans.extractPatterns(points, lower_bounds, upper_bounds);

	for (int i = 0; i < problem_dimension; i++) {
		toFix[i] = 0;
	}

	double* average = calculateAverageInElite(elite, elite_size, problem_dimension);
	deviations = calculateDeviationsInElite(elite, elite_size, problem_dimension, l, u, average);
	
	int r = urand_ind(patterns.size());
	//cout << " -  r: " << r << endl;
	map<int, double> pattern = patterns[r];
	std::map<int, double>::iterator it = pattern.begin();
    while (it != pattern.end())
    {
		elite->_template[it->first] = it->second;

		if (StdDevThresold == 0 || deviations[it->first] < StdDevThresold) {
			// selecionar somente variaveis menos dispersas
			toFix[it->first] = 1; 
		}

		it++;
	}	
}


vector<map<int, double>> ExtractPatterns(EliteSet *elite, int problem_dimension, int *toFix, double *l, double *u, double h, bool clustering_before, double &max_radius)
{
	Xmeans xmeans;
	vector<vector<double>> points;
	vector<double> lower_bounds;
	vector<double> upper_bounds;
	int elite_size = CURRENTSIZE;
	double StdDevThresold = std_threshold;

	for (int i = 0; i < elite_size; i++)
	{
		vector<double> point;
		for (int j = 0; j < problem_dimension; j++)
		{
			point.push_back(elite->solutions[i].x[j]);
		}
		points.push_back(point);
		
	}

    for (int i = 0; i < problem_dimension; i++)
	{
		lower_bounds.push_back(l[i]);
		upper_bounds.push_back(u[i]);
	}

	vector<map<int, double>> patterns = xmeans.extractPatterns(points, lower_bounds, upper_bounds);
	max_radius = xmeans.getMaxRadius();
	return patterns;
}


void FindSubSetPatterns(EliteSet *elite, int N, int *toFix, double *l, double *u, double h, double threshold)
{
	int histogram[N];
	int *Fixed;
	double StdDevThresold = std_threshold;
	int i, j;

	int elite_size = CURRENTSIZE;

	// calcular m�dia;
	double *media = (double *)malloc(N * sizeof(double));
	double *pattern = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		pattern[i] = 0.0;
		media[i] = 0.0;
		for (j = 0; j < elite_size; j++)
		{
			Solution s = elite->solutions[j];
			pattern[i] += s.x[i];
		}

		pattern[i] = pattern[i] / elite_size;
	}

	// calcular desvio padr�o
	double *stdDeviation = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		stdDeviation[i] = 0.0;
		for (j = 0; j < elite_size; j++)
		{ 
			Solution s = elite->solutions[j];

			double normalized_value = (s.x[i] - l[i]) / (u[i] - l[i]);
			double normalized_mean = (pattern[i] - l[i]) / (u[i] - l[i]);
			
			stdDeviation[i] += (normalized_value - normalized_mean) * (normalized_value - normalized_mean);
			//stdDeviation[i] += (s->x[i] - pattern[i]) * (s->x[i] - pattern[i]);
		}

		stdDeviation[i] = sqrt(stdDeviation[i] / elite_size);
	}

	// construir o conjunto de posi��es de fixas
	int k = 0;
	for (i = 0; i < N; i++)
	{	
		double normalized_mean = (pattern[i] - l[i]) / (u[i] - l[i]);
		//printf(" %lf ", stdDeviation[i]/normalized_mean);
		if (stdDeviation[i]/normalized_mean < threshold)
		{
			toFix[i] = 1;
			elite->_template[i] = pattern[i];
		} 
		else
		{
			toFix[i]=0;
		}
	}
	//printf("\n");
}


void FindPatternsThr(EliteSet *elite, int N, int *toFix, double *l, double *u, double h)
{
	int histogram[N];
	int *Fixed;
	double StdDevThresold = std_threshold;
	int i, j;

	int elite_size = CURRENTSIZE;

	// calcular média;
	double *media = (double *)malloc(N * sizeof(double));
	double *pattern = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		media[i] = 0.0;
		pattern[i] = 0.0;
		for (j = 0; j < elite_size; j++)
		{
			Solution s = elite->solutions[j];
			double normalized = (s.x[i] - l[i]) / (u[i] - l[i]);
			media[i] += normalized;
			pattern[i] += s.x[i];
		}

		media[i] = media[i] / elite_size;
		pattern[i] /= elite_size;
	}

	// calcular desvio padrão
	double *stdDeviation = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		stdDeviation[i] = 0.0;
		for (j = 0; j < elite_size; j++)
		{ 
			Solution s = elite->solutions[j];
			double normalized = (s.x[i] - l[i]) / (u[i] - l[i]);
			stdDeviation[i] += (normalized - media[i]) * (normalized - media[i]);
		}

		stdDeviation[i] = sqrt(stdDeviation[i] / elite_size);
		printf(" %.9f ", stdDeviation[i]);
	}
	printf("\n");

	// definir posições que serão fixas
	int qtdToFix = 0;
	for (i = 0; i < N; i++)
	{
		if (stdDeviation[i] < StdDevThresold)
		{
			toFix[i] = 1;
		}
		else
		{
			toFix[i] = 0;
		}
	}

	// constroir o conjunto de posições de fixas
	int k = 0;
	for (i = 0; i < N; i++)
	{
		if (toFix[i])
		{
			elite->_template[i] = pattern[i];
		}
		else
		{
			elite->_template[i] = l[i] + genrand_real1() * (u[i] - l[i]);
		}
	}
}

void FindPatternsPeak(EliteSet *elite, int N, int *toFix, double *l, double *u, double h)
{
	int *histogram;
	int *Fixed;

	// double StdDevThresold = 0.19;
	// double StdDevThresold = 0.4; // ros

	double StdDevThresold = 0.1; // ack

	int i, j;

	double *media = (double *)malloc(N * sizeof(double));
	int elite_size = CURRENTSIZE;

	// calcular média;
	double *peaks = (double *)malloc(N * sizeof(double));

	for (i = 0; i < N; i++)
	{
		int nbins = (u[i] - l[i]) / h;
		histogram = (int *)malloc(nbins * sizeof(int));

		for (j = 0; j < nbins; j++)
			histogram[j] = 0;

		int peak = -1;
		int maior = 0;
		media[i] = 0.0;
		for (j = 0; j < elite_size; j++)
		{
			Solution s = elite->solutions[j];
			media[i] += s.x[i] / elite_size;
			int pos = (s.x[i] - l[i]) / h;
			histogram[pos]++;

			if (histogram[pos] > maior)
			{
				maior = histogram[pos];
				peak = pos;
			}
		}

		peaks[i] = peak * h + l[i];
		//printf("Qtde x%d: %d\n", i, maior);
	}

	// calcular desvio padrão
	double *stdDeviation = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		stdDeviation[i] = 0.0;
		for (j = 0; j < elite_size; j++)
		{
			Solution s = elite->solutions[j];
			stdDeviation[i] += (s.x[i] - media[i]) * (s.x[i] - media[i]);
		}

		stdDeviation[i] = sqrt(stdDeviation[i] / N);
	}

	// definir posições que serão fixas
	int qtdToFix = 0;
	for (i = 0; i < N; i++)
	{
		//printf(" %f ", stdDeviation[i]/media[i]);
		// 0.0003
		// 0.1
		toFix[i] = 1;
	}

	// constroir o conjunto de posições de fixas
	int k = 0;
	for (i = 0; i < N; i++)
	{
		//if (toFix[i])
		//{
		//        printf(" %d ", i);
		elite->_template[i] = peaks[i];
		//}
		//else
		//{
		//        elite->_template[i] = -1;
		//}
	}
}





double average(double *data, int n)
{
	if (n == 0)
		return 0.0;

	double average = 0.0;
	int i;
	for (i = 0; i < n; i++)
	{
		average += data[i];
	}

	average = average / n;

	return average;
}

double standard_deviation(double *data, int n, double media)
{
	if (n == 0)
		return 0.0;

	double soma = 0.0;
	int i;
	for (i = 0; i < n; i++)
	{
		soma += (data[i] - media) * (data[i] - media);
	}

	return sqrt(soma / n);
}

double modulo(double value)
{
	if (value < 0)
		return -1 * value;
	return value;
}

void line_search(int n, double h, double *l, double *u,
				  double *x, int i, double *z_i, double *g_i, float t2, double *x_s, double f,
				  double opt, FILE *outFile, int *count_calc_cost, int *fecount, int stp, int fes, Funcao* func)
{
	//cout << "line_search" << endl;
	int j, start, end, z;
	double x_i, c, s_domain, e_domain;

	double *x_ = (double *) malloc (n * sizeof(double));
	copy_double_vector(x_, x, n);

	x_i = *(x_ + i);
	s_domain = *(l + i);
	e_domain = *(u + i);

	*g_i = func->calc(x_);
	*z_i = x_i;

	// int nSteps = (int) floor ((e_domain - s_domain) / h);
	// double xCurr = s_domain;

	// for (j=0; j<nSteps; j++) {
	// 	*(xAux + i) = s_domain + j*h;
	// 	c = func->calc(xAux);

	// 	if (c < *g_i) {
	// 		*g_i = c;
	// 		*z_i = *(xAux + i);
	// 	}

	// 	xCurr+=h;
	// }

	start = (int)floor((x_i - s_domain) / h);

	for (j = -start; j < 0; j++)
	{
		*(x_ + i) = x_i + j * h;
		c = func->calc(x_);

		if (c < *g_i)
		{
			*g_i = c;
			*z_i = *(x_ + i);
		}
	}

	end = (int)floor((e_domain - x_i) / h);
	for (j = 1; j <= end; j++)
	{
		*(x_ + i) = x_i + j * h;
		c = func->calc(x_);

		if (c < *g_i)
		{
			*g_i = c;
			*z_i = *(x_ + i);
		}
	}
	
	//*(x + i) = *z_i;
	
	free(x_);
}




void construct_greedy_randomized(int n, double ro, double h, double *l, double *u,
								 double *x, int *imprc, float t2, double *x_s, double f, double opt, FILE *outFile, 
								 int *count_calc_cost, int *fecount, int stp, int fes, int *fixed, Funcao *func)
{
	//cout << "construct_greedy_randomized " << h << endl;
	std::vector<int> unfixed (n);
	std::iota(unfixed.begin(), unfixed.end(), 0);
	
	double alpha = dis(gen);
	bool reuse = false;
	double g_min, g_max;

	double *z = (double*) malloc(n *sizeof(double));
	double *g =  (double*) malloc(n *sizeof(double));

	// for (int i = 0; i < n; i++) z[i] = INFINITY;
	// for (int i = 0; i < n; i++) g[i] = INFINITY;

	std::vector<int> rcl;

	double cost = func->calc(x);

	//Util::printX(x, n);
	//cout << ", cost: " << cost << ", h = " << h << endl;

	*(imprc) = 0;
	//cout << "construção - unfixed size: " << unfixed.size() << endl;
	while(unfixed.size()) {
		//cout << "construção - while" << endl;
		g_min = INFINITY;
		g_max = -INFINITY;

		// cout << "unfixed: " ;
		// for (int i : unfixed) {
		// 	cout << i << " ";
		// } cout << endl;

		for (int i : unfixed)
		{
			if (!reuse) {
				line_search(n, h, l, u, x, i, &z[i], &g[i], t2, x_s, f, opt, outFile, count_calc_cost, fecount, stp, fes, func);
			}
			//cout << " z_" << i << " " << z[i] << " g_" << i << " " << g[i] << endl;

			if (g_min > g[i]) g_min = g[i];
			if (g_max < g[i]) g_max = g[i];
		}

		//cout << "construção - g_min: " << g_min << " g_max: " << g_max << endl;
		
		rcl.clear();
		double threshold = g_min + alpha * (g_max - g_min);
		for (int i : unfixed)
		{
			if (g[i] <= threshold){
				rcl.push_back(i);
			}
		}
		
		// cout << "rcl: ";
		// for (int i : rcl) {
		// 	cout << i << " ";
		// } cout << endl;
		int j = urand_ind(rcl.size());
		//cout << "construção - rcl size: " << rcl.size() << " j: " << j << endl;
		int selected = rcl[j];

		// cout << "j = " << j << " selected = " << selected << endl; 

		// cout << "x: ";
		// Util::printX(x, n); cout << endl;

		// cout << "z: ";
		// Util::printX(z, n); cout << endl;

		//cout << x[selected] << " - " << z[selected];
		if (Util::equals(x[selected],z[selected])) {
			reuse = true;
		} else {
			x[selected] = z[selected];
			reuse = false;
			*(imprc) = 1;
		}

		unfixed.erase(std::find(unfixed.begin(), unfixed.end(), selected));	
		//cout << unfixed.size() << " ";	
	}
	//cout << endl;

	// double cost_new = func->calc(x);
	// if	(cost_new < cost) {
	// 	*(imprc) = 1;
	// }

	//cout << "fim construct_greedy_randomized" << endl;

	delete[] z;
	delete[] g;
}

double calcNorma(double *x, double *xGrid, int n)
{
	double aux = 0.0, norma;

	int i;
	for (i = 0; i < n; i++){
		aux += pow((xGrid[i]-x[i]), 2);
	}

	norma = sqrt(aux);
	return norma;
}

void randSelectElementBh(double *xCurr, double *xBestAux, double *l, double *u, int n, double h){
	int auxI, i;
	double gridI; 
	double numPointsPos, numPointsNeg;
	double numPoints = 0.0;
	double r, normaVector;
	double *xGrid = (double*) malloc(n *sizeof(double));
	double *bhSelected = (double*) malloc(n *sizeof(double));

	double aux = 0.0, distancia = 0.0;
	double sum = 0.0;
	//printf("\trandSelectElementBh [1]\n");
	// Escolhe aleatoriamente um ponto do Grid de tamanho de passo h.
	do {
		for (i = 0; i < n; i++){
			// Calcula o numero de pontos naquela direcao com tamanho de passo h.
			numPointsPos = floor((u[i] - xBestAux[i])/h);
			numPointsNeg = floor((xBestAux[i] - l[i])/h);
					
			numPoints = numPointsPos + numPointsNeg;

			//numPoints = floor((u[i]-l[i])/h);
			// Escolhe aleatoriamente um dos indices do Grid na direcao i.
			//r = dRand();
			
			r = genrand_real1();
			auxI = (int)(r * (double)numPoints);
			sum += auxI;	
			gridI = ((double)auxI - numPointsNeg); 
						
			// Calcula a posicao da direcao no eixo i.
			xGrid[i] = xBestAux[i] + ((double)gridI*h);
		}
	} while (sum == 0.0);
		//printf("sum; %lf\n", sum);


	//printf("\trandSelectElementBh [2]\n");
	// Calcula a projecao do vetor de x a xGrid.
	normaVector = calcNorma(xBestAux, xGrid, n); 
	// if (normaVector == 0)
	//  	normaVector = 1.0;
	for (i = 0; i < n; i++){
		bhSelected[i] = xBestAux[i] + ((h* (xGrid[i]- xBestAux[i]))/normaVector);
	}

	
	copy_double_vector(xCurr, bhSelected, n);			
	free(xGrid);
	free(bhSelected);

	//printf("\trandSelectElementBh [3]\n");	
 }


void local_search_1(int n, double ro, double h, double *l, double *u, double max_points,
					double *x, int *impr, float t2, double *x_s, double f, double opt, 
					FILE *outFile, int *count_calc_cost, int *fecount, int stp, int fes, Funcao * func)
{
	int ok, i, num_points_examined;
	double f_str, c, norm, s_domain = l[0], e_domain = u[0];

	PyObject *pValue, *pArgs;

	f_str = func->calc(x);

	double *xCurr= (double *)malloc(n * sizeof(double));
	double *xBest = (double *)malloc(n * sizeof(double));

	copy_double_vector(xBest, x, n);
	
	*(impr) = 0;

	s_domain = l[0];
	e_domain = u[0];

	int num_grid_points =(int) pow(floor((e_domain - s_domain)/h),n);
	int max_num_pts = floor(ro*num_grid_points);

	max_points = 1000;
	if ( max_num_pts > max_points ) {
			max_num_pts = max_points;
	}

	//printf("LOCAL SEARCH\n");
	num_points_examined = 0; 

	while (num_points_examined <= max_num_pts) {
		num_points_examined++;
		randSelectElementBh(xCurr, xBest, l, u, n, h);
		double fXTmp = func->calc(xCurr);
		if (Util::feasible(x, l, u, n) && fXTmp < f_str) {
			copy_double_vector(xBest, xCurr, n);
			f_str = fXTmp;
			*(impr) = 1;
			num_points_examined = 0;
		}
	}

	copy_double_vector(x, xBest, n);
	free(xBest);
	free(xCurr);
}


void local_search(int n, double ro, double h, double *l, double *u, double max_points,
					double *x, int *impr, float t2, double *x_s, double f, double opt, 
					FILE *outFile, int *count_calc_cost, int *fecount, int stp, int fes, Funcao *func)
{
	//cout << "local_search" << endl;
	int ok, i, num_points_examined;
	long double f_str, c, norm, s_domain = l[0], e_domain = u[0];

	double *y = (double *)malloc(n * sizeof(double));
	double *z = (double *)malloc(n * sizeof(double));
	for (i = 0; i < n; i++) z[i] = 0.0;
	
	*(impr) = 0;
	f_str = func->calc(x);
	

	long long unsigned num_grid_points = pow(ceil((e_domain - s_domain)/h),n);
	//cout << "e_domain " << e_domain << " s_domain " << s_domain << " h " << h << " n " << n << " num_grid_points " << num_grid_points << endl; 

	int max_num_pts = floor(ro*num_grid_points);

	//max_points = 100000;
	max_points = 1.0e+3;
	if ( max_num_pts > max_points ) {
			max_num_pts = max_points;
	}

	num_points_examined = 0;
	
	//cout << "ro: " << ro << endl;
	//cout << "max_num_pts: " << max_num_pts << endl;
	while (num_points_examined <= max_num_pts)
	{
		//cout << "num_points_examined: " << num_points_examined << endl;
		ok = 0;
		do
		{
			norm = 0.0;
			for (i = 0; i < n; i++)
			{

				s_domain = (*(l + i));
				e_domain = (*(u + i));
				if (dis(gen) <= 0.5)
				{
					*(z + i) = (double)urand_ind((int)ceil((e_domain - (*(x + i))) / h));
				}
				else
				{
					*(z+i) = (double) -1.0*urand_ind( (int)ceil( ((*(x+i))-s_domain)/h) );
				}
				if ((*(z + i)) != 0)
				{
					ok = 1;
				}
				norm = norm + pow((*(z + i)) * h, 2);
			}

		} while (!ok);

		norm = sqrt(norm);

		for (i = 0; i < n; i++)
		{
			*(y + i) = *(x + i) + (1.0 / norm) * h * (*(z + i)) * h;
			
		}

		c = func->calc(y);

		if (f_str > c)
		{
			copy_double_vector(x, y, n);
			f_str = c;
			num_points_examined = 0;
			*(impr) = 1;

			//cout << "f_str: " << f_str << endl;
		}

		num_points_examined++;
	}

	//f_str = func->calc(x);

	free(y);
	free(z);
}

AlgorithmBehaviour createAlgorithmBehaviour(int problem_dimension, int num_of_iterations, int elite_size)
{
	AlgorithmBehaviour alg_behaviour;
	//alg_behaviour.costs_by_iteration = new double[num_of_iterations];
	alg_behaviour.elite_costs = new double[elite_size];
	alg_behaviour.elite_dataset = new double*[elite_size];

	for (int i = 0; i < elite_size; i++)
	{
		alg_behaviour.elite_dataset[i] = new double[problem_dimension];
	}
	

	return alg_behaviour;
}

EliteSet *createEliteSet(int problem_dimension, int elite_size) {
	EliteSet *elite;
	elite = (EliteSet *)malloc(sizeof(EliteSet));
	CURRENTSIZE = 0;
	elite->inf_threshold = INFINITY;

	elite->solutions = (Solution *)malloc(elite_size * sizeof(Solution));
	elite->_template = (double *)malloc(problem_dimension * sizeof(double));

	for (int i = 0; i < elite_size; i++) {
		//elite->solutions[i] = (Solution *)malloc(sizeof(Solution));
		elite->solutions[i].x = (double *)malloc(problem_dimension * sizeof(double));
	}

	return elite;
}

double isNearToOptimum(Funcao *func, double best, double precision)
{
	double optimum = func->getMinValue();
	double GAP = fabs(optimum - best);
	
	if (!optimum)
		return GAP <= precision;
	else 
		return GAP <= precision * fabs(optimum);
}

bool IsStopCriterionReached(int current_itr, double best_fo, ExperimentType exp_type, Funcao *func, bool &isSignClose) {
	switch (exp_type)
	{
	case BY_ITERATION: 
		return current_itr >= MAX_ITERATION;
		break;
	
	case BEFORE_AFTER:
		/* code */
		break;
	
	case BY_CFOs:

		return func->getFnEvals() >= MAX_nCFOs;
		break;
	
	case BY_GAP:
		if (isNearToOptimum(func, best_fo, 1.0e-3)) {
			//cout << "NEAR OPTIMUM" << endl;
			isSignClose = true;
			return true;
		}
		
		isSignClose = false;
		return false;
		//return current_itr >= MAX_ITERATION;
	default:
		break;
	}
}


bool phaseI_IsEnded(int current_itr, double best_fo, ExperimentType exp_type, Funcao *func, double percent) {
	switch (exp_type)
	{
	case BY_ITERATION: 
		return current_itr >= percent * MAX_ITERATION;
		break;
	
	case BEFORE_AFTER:
		/* code */
		break;
	
	case BY_CFOs:
		return func->getFnEvals() >= percent * MAX_nCFOs;
		break;
	
	case BY_GAP:
		return func->isNearOptimum(best_fo);
		break;
	
	
	default:
		break;
	}
}

double distance(map<int, double> & v1, map<int, double> & v2) {
	double sum = 0;
	map<int, double>::iterator it1 = v1.begin();
	map<int, double>::iterator it2 = v2.begin();
	while (it1 != v1.end()) {
		sum += (it1->second - it2->second) * (it1->second - it2->second);

		it1++;
		it2++;
	}
	
	return sqrt(sum);
}



double cgrasp(const char * mining_strategy, int isContinuousMining, double dmStartMoment, double patternPercentUsed, int eliteSize, int dimension, double *lower_bounds, double *upper_bounds, 
        Funcao *func, double hs, double he, double plo, int number_of_iterations, int max_functions_calls, unsigned long seed, bool &success, double standardDeviation)
{
	int i, n, count_calc_cost, imprc, imprl, ls_opt, size, iterations, function_evaluations, count = 0, countin, loop = 1, fecount = 0, stp = 0;
	double opt, h_s, h_e, ro, max_points, ep, cost_x, h, f_star, threshold, *x, *x_star, *best_solution, *l, *u;
	float t2;

	gen.seed(seed);

	success = false;
	int xmeansStrategy = 0; // 0 means dont use xmeans (clustering)
	
	FILE *outFile;

	string _mining_strategy(mining_strategy);
	if (_mining_strategy.find("dm") != string::npos)
		useDataMining = 1;
	else useDataMining = 0; // runs standard C-grasp 

	if (_mining_strategy.find("amx") != string::npos)
		xmeansStrategy = 4;
	else if (_mining_strategy.find("rmx") != string::npos)
		xmeansStrategy = 3;
	else if (_mining_strategy.find("mx") != string::npos)
		xmeansStrategy = 2;
	else if (_mining_strategy.find("x") != string::npos)
		xmeansStrategy = 1;



	//cout << _mining_strategy << " - " << xmeansStrategy << endl;

 	//--- initializing parameters ----------------------------------
	l = lower_bounds;
	u = upper_bounds;	
	
	h_s = hs;
	h_e = he;
	ro = plo;
	max_points = 100;
	double percent = 0.7;
	std_threshold = standardDeviation;

	success = false;
	iterations = number_of_iterations;
	MAX_ITERATION = iterations;
	g_elite_size = eliteSize;
	threshold = 1.0e-4;
	n = dimension;

	extraiu = 0;
	int print = 0; // Flag para controlar impressão dos resultados

	//--- initializing data structures ----------------------------------
	x = (double *)malloc(dimension * sizeof(double));
	x_star = (double *)malloc(dimension * sizeof(double));
	f_star = 1.0e+20;
	cost_x = 1.0e+20;


	EliteSet* elite;
	if (useDataMining)
		elite = createEliteSet(dimension, g_elite_size);
	
	int* fixedPositions = new int[dimension];
	_template = new double[dimension];
	for (int k = 0; k < n; k++) fixedPositions[k] = 0;

	int deviations_register_count = 0;

	//init_genrand(seed);
	//out << "seed: " << seed << endl;

	ExperimentType expType = BY_GAP;
	ExperimentType phaseIStopCriterion = BY_CFOs;
	MAX_ITERATION = number_of_iterations;
	MAX_nCFOs = max_functions_calls;

	bool eliteUpdated = false;
	vector<map<int, double>> patternsList, newPatterns;
	int currentPatternIndex = 0;
	double max_radius = 0;
	
	//debugLevel = DEBUG_LEVEL1_;
	int n_its = 0;
	while (!IsStopCriterionReached(count, f_star, expType, func, success)) 
	{
		n_its++;
		if (useDataMining) {
			if (phaseI_IsEnded(count, f_star, phaseIStopCriterion, func, dmStartMoment) && (!extraiu || isContinuousMining)) {
				if (!xmeansStrategy) {
					FindPatterns(elite, n, fixedPositions, l, u, h);
				} else {
					if (xmeansStrategy == 1) {
						FindPatternsWithXmeans(elite, n, fixedPositions, l, u, h, false);
					} else if (xmeansStrategy == 2)
						FindPatternsWithXmeans(elite, n, fixedPositions, l, u, h, true);
					else if (xmeansStrategy == 3)
						ChooseRandomPattern(elite, n, fixedPositions, l, u, h, true);
					else if (xmeansStrategy == 4) {
						if (eliteUpdated) {
							patternsList = ExtractPatterns(elite, n, fixedPositions, l, u, h, true, max_radius);
						}

						if (patternsList.size()) {

							if (currentPatternIndex >= patternsList.size())
								currentPatternIndex = 0;
							
							map<int, double> pattern = patternsList[currentPatternIndex];
							for (int i = 0; i < dimension; i++) {
								fixedPositions[i] = 0;
							}

							double* average = calculateAverageInElite(elite, g_elite_size, dimension);
							deviations = calculateDeviationsInElite(elite, g_elite_size, dimension, l, u, average);
							
							std::map<int, double>::iterator it = pattern.begin();
							while (it != pattern.end())
							{
								elite->_template[it->first] = it->second;

								if (std_threshold == 0 || deviations[it->first] < std_threshold) { // selecionar somente variaveis menos dispersas
									fixedPositions[it->first] = 1; 
								}

								it++;
							}
							currentPatternIndex++;
							
								
						}
							
					}
				}	
				
				if 	(eliteUpdated)
					copy_double_vector(_template, elite->_template, dimension);
				
				extraiu = 1;
			}
		}
		
		// if (!useDataMining)
		
		//cout << "Iteração " << count + 1 << endl;
		//cout << "Func Evals = " << func->getFnEvals() << " x Max CFOs = " <<  MAX_nCFOs << endl;
		//cout << "Best: " << f_star << endl;

		UnifRand(dimension, lower_bounds, upper_bounds, 
			x, fixedPositions, _template, patternPercentUsed);	
		

		h = h_s;
		//cout << "h_s: " << h_s << " - h_e: " << h_e << endl;
		while (h >= h_e) 
		{

			if (IsStopCriterionReached(count, f_star, expType, func, success))
				break;

			res = f_star;
			imprc = 0;
			imprl = 0;

			if (debugLevel == DEBUG_LEVEL2_) 
			{
				cost_x =  func->calc(x);
				printf("begin: ");
				Util::printX(x, n); printf(" - %lf ", cost_x);
				printf("\n");
			}


			construct_greedy_randomized(n, ro, h, l, u, x, &imprc, t2, x_star, f_star, opt, outFile, &count_calc_cost, &fecount, stp, function_evaluations, NULL, func);
			cost_x = func->calc(x);
			//cout << "build: " << cost_x << endl;

			if (debugLevel == DEBUG_LEVEL1_) 
			{
				printf("build: ");
				Util::printX(x, n); printf(" - %lf ", cost_x);
				printf("\n");
			}

			local_search(n, ro, h, l, u, max_points, x, &imprl,	t2, x_star, f_star, opt, outFile, &count_calc_cost, &fecount, stp, function_evaluations, func);
			cost_x = func->calc(x) ;
			//cout << "local_search: " << cost_x << endl;

			if (debugLevel == DEBUG_LEVEL1_) 
			{
				printf("local search: ");
				Util::printX(x, n); printf(" - %lf ", cost_x);
				printf("\n");
			}
			if (cost_x < f_star) {
				f_star = cost_x;
				copy_double_vector(x_star, x, dimension);
			}

			//if (debugLevel == DEBUG_LEVEL1_) 
			//	printf("%d - %d \n", imprc, imprl);

			if ((imprc == 0) && (imprl == 0)) {
				h = h / 2.0;
			}

			//cout << "h = " << h << endl;
			
			if (useDataMining) {
				eliteUpdated = updateEliteSet(x, cost_x, elite, n);
			}
		}

		//cout << "iteration " << count + 1 << " - best: " << f_star  << " - founded: " << cost_x << endl;
		//Util::printX(x_star, n); cout << endl;
		count++;
	}

	cout << "Nº of Iterations: " << n_its << endl;
	//cout << "best: " << f_star;
	free(x);
	free(x_star);
	free(fixedPositions);
	free(_template);

	return f_star;
}
