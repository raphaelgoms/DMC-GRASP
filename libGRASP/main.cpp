#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <string>

#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

//#include "DMCGrasp.h"
#include "Funcao.h"
#include "MersenneTwister.h"
#include "Util.h"

#include "Rosenbrock2.h"
#include "Zakharov.h"
#include "SumSquares.h"
#include "Branin.h"
#include "Easom.h"
#include "GoldsteinPrice.h"
#include "Shekel.h"
#include "Hartmann.h"
#include "Shubert.h"
#include "Beale.h"
#include "Booth.h"
#include "Bohachevsky.h"
#include "Hump.h"
#include "Matyas.h"
#include "Schwefel.h"

#include "Colville.h"
#include "Perm.h"
#include "Perm0.h"
#include "PowerSum.h"
#include "Griewank.h"
#include "Rastrigin.h"
#include "Trid.h"
#include "Powell.h"
#include "DixonPrice.h"
#include "Ackley.h"
#include "Levy.h"
#include "Sphere.h"

#include "cgrasp.h"
#include <sys/resource.h>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <map>

#include <sys/stat.h>
int exist(const char *name)
{
  struct stat   buffer;
  return (stat (name, &buffer) == 0);
}

using namespace std;
#define PURO 	0
#define HIBRIDO 1
#define BFGS 	2
#define DATAMINING 3


 int64_t getMilisegundos() {
	struct timeval tempo;
	gettimeofday(&tempo, NULL);
	return (int64_t) tempo.tv_sec * 1000 + tempo.tv_usec / 1000; 
 }


typedef struct {
	double *u, *l;
	double hs, he, plo;
} Parameters;


Parameters getParameters(int iFuncNum, int n, Funcao **func){
	int m;
	double *u, *l;
	double hs, he, plo;
	l = new double[n];
	u = new double[n];
	int exp_number = 1;

	plo = 0.7;
				
	switch(iFuncNum)
	{
		case Funcao::ROSENBROCK: 
				hs = 0.5;
				he = 0.0001;
				if (n == 20) {
					hs = 0.1;
					he = 0.05;
				}

				for (int i =0; i < n; i++){
					l[i] = -10.0;
					u[i] = 10.0;
				}
				*func = new Rosenbrock2(n);
				break;			
		case Funcao::ZAKHAROV:
				if (n == 5) {
					hs = 1.0;
					he = 0.001;
				}			
				else if (n == 10) {
					hs = 1.0;
					he = 0.001;
				}			
				else if (n == 20){			
					hs = 1.0;
					he = 0.0001;
				} else {
					hs = 1.0;
					he = 0.01;
				}

				for (int i =0; i < n; i++){
					l[i] = -5.0;
					u[i] = 10.0;
				}
				*func = new Zakharov(n);
				break;
		case Funcao::SUMSQUARES:
				hs = 1.0;
				he = 0.001;
				for (int i =0; i < n; i++){
					l[i] = -5.0;
					u[i] = 10.0;
				}
				*func = new SumSquares(n);
				break;	
		case Funcao::BRANIN:
				hs = 1.0;
				he = 0.02;

				l[0] = -5; l[1] = -5;
			    u[0] = 15; u[1] = 15;
				*func = new Branin(n);
				break;
	
		case Funcao::GOLDSTEINPRICE:
			    hs = 1.0;
				he = 0.001;
				for (int i =0; i < n; i++){
					l[i] = -2.0;
					u[i] = 2.0;
				}
				*func = new GoldsteinPrice(n);
				break;

		case Funcao::EASOM:
			 	hs = 1;
				he = 0.01;
				
			    for (int i =0; i < n; i++){
					l[i] = -100.0;
					u[i] = 100.0;
				}
				*func = new Easom(n);
				break;
	
		case Funcao::SHEKEL:
				m = n;			    
				n = 4;

				//cout << "m = " << m << " n = " << n << endl;       

				delete []l;
				delete []u;
					
				l = new double[n];				
				u = new double[n];				

				hs = 1;
				he = 0.005;

				if (m==10){
					hs = 1;
					he = 0.001;
					//cout << "m = " << m << " n = " << n << endl; 
				}

				//hs = 0.1;
				//he = 0.05;
				for (int i=0; i < n; i++){
					l[i] = 0.0;
					u[i] = 10.0;
				}
				*func = new Shekel(n,m);
				break;

		case Funcao::HARTMANN:
				hs = 0.5;
				//he = 0.001;
				if (n == 3)
					he = 0.005;
				else 
					he = 0.001;		
 
				for (int i =0; i < n; i++){
					l[i] = 0.0;
					u[i] = 1.0;
				}
				*func = new Hartmann(n,4);
				break;

		case Funcao::SHUBERT:
				hs = 1.0;
				he = 0.01;
				//hs = 0.1;
				//he = 0.05;
				for (int i =0; i < n; i++){
					l[i] = -10.0;
					u[i] = 10.0;
				}
				*func = new Shubert(n);
				break;
	
		case Funcao::BEALE:
				hs = 1.0;
				he = 0.05;
				for (int i =0; i < n; i++){
					l[i] = -4.5;
					u[i] = 4.5;
				}
				*func = new Beale(n);
				break;

		case Funcao::BOOTH:
				hs = 1.0;
				he = 0.05;
				for (int i =0; i < n; i++){
					l[i] = -10.0;
					u[i] = 10.0;
				}
				*func = new Booth(n);
				break;

		case Funcao::BOHACHEVSKY:
				hs = 1.0;
				he = 0.01;
				for (int i =0; i < n; i++){
					l[i] = -50.0;
					u[i] = 100.0;
				}
				*func = new Bohachevsky(n);
				break;

		case Funcao::HUMP:
				hs = 1.0;
				he = 0.01;
				for (int i =0; i < n; i++){
					l[i] = -5.0;
					u[i] = 5.0;
				}
				*func = new Hump(n);
				break;

		case Funcao::MATYAS:
				hs = 1.0;
				he = 0.1;
				for (int i =0; i < n; i++){
					l[i] = -5.0;
					u[i] = 10.0;
				}
				*func = new Matyas(n);
				break;

		case Funcao::SCHWEFEL:
				if (n == 2){	
					hs = 5.0;
					he = 0.001;;
				}
				else if (n==6){
					hs = 50.0;
					he = 0.001;;
				}

				for (int i =0; i < n; i++){
					l[i] = -500.0;
					u[i] = 500.0;
				}
				*func = new Schwefel(n);
				break;
	
		case Funcao::COLVILLE:
				hs = 1.0;
				he = 0.0005;
				for (int i =0; i < n; i++){
					l[i] = -10.0;
					u[i] = 10.0;
				}
				*func = new Colville(n);
				break;

		case Funcao::PERM:
				//hs = 1.0;
				//he = 0.01;
				hs = 0.1;
				he = 0.0001;
				for (int i =0; i < n; i++){
					l[i] = (double) -n;
					u[i] = (double) n;
				}
				*func = new Perm(n);
				break;

		case Funcao::PERM0:
				hs = 1.0;
				he = 0.0001;
				//hs = 0.1;
				//he = 0.0125;
				for (int i =0; i < n; i++){
					l[i] = (double) -n;
					u[i] = (double) n;
				}
				*func = new Perm_0(n);
				break;

		case Funcao::POWERSUM:
				hs = 1.0;
				he = 0.01;
				for (int i =0; i < n; i++){
					l[i] = (double) 0.0;
					u[i] = (double) n;
				}
				*func = new PowerSum(n);
				break;


		case Funcao::GRIEWANK:
				hs = 10.0;
				he = 0.001;
				for (int i =0; i < n; i++){
					l[i] = -300.0;
					u[i] = 600.0;
				}
				*func = new Griewank(n);
				break;

		case Funcao::RASTRIGIN:
				//hs = 1.0;
				//he = 0.1;
				hs = 0.5;
				he = 0.0005;
				for (int i =0; i < n; i++){
					l[i] = -2.56;
					u[i] = 5.12;
				}
				*func = new Rastrigin(n);
				break;


		case Funcao::TRID:
				if (n == 10){
					hs = 5.0;
					he = 0.1;
				}
				else{
					hs = 1.0;
					he = 0.1;
				}

				for (int i =0; i < n; i++){
					l[i] = (double) -n*n;
					u[i] = (double) n*n;
				}

				*func = new Trid(n);
				break;

		case Funcao::POWELL:
				hs = 1.0;
				he = 0.0005;
				for (int i =0; i < n; i++){
					l[i] = -4.0;
					u[i] = 5.0;
				}
				*func = new Powell(n);
				break;

		case Funcao::DIXONPRICE:
				hs = 5.0;
				he = 0.0005;//25;
				for (int i =0; i < n; i++){
					l[i] = -10.0;
					u[i] = 10.0;
				}
				*func = new DixonPrice(n);
				break;

		case Funcao::ACKLEY:
				hs = 5.0;
				he = 0.0005;
				for (int i =0; i < n; i++){
					l[i] = -15.0;
					u[i] = 30.0;
				}
				*func = new Ackley(n);
				break;


		case Funcao::LEVY:
				hs = 1.0;
				he = 0.0005;
				for (int i =0; i < n; i++){
					l[i] = -10.0;
					u[i] = 10.0;
				}
				*func = new Levy(n);
				break;

		case Funcao::SPHERE:
				hs = 2.0;
				he = 0.005;
				for (int i =0; i < n; i++){
					l[i] = -2.56;
					u[i] = 5.12;
				}
				*func = new Sphere(n);
				break;


	}

	hs = 0.5;

	// if (n < 10) {
	// 	he = 0.01;
	// } else {
	// 	he = 0.001;
	// }

	//he = 0.01;
	he = 0.0005;

	

	Parameters parameters;
	parameters.l = l;
	parameters.u = u;
	parameters.hs = hs;
	parameters.he = he;
	parameters.plo = plo;

	return parameters;
 }


 void usage(){
	printf("\nForma de uso: ");
	printf(" CGrasp FN N MAX_EVALS RUNS TYPE M\n");
	printf(" FN: nome da funcao (ex: ROSENBROCK, ZAKHAROV, etc)\n");
	printf(" N: numero de dimensoes\n");
	printf(" RUNS: numero de execucoes do CGrasp p\n");
	printf(" MAX_EVALS: numero maximo de avaliacoes da funcao. Caso 0 rodar até encerrar o algoritmo\n");
	printf(" TYPE: tipo (PURO, HIBRIDO)\n");
 	printf(" M: numero de correcoes da matriz do L-BFGS (m < n)\n\n");
 }

 int getFuncNumb(char *funcName){
	
	if (!strcmp(funcName, "ROSENBROCK")) 	return Funcao::ROSENBROCK;
	if (!strcmp(funcName, "ZAKHAROV")) 		return Funcao::ZAKHAROV;
	if (!strcmp(funcName, "SUMSQUARES")) 	return Funcao::SUMSQUARES;
	if (!strcmp(funcName, "BRANIN")) 		return Funcao::BRANIN;
	if (!strcmp(funcName, "EASOM")) 		return Funcao::EASOM;
	if (!strcmp(funcName, "GOLDSTEINPRICE"))return Funcao::GOLDSTEINPRICE;
	if (!strcmp(funcName, "SHEKEL")) 		return Funcao::SHEKEL;
	if (!strcmp(funcName, "HARTMANN")) 		return Funcao::HARTMANN;
	if (!strcmp(funcName, "SHUBERT")) 		return Funcao::SHUBERT;
	if (!strcmp(funcName, "BEALE"))			return Funcao::BEALE;
	if (!strcmp(funcName, "BOOTH"))			return Funcao::BOOTH;
	if (!strcmp(funcName, "BOHACHEVSKY"))	return Funcao::BOHACHEVSKY;
	if (!strcmp(funcName, "HUMP"))			return Funcao::HUMP;
	if (!strcmp(funcName, "MATYAS"))		return Funcao::MATYAS;
	if (!strcmp(funcName, "SCHWEFEL"))		return Funcao::SCHWEFEL;
	if (!strcmp(funcName, "COLVILLE"))		return Funcao::COLVILLE;
	if (!strcmp(funcName, "PERM"))			return Funcao::PERM;
	if (!strcmp(funcName, "PERM0"))			return Funcao::PERM0;
	if (!strcmp(funcName, "POWERSUM"))		return Funcao::POWERSUM;
	if (!strcmp(funcName, "GRIEWANK"))		return Funcao::GRIEWANK;
	if (!strcmp(funcName, "RASTRIGIN"))		return Funcao::RASTRIGIN;
	if (!strcmp(funcName, "TRID"))			return Funcao::TRID;
	if (!strcmp(funcName, "POWELL"))		return Funcao::POWELL;
	if (!strcmp(funcName, "DIXONPRICE"))	return Funcao::DIXONPRICE;
	if (!strcmp(funcName, "ACKLEY"))		return Funcao::ACKLEY;
	if (!strcmp(funcName, "LEVY"))			return Funcao::LEVY;
	if (!strcmp(funcName, "SPHERE"))		return Funcao::SPHERE;

	return -1;
 }

 int getType(char *typeName){
	
	if (!strcmp(typeName, "PURO")) 		return PURO;
	if (!strcmp(typeName, "BFGS")) 		return BFGS;
	if (!strcmp(typeName, "HIBRIDO")) 	return HIBRIDO;
	if (!strcmp(typeName, "DATAMINING")) 	return DATAMINING;
	
	return -1;
 }
 
 void saveGaps(double *mediaGaps, int numIter){
		FILE *arqlog100 = fopen("resultados/gap100.txt", "a+");		
		FILE *arqlog500 = fopen("resultados/gap500.txt", "a+");		
		FILE *arqlog1000 = fopen("resultados/gap1000.txt", "a+");		
		FILE *arqlog5000 = fopen("resultados/gap5000.txt", "a+");		
		FILE *arqlog10000 = fopen("resultados/gap10000.txt", "a+");		
		FILE *arqlog20000 = fopen("resultados/gap20000.txt", "a+");		
		FILE *arqlog50000 = fopen("resultados/gap50000.txt", "a+");		

		/*for (int j = 0; j < 7; j++){		
			printf("Media do GAP = %lf \n", mediaGaps[j]/(double)numIter);
			fprintf(arqlog, "%lf (%s-%d) - num_aval = %d \n", mediaGaps[j]/(double)numIter, iFuncName, n, j);
		}*/		
		printf("Media do GAP 100 = %lf \n", mediaGaps[0]/(double)numIter);
		fprintf(arqlog100, "%lf\n", mediaGaps[0]/(double)numIter);
	
		printf("Media do GAP 500 = %lf \n", mediaGaps[1]/(double)numIter);
		fprintf(arqlog500, "%lf\n", mediaGaps[1]/(double)numIter);

		printf("Media do GAP 1000 = %lf \n", mediaGaps[2]/(double)numIter);
		fprintf(arqlog1000, "%lf\n", mediaGaps[2]/(double)numIter);

		printf("Media do GAP 5000 = %lf \n", mediaGaps[3]/(double)numIter);
		fprintf(arqlog5000, "%lf\n", mediaGaps[3]/(double)numIter);

		printf("Media do GAP 10000 = %lf \n", mediaGaps[4]/(double)numIter);
		fprintf(arqlog10000, "%lf\n", mediaGaps[4]/(double)numIter);

		printf("Media do GAP 20000 = %lf \n", mediaGaps[5]/(double)numIter);
		fprintf(arqlog20000, "%lf\n", mediaGaps[5]/(double)numIter);

		printf("Media do GAP 50000 = %lf \n", mediaGaps[6]/(double)numIter);
		fprintf(arqlog50000, "%lf\n", mediaGaps[6]/(double)numIter);


		fclose(arqlog100);
		fclose(arqlog500);
		fclose(arqlog1000);
		fclose(arqlog5000);
		fclose(arqlog10000);
		fclose(arqlog20000);
		fclose(arqlog50000);
}



char *build_path(char * iFuncName, char *dimension, char *algorithm)
{
	char* pathname = "";

	strcat(pathname, algorithm);
	strcat(pathname, "-");
	strcat(pathname, iFuncName);
	strcat(pathname, "-");
	strcat(pathname, dimension);
	strcat(pathname, "-");
	
	return pathname;
}

void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total)
{
  long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
  struct rusage ptempo;

  getrusage(0,&ptempo);

  seg_CPU = ptempo.ru_utime.tv_sec;
  mseg_CPU = ptempo.ru_utime.tv_usec;
  seg_sistema = ptempo.ru_stime.tv_sec;
  mseg_sistema = ptempo.ru_stime.tv_usec;

 *seg_CPU_total     = (seg_CPU + 0.000001 * mseg_CPU);
 *seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
}

void createFolder(string path) {
	bool success = mkdir(path.c_str(), 0777);
}

ofstream openFile(string path) {
	ofstream file;
	if (!exist(path.c_str())) 
		file.open(path);
	else
		file.open(path, ios::app);

	return file;
}

void saveElite(string rootpath, int numberofiterations, int problemdimension, int elitesize, string funcname, AlgorithmBehaviour &behaviour){

	// Save elite by iteration:
	string elite_by_iter_filepath(rootpath + "/" + funcname + "/");
	createFolder(elite_by_iter_filepath);

	elite_by_iter_filepath += "BY_ITERTATION/";
	createFolder(elite_by_iter_filepath);
	
	ofstream elite_by_iter_file;

	for (int i = 0; i < behaviour.elite_by_iteration.size(); i++)
	{
		elite_by_iter_file = openFile(elite_by_iter_filepath + "ITERATION_" + to_string(i+1) + ".csv");

		for (int i = 0; i < problemdimension; i++) 
			elite_by_iter_file << "x" << i+1 << ";";  

		elite_by_iter_file << "Cost" << endl;
		for (int j = 0; j < elitesize; j++)
		{
			for (int l = 0; l < problemdimension; l++)
			{
				elite_by_iter_file  << behaviour.elite_by_iteration[i][j][l] << ";";
			}
			elite_by_iter_file << behaviour.elite_costs_by_iteration[i][j] << endl;
		}
		elite_by_iter_file.close();
	}
	
	
	// Save elite whole history 
	string elite_history_filepath(rootpath + "/" + funcname + "/HISTORY/");
	createFolder(elite_history_filepath);

	ofstream elite_history_file;

	for (int i = 0; i < behaviour.elite_history.size(); i++)
	{
		elite_history_file = openFile(elite_history_filepath + "ITERATION_" + to_string(i+1) + ".csv");

		for (int i = 0; i < problemdimension; i++) 
			elite_history_file << "x" << i+1 << ";";  

		elite_history_file << "Cost" << endl;
		for (int j = 0; j < elitesize; j++)
		{
			for (int l = 0; l < problemdimension; l++)
			{
				elite_history_file  << behaviour.elite_history[i][j][l] << ";";
			}
			elite_history_file << behaviour.elite_costs_history[i][j] << endl;
		}
	}

	elite_history_file.close();

	//save elite just before mining
	string elite_bf_mining_filepath(rootpath + "/" + funcname + "/");
	ofstream elite_bf_mining_file;

	elite_bf_mining_file = openFile(elite_bf_mining_filepath + "ELITE_BEFORE_MINING.csv");

	for (int i = 0; i < problemdimension; i++) 
		elite_bf_mining_file << "x" << i+1 << ";";  

	elite_bf_mining_file << "Cost" << endl;
	for (int j = 0; j < elitesize; j++)
	{
		for (int l = 0; l < problemdimension; l++) {
			elite_bf_mining_file  << behaviour.elite_dataset[j][l] << ";";
		}

		elite_bf_mining_file << behaviour.elite_costs[j] << endl;
	}

	elite_bf_mining_file.close();
}

void saveBehaviour(string EXPTYPE, AlgorithmBehaviour behavior, int number_of_iterations, int elite_size, string func, int dimension, string alg_prefix, int seed) 
{
	ofstream before_after_file;
	string before_after_path("./results/"+EXPTYPE+"/"+ alg_prefix +"/before_after_" + to_string(seed) + ".csv");
	if (!exist(before_after_path.c_str())) {
		before_after_file.open("./results/"+EXPTYPE+"/"+ alg_prefix +"/before_after_" + to_string(seed) + ".csv");
		before_after_file << "Func;Dim;" ;
		for (int j=0; j < number_of_iterations-1; j++)
			before_after_file << "ITR_" << j+1 << ";";
		before_after_file << "ITR_" << number_of_iterations << endl;
	} else {
		before_after_file.open("./results/"+EXPTYPE+"/"+ alg_prefix +"/before_after_" + to_string(seed) + ".csv", ios::app);
	}

	createFolder("./results/"+EXPTYPE+"/"+ alg_prefix +"/solutions");
	createFolder("./results/"+EXPTYPE+"/"+ alg_prefix +"/solutions/" + to_string(seed));
	createFolder("./results/"+EXPTYPE+"/"+ alg_prefix +"/solutions/" + to_string(seed) + "/HISTORY");

	ofstream solutions_history_file;
	solutions_history_file = openFile("./results/"+EXPTYPE+"/"+ alg_prefix +"/solutions/" + to_string(seed) + "/HISTORY/" + func + "_" + to_string(dimension) + ".csv");

	for (int i = 0; i < dimension; i++) 
		solutions_history_file << "x" << i+1 << ";";  

	solutions_history_file << "Cost" << endl;


	for (int i = 0; i < behavior.solutions_history.size(); i++)
	{
		for (int l = 0; l < dimension; l++)
		{
			solutions_history_file  << behavior.solutions_history[i][l] << ";";
		}
		solutions_history_file << behavior.costs_history[i] << endl;
	}
	solutions_history_file.close();

	createFolder("./results/"+EXPTYPE+"/"+ alg_prefix +"/solutions/" + to_string(seed) + "/BY_ITERATION");

	ofstream solutions_by_iteration_file;
	solutions_by_iteration_file = openFile("./results/"+EXPTYPE+"/"+ alg_prefix +"/solutions/" + to_string(seed) + "/BY_ITERATION/" + func + "_" + to_string(dimension) +  ".csv");
	
	for (int i = 0; i < dimension; i++) 
		solutions_by_iteration_file << "x" << i+1 << ";";  

	solutions_by_iteration_file << "Cost" << endl;

	for (int i = 0; i < behavior.solutions_by_iteration.size(); i++)
	{
		for (int l = 0; l < dimension; l++) {
			solutions_by_iteration_file  << behavior.solutions_by_iteration[i][l] << ";";
		}
		solutions_by_iteration_file << behavior.costs_by_iteration[i] << endl;
	}

	solutions_by_iteration_file.close();

	before_after_file << func << ";" << dimension << ";";
	for (int j=0; j < number_of_iterations-1; j++)
		before_after_file << behavior.costs_by_iteration[j] << ";";
	before_after_file << behavior.costs_by_iteration[number_of_iterations-1] << endl;
	before_after_file.close();

	ofstream elite_file;
	string elite_file_name = "./results/"+EXPTYPE+"/"+ alg_prefix +"/elite";
	bool success = mkdir(elite_file_name.c_str(), 0777);
	elite_file_name += "/" +to_string(seed);
	success = mkdir(elite_file_name.c_str(), 0777);

	elite_file_name += "/DEV_" + func + "-" + to_string(dimension) + ".csv";
	elite_file.open(elite_file_name);

	for (int i = 0; i < dimension; i++) {
		elite_file << "x" << i+1 << ";"; 
	}
	
	elite_file << "BestCost;WorstCost;AvgCost" << endl;
	for (int i = 0; i < behavior.deviations_in_elite_by_dimension.size(); i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			elite_file  << behavior.deviations_in_elite_by_dimension[i][j] << ";";
		}

		elite_file << behavior.elite_best_costs[i] << ";";
		elite_file << behavior.elite_worst_costs[i] << ";";
		elite_file << behavior.elite_avg_costs[i] << endl;
	}

	elite_file.close();	

	saveElite("./results/"+EXPTYPE+"/"+ alg_prefix +"/elite/" + to_string(seed), number_of_iterations, dimension, elite_size, func, behavior);
}

string getLastLineFromFile(string path)
{
	fstream fs;
	fs.open(path);
	std::string lastline;
	if(fs.is_open())
	{
		fs.seekg(-1, std::ios_base::end);
		if(fs.peek() == '\n')
		{
			fs.seekg(-1, std::ios_base::cur);
			int i = fs.tellg();
			for(i;i > 0; i--)
			{
				if(fs.peek() == '\n')
				{
				//Found
				fs.get();
				break;
				}
				fs.seekg(i, std::ios_base::beg);
			}
		}
	
		getline(fs, lastline);
	}

	return lastline;
}

std::vector<int> getVectorFromFile(string path)
{
	fstream fs;
	fs.open(path);
	std::vector<int> vector;

	if(fs.is_open())
	{
		std::string line;
		while(getline(fs, line))
		{
			vector.push_back(atoi(line.c_str()));
		}
	}		
	return vector;
}

map<std::string, std::vector<int>> createMaxCFObySeedMap(std::string funcName) 
{
	map<std::string, std::vector<int>> maxCFObySeedMap;
	std::vector<int> maxCfos = getVectorFromFile("./CFOs/cfos_+" + funcName);
	maxCFObySeedMap.insert(pair<std::string,std::vector<int>>("funcName", maxCfos));

	return maxCFObySeedMap;
}

map<std::string, int> createMaxCFOMap() {
	map<std::string, int> *max_cfo_map = new map<std::string, int>();
	max_cfo_map->insert(pair<std::string, int>("BEALE2", 504005));
	max_cfo_map->insert(pair<std::string, int>("BOHACHEVSKY2", 190894));
	max_cfo_map->insert(pair<std::string, int>("BOOTH2", 15044));
	max_cfo_map->insert(pair<std::string, int>("BRANIN2", 14318));
	max_cfo_map->insert(pair<std::string, int>("EASOM2", 2.02724e+07));
	max_cfo_map->insert(pair<std::string, int>("GOLDSTEINPRICE2", 71087));
	max_cfo_map->insert(pair<std::string, int>("MATYAS2", 1850));
	max_cfo_map->insert(pair<std::string, int>("HUMP2", 8314));
	max_cfo_map->insert(pair<std::string, int>("ROSENBROCK2", 6978));
	max_cfo_map->insert(pair<std::string, int>("SCHWEFEL2", 166836));
	max_cfo_map->insert(pair<std::string, int>("SHUBERT2", 19425));
	max_cfo_map->insert(pair<std::string, int>("ZAKHAROV2", 10332));
	max_cfo_map->insert(pair<std::string, int>("SPHERE3", 8793));
	max_cfo_map->insert(pair<std::string, int>("HARTMANN3", 28678));
	max_cfo_map->insert(pair<std::string, int>("COLVILLE4",272563));
	max_cfo_map->insert(pair<std::string, int>("PERM4", 3.746e+06));
	max_cfo_map->insert(pair<std::string, int>("PERM04", 3.15803e+06));
	max_cfo_map->insert(pair<std::string, int>("POWERSUM4", 1.92637e+06));
	max_cfo_map->insert(pair<std::string, int>("SHEKEL45", 4.9803e+06));
	max_cfo_map->insert(pair<std::string, int>("SHEKEL47", 6.47833e+06));
	max_cfo_map->insert(pair<std::string, int>("SHEKEL410", 8.90985e+06));
	max_cfo_map->insert(pair<std::string, int>("ROSENBROCK5", 1.09573e+07));
	max_cfo_map->insert(pair<std::string, int>("ZAKHAROV5", 27199));
	max_cfo_map->insert(pair<std::string, int>("HARTMANN6", 10356));
	max_cfo_map->insert(pair<std::string, int>("SCHWEFEL6", 2.22138e+06));
	max_cfo_map->insert(pair<std::string, int>("TRID6", 15865));
	max_cfo_map->insert(pair<std::string, int>("GRIEWANK10", 1.31649e+07));
	max_cfo_map->insert(pair<std::string, int>("RASTRIGIN10",600119));
	max_cfo_map->insert(pair<std::string, int>("ROSENBROCK10", 2.21233e+07));
	max_cfo_map->insert(pair<std::string, int>("SUMSQUARES10", 249864));
	max_cfo_map->insert(pair<std::string, int>("TRID10", 110000));
	max_cfo_map->insert(pair<std::string, int>("ZAKHAROV10", 91973));
	max_cfo_map->insert(pair<std::string, int>("GRIEWANK20", 4.48394e+07));
	max_cfo_map->insert(pair<std::string, int>("RASTRIGIN20", 3.78949e+06));
	max_cfo_map->insert(pair<std::string, int>("ROSENBROCK20", 1.11063e+07));
	max_cfo_map->insert(pair<std::string, int>("SUMSQUARES20", 1.62598e+06));
	max_cfo_map->insert(pair<std::string, int>("ZAKHAROV20", 654126));
	max_cfo_map->insert(pair<std::string, int>("POWELL24", 18808));
	max_cfo_map->insert(pair<std::string, int>("DIXONPRICE25", 1.52986e+06));
	max_cfo_map->insert(pair<std::string, int>("ACKLEY30", 5.50139e+07));
	max_cfo_map->insert(pair<std::string, int>("LEVY30", 4.82452e+06));
	max_cfo_map->insert(pair<std::string, int>("SPHERE30", 3.79675e+06));

	return *max_cfo_map;
}


bool isAlgorithmCode(const char *str, string &algorithmName)
{
	
	if (!strcmp("--dmc", str)) {
		algorithmName = "dmc";
		return true;
	} else if (!strcmp("--xdmc", str)) {
		algorithmName = "xdmc";
		return true;
	} else if (!strcmp("--mxdmc", str)) {
		algorithmName = "mxdmc";
		return true;
	}  else if (!strcmp("--c", str)) {
		algorithmName = "c";
		return true;
	} else if (!strcmp("--rmxdmc", str)) {
		algorithmName = "rmxdmc";
		return true;
	} else if (!strcmp("--amxdmc", str)) {
		algorithmName = "amxdmc";
		return true;
	}

	algorithmName = "";
	return false;

}

int getMaxCFOsNumber(std::string filename, std::string func_name) {

	std::ifstream myFile(filename);
    if(!myFile.is_open()) cout << "Could not open file" << endl;

	int CFOsColumn = -1;
	std::string line, colname, val;
    int col_number;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);
        std::stringstream ss(line);

        // Extract each column name
		col_number = 0;
        while(std::getline(ss, colname, ';')){
			if (colname == "CFOs")
				CFOsColumn = col_number;

			col_number++;
        }
    }

	while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        while(std::getline(ss, val, ';')){
			//cout << "val: "<<  val << endl;
			if (colIdx == 0 && val != func_name) {
				break;
			}	
			
			if (colIdx == CFOsColumn) {
				return floor(stof(val));
			}

			//if(ss.peek() == ';') ss.ignore();
            colIdx++;
        }
    }

	myFile.close();
	return -1;
}




int main(int argc, char **argv)
{	
	string algCode;
	char *funcName;
	int iFuncNumb;
	int problemDim;
	int seed;
	int eliteSize;
	double dmStartMoment;
	int dmFreqStrategy; 
	double patternPercentUsed = 0.0;
	double standardDeviation = 0.0;

	string arguments("");

	bool print = false;

	map<std::string, int> max_cfo_map = createMaxCFOMap();	
	
	int i = 0;
	bool alg_found = false;
	while (++i < argc) {
		
		if (!alg_found && isAlgorithmCode(argv[i], algCode)) {
			alg_found = true;
			if (print) cout << "Algorithm: " << algCode << endl;
			arguments += algCode + " ";
		}
		else
		if (!strcmp("-i", argv[i])) {
			funcName = argv[++i];
			boost::to_upper(funcName);
			iFuncNumb = getFuncNumb(funcName);
			if (print) cout << "Instance: " << iFuncNumb << endl;
			//arguments += "-i " + string(argv[i])  + " ";
		}
		else
		if (!strcmp("--nvar", argv[i])) {
			problemDim = atoi(argv[++i]);
			if (print) cout << "Number of variables: " << problemDim << endl;
			//arguments += "--nvar " + string(argv[i])  + " ";
		} 
		else
		if (!strcmp("--seed", argv[i])) {
			seed = atoi(argv[++i]);
			if (print) cout << "Seed: " << seed << endl;
		} 
		else
		if (!strcmp("--elsz", argv[i])) {
			eliteSize = atoi(argv[++i]);
			if (print) cout << "Elite size: " << eliteSize << endl;
			arguments += "--elsz " +  string(argv[i])  + " ";
		}

		if (!strcmp("--dmstart", argv[i])) {
			dmStartMoment = atof( argv[++i]);
			if (print) cout << "DM start moment: " << dmStartMoment << endl;
			arguments += "--dmstart " +  string(argv[i])  + " ";
		}

		if (!strcmp("--dmfreq", argv[i])) {
			
			if (!strcmp("jot", argv[++i])) {
				dmFreqStrategy = 0; // just one time
			} else {
				dmFreqStrategy = 1; // continuous
			}

			if (print) cout << "DM frequency: " << argv[i] << endl;
			arguments += "--dmfreq " + string(argv[i])  + " ";
		}

		if (!strcmp("--ptsz", argv[i])) {
			patternPercentUsed = atof( argv[++i]);
			if (print) cout << "Percent of pattern used: " << patternPercentUsed << endl;
			arguments += "--ptsz " + string(argv[i]) + " ";
		}

		if (!strcmp("--sd", argv[i])) {
			
			standardDeviation = atof( argv[++i]);
			if (print) cout << "Standard Deviation used: " << standardDeviation << endl;
			arguments += "--sd " +  string(argv[i])  + " ";
		}

	}

	string funcCode = funcName;
	if (iFuncNumb == Funcao::SHEKEL) {
		funcCode += "4" + to_string(problemDim);
		//problemDim = 4;
	} else {
		funcCode += to_string(problemDim);
	}
	

	int max_cfos = 0;
	if (algCode.compare("c"))
		max_cfos = getMaxCFOsNumber("./out(c ).csv", funcCode);

	//cout << funcCode << ": " << max_cfos << endl;
	//return 0;
	//int max_cfos = max_cfo_map[funcCode];

	bool success = 0;
	int number_of_iterations = 20;
	
	double soma = 0.0;
	double min = 1.0e+30;
	double result;
	double avg_cfos = 0;
	double avg_time = 0;
	double success_rate = 0;
	double s_CPU_inicial, s_CPU_final;
  	double s_total_inicial, s_total_final;
	
	bool save_results = true;
	//std::ofstream ttt_time_file("./tttplots_data/" + algCode +"_" + funcCode + "_time_file.dat");	
	//std::ofstream ttt_cfo_file("./tttplots_data/" +algCode +"_" + funcCode + "_cfo_file.dat");	
	std::ofstream fo_file("./fobj/" +algCode +"_" + funcCode + "_fo_file.dat");
	int num_runs = 30;
	for (int i = 0; i < num_runs; i++)
	{
		seed = 270000 + i + 1;

		Funcao *func;
		Parameters parameters = getParameters(iFuncNumb, problemDim, &func);

		int dimension = problemDim;
		if (iFuncNumb == Funcao::SHEKEL) 
			dimension = 4;

		//double x[5] = { 1, 1, 1, 1, 1 };
		//cout << "Best: " << func->calc(x) << endl;

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		result = cgrasp(funcCode.c_str(), algCode.c_str(), dmFreqStrategy, dmStartMoment, patternPercentUsed, eliteSize, dimension, parameters.l, parameters.u, func, parameters.hs, 
			parameters.he, parameters.plo, number_of_iterations, max_cfos, seed, success, standardDeviation);
		Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);

		if (result < min)
			min = result;

		if (save_results) {
			// if (success) {
			//  	ttt_time_file << (s_CPU_final - s_CPU_inicial) << endl;
			//  	ttt_cfo_file << func->getFnEvals() << endl;
			// }

			fo_file << result << endl;
		}

		//cout << "result: " << result << endl;
		
		
		soma += result;
		success_rate += success;
		avg_cfos += func->getFnEvals();
		avg_time += (s_CPU_final - s_CPU_inicial);
		delete func;
	}
	
	double average = soma / num_runs;
	avg_cfos /= num_runs;
	avg_time /= num_runs;
	//success_rate = success_rate * 100 / num_runs;
	avg_time /= num_runs;
	
	string filepath("out(" + arguments + ").csv");

	std::ofstream outfile;
	if (!exist(filepath.c_str())) {
		outfile.open(filepath);
		outfile << "Func;Dim;BestFO;AvgFO;CFOs;Time;SC" << endl;
	} else {
		outfile.open(filepath, std::ofstream::out | std::ofstream::app);
	}

	outfile << funcCode << ";" << problemDim << ";" << min << ";" << average << ";" << avg_cfos << ";" << avg_time << ";" << success_rate << endl;
	outfile.close();
	// ttt_time_file.close();
	// ttt_cfo_file.close();
	fo_file.close();

	//if (print) 
	//	cout << "Result: " << result << endl;

	return 0;	
}

