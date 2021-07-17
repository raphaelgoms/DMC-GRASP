#ifndef UTIL_H_
#define UTIL_H_

#include "MersenneTwister.h"

class Util{
	private:
		MTRand *mtRand;
				
	public:
		Util();
		virtual ~Util();

		//static void setMtRand(MTRand *mtRand2);
		static double dRand();
		
		static void copy(double *xAux, double *x, int n);
		static double calcNorma(double *x, double *xGrid, int n);
		static double calcNorma(double *x, int n);
		static double dist(double *x1, double *x2, int n);		

		static bool feasible(double *x, double *l, double *u, int n);
		static bool equals(double x1, double x2);
		static void printX(double *x, int n);

		static double maxreal(double x1, double x2);
		static double minreal(double x1, double x2);
		static double dotproduct(double *x1, double *x2, int n);
		static double sqr(double x1);
		static void addvector(double *x1, double *x2, int n);
		static void subvector(double *x1, double *x2, int n);


};

#endif
