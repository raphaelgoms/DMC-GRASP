#ifndef HARTMANN_H_
#define HARTMANN_H_

#include "Funcao.h"

class Hartmann: public Funcao{
	private:
		int m, n;
		int cont;
		
		double func43(double *x);
		double func46(double *x);
		
	public:
		Hartmann(int n, int m);
		virtual ~Hartmann();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);
		virtual double getGap();

};

#endif /*Hartmann_H_*/
