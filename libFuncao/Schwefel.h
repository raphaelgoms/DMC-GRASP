#ifndef SCHWEFEL_H_
#define SCHWEFEL_H_

#include "Funcao.h"

class Schwefel: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Schwefel(int n);
		virtual ~Schwefel();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /* SCHWEFEL_H_*/
