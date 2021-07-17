#ifndef SHEKEL_H_
#define SHEKEL_H_

#include "Funcao.h"

class Shekel: public Funcao{
	private:
		int m, n;
		int cont;
		
	public:
		Shekel(int n, int m);
		virtual ~Shekel();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);
		virtual double getGap();

};

#endif /*SHEKEL_H_*/
