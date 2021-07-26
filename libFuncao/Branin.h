#ifndef BRANIN_H_
#define BRANIN_H_

#include "Funcao.h"

class Branin: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Branin(int n);
		virtual ~Branin();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*BRANIN_H_*/
