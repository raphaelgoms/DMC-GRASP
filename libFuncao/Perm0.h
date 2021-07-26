#ifndef PERM_0_H_
#define PERM_0_H_

#include "Funcao.h"

class Perm_0: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Perm_0(int n);
		virtual ~Perm_0();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*PERM_0_H_*/
