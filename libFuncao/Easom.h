#ifndef EASOM_H_
#define EASOM_H_

#include "Funcao.h"

class Easom: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Easom(int n);
		virtual ~Easom();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

		virtual double getGap();
};

#endif /*EASOM_H_*/
