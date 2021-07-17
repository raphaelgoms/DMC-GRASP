#ifndef RASTRIGIN_H_
#define RASTRIGIN_H_

#include "Funcao.h"

class Rastrigin: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Rastrigin(int n);
		virtual ~Rastrigin();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*RASTRIGIN_H_*/
