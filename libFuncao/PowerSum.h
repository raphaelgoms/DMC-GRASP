#ifndef POWER_SUM_H_
#define POWER_SUM_H_

#include "Funcao.h"

class PowerSum: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		PowerSum(int n);
		virtual ~PowerSum();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*POWER_SUM_H_*/
