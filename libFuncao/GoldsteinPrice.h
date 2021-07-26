#ifndef GOLDSTEINPRICE_H_
#define GOLDSTEINPRICE_H_

#include "Funcao.h"

class GoldsteinPrice: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		GoldsteinPrice(int n);
		virtual ~GoldsteinPrice();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);
		virtual double getGap();
};

#endif /*GOLDSTEINPRICE_H_*/
