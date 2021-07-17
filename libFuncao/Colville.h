#ifndef COLVILLE_H_
#define COLVILLE_H_

#include "Funcao.h"

class Colville: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Colville(int n);
		virtual ~Colville();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*COLVILLE_H_*/
