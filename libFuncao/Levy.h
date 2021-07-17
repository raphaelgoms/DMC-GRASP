#ifndef LEVY_H_
#define LEVY_H_

#include "Funcao.h"

class Levy: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Levy(int n);
		virtual ~Levy();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*LEVY_H_*/
