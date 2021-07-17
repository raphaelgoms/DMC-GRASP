#ifndef BOOTH_H_
#define BOOTH_H_

#include "Funcao.h"

class Booth: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Booth(int n);
		virtual ~Booth();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*BOOTH_H_*/
