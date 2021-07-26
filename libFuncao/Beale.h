#ifndef BEALE_H_
#define BEALE_H_

#include "Funcao.h"

class Beale: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Beale(int n);
		virtual ~Beale();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*BEALE_H_*/
