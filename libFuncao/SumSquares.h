#ifndef SUMSQUARES_H_
#define SUMSQUARES_H_

#include "Funcao.h"
#include "Util.h"

class SumSquares: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		SumSquares(int n);
		virtual ~SumSquares();

		void setFnEvals(int c);
		virtual int getFnEvals();
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(double *x);
		virtual void calcGrad(double *x, double *g);
		virtual double getGap();
};

#endif /*SumSquares_H_*/
