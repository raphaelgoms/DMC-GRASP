#ifndef SHUBERT_H_
#define SHUBERT_H_

#include "Funcao.h"

class Shubert: public Funcao{
	private:
		int n;
		int cont;
		int contGrad;
		
	public:
		Shubert(int n);
		virtual ~Shubert();

		void setFnEvals(int c);
		virtual int getFnEvals();
		virtual int getGradEvals();
		virtual double getGap();

		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);
};

#endif /*Shubert_H_*/
