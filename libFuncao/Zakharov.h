#ifndef ZAKHAROV_H_
#define ZAKHAROV_H_

#include "Funcao.h"

class Zakharov: public Funcao{
	private:
		int n;
		int cont;
		int contGrad;
			
	public:
		Zakharov(int n);
		virtual ~Zakharov();
	

		virtual int getFnEvals();
		virtual int getGradEvals();	
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);
		virtual double getGap();
};

#endif /*ZAKHAROV_H_*/
