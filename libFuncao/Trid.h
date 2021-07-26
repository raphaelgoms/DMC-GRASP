#ifndef TRID_H_
#define TRID_H_

#include "Funcao.h"

class Trid: public Funcao{
	private:
		int n;
		int cont;
		
	public:
		Trid(int n);
		virtual ~Trid();
		virtual int getFnEvals();
		void setFnEvals(int c);
		virtual bool isNearOptimum(double fBest);				
		virtual double calc(double *x);
		virtual double getGap();
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);

};

#endif /*TRID_H_*/
