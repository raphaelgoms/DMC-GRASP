#ifndef FUNCAO_H_
#define FUNCAO_H_

#include "ap.h"

class Funcao
{
	protected:
		int cont;
		int contGrad;
		double minValue;
		double bestValue;

	public:
		static const int ROSENBROCK 	= 1;
		static const int ZAKHAROV		= 2;
		static const int SUMSQUARES		= 3;
		static const int BRANIN			= 4;
		static const int EASOM			= 5;
		static const int GOLDSTEINPRICE = 6;
		static const int SHEKEL 		= 7;
		static const int HARTMANN 		= 8;
		static const int SHUBERT 		= 9;
		static const int BEALE 			= 10;
		static const int BOOTH 			= 11;
		static const int BOHACHEVSKY 	= 12;
		static const int HUMP 			= 13;
		static const int MATYAS 		= 14;
		static const int SCHWEFEL 		= 15;
		static const int COLVILLE		= 16;
		static const int PERM			= 17;
		static const int PERM0			= 18;
		static const int POWERSUM		= 19;
		static const int GRIEWANK		= 20;
		static const int RASTRIGIN		= 21;
		static const int TRID			= 22;
		static const int POWELL			= 23;
		static const int DIXONPRICE		= 24;
		static const int ACKLEY			= 25;
		static const int LEVY			= 26;
		static const int SPHERE			= 27;

		Funcao();
		virtual ~Funcao();
		virtual double calc(double *x);
		virtual double calc2(ap::real_1d_array x);
		virtual void calcGrad(ap::real_1d_array &x, ap::real_1d_array &g);
		virtual bool isNearOptimum(double fBest);				
		
		virtual int getFnEvals();
		virtual int getGradEvals();
		virtual double getGap();

		void setBestValue(double fX);		
};

#endif /*FUNCAO_H_*/
