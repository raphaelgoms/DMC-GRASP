#include <stdio.h>
#include <math.h>

#include "Branin.h"
#include "Util.h"

Branin::Branin(int n){
	cont = 0;
	this->n = n;
	minValue = 0.397887;
}

Branin::~Branin(){

}

void Branin::setFnEvals(int c){
	cont = c;	
}

int Branin::getFnEvals(){
	return cont;	
}

double Branin::getGap(){
	return (bestValue - minValue);
}

bool Branin::isNearOptimum(double fBest){
	//double bestValue = 0.397887;
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Branin::calc(double *x){
	cont++;
	long double fX,pi=3.14159265;

  	fX = pow(x[1] - (5.1*pow(x[0],2))/(4*pow(pi,2)) + 5*x[0]/pi - 6,2) +
          10*(1 - 1/(8*pi)) * cos(x[0]) + 10;
	
   return fX;
}


double Branin::calc2(ap::real_1d_array x){
	cont++;
	long double fX,pi=3.14159265;

  	fX = pow(x(2) - (5.1*pow(x(1),2))/(4*pow(pi,2)) + 5*x(1)/pi - 6,2) +
          10*(1 - 1/(8*pi)) * cos(x(1)) + 10;
	
   return fX;
}


void Branin::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	long double pi=3.14159265;

	// g(1) = 2((-5/4pi^2) + (5/pi))*((5x/pi) - (5x/4pi^2)+y-6) - 10(1-1/8pi)sin(x)
	// g(2) = 2((5x/pi)-(5x/4pi^2) + y - 6)

	contGrad++;
	g(1) = 2*(-5.0/(4.0*pow(pi,2)) + (5.0/pi))* ((5.0*x(1)/pi) - (5.0*x(1)/(4*pow(pi,2))) + x(2) - 6.0) - 
			10*(1 - 1/(8*pi))*sin(x(1));

	g(2) =  2*((5.0*x(1)/pi) - (5.0*x(1)/(4*pow(pi,2))) + x(2) - 6.0);

}
