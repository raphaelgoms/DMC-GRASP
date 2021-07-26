#include <stdio.h>
#include <math.h>

#include "Shubert.h"
#include "Util.h"

Shubert::Shubert(int n){
	contGrad = 0;	
	cont = 0;
	this->n = n;
	minValue = -186.7309;
}

Shubert::~Shubert(){

}

void Shubert::setFnEvals(int c){
	cont = c;	
	contGrad = 0;
}

int Shubert::getFnEvals(){
	return cont;	
}

int Shubert::getGradEvals(){
	return contGrad;	
}

double Shubert::getGap(){
	return (bestValue - minValue);
}

bool Shubert::isNearOptimum(double fBest){
	//double bestValue = -186.7309;
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = fabs(minValue)*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Shubert::calc(double *x){
	cont++;
	long double value;

  	value = (cos(2*x[0] +1) + 2*cos(3*x[0] +2) + 3*cos(4*x[0] +3) +
    	       4*cos(5*x[0] +4) + 5*cos(6*x[0] +5)) * (cos(2*x[1] +1) + 2*cos(3*x[1] +2) +
    	       3*cos(4*x[1] +3) + 4*cos(5*x[1] +4) + 5*cos(6*x[1] +5));
  	return value;
}


double Shubert::calc2(ap::real_1d_array x){
	cont++;
	long double value;

  	value = (cos(2*x(1) +1) + 2*cos(3*x(1) +2) + 3*cos(4*x(1) +3) +
    	       4*cos(5*x(1) +4) + 5*cos(6*x(1) +5)) * (cos(2*x(2) +1) + 2*cos(3*x(2) +2) +
    	       3*cos(4*x(2) +3) + 4*cos(5*x(2) +4) + 5*cos(6*x(2) +5));
  	return value;
}


void Shubert::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	contGrad++;

  	g(1) = -2 *(cos(2*x(2) +1) + 2*cos(3*x(2) +2) + 3*cos(4*x(2) +3) +
    	       4*cos(5*x(2) +4) + 5*cos(6*x(2) +5)) * (sin(2*x(1) +1) + 3*sin(3*x(1) +2) +
    	       6*sin(4*x(1) +3) + 10*sin(5*x(1) +4) + 15*sin(6*x(1) +5));
  

  	g(2) = -2 *(cos(2*x(1) +1) + 2*cos(3*x(1) +2) + 3*cos(4*x(1) +3) +
    	       4*cos(5*x(1) +4) + 5*cos(6*x(1) +5)) * (sin(2*x(2) +1) + 3*sin(3*x(2) +2) +
    	       6*sin(4*x(2) +3) + 10*sin(5*x(2) +4) + 15*sin(6*x(2) +5));

}
