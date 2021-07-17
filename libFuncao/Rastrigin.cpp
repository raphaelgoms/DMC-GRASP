#include <stdio.h>
#include <math.h>

#include "Rastrigin.h"
#include "Util.h"

Rastrigin::Rastrigin(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Rastrigin::~Rastrigin(){

}

void Rastrigin::setFnEvals(int c){
	cont = c;	
}

int Rastrigin::getFnEvals(){
	return cont;	
}

double Rastrigin::getGap(){
	return (bestValue - minValue);
}

bool Rastrigin::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Rastrigin::calc(double *x){
	cont++;
	long double value = 0, pi=3.14159265;
	int i;

	value = 10*n;	for (i=0; i<n; i++)
		value += (pow(x[i],2) - 10*cos(2*pi*x[i]));

	return value;
}


double Rastrigin::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Rastrigin::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
