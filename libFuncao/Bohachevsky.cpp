#include <stdio.h>
#include <math.h>

#include "Bohachevsky.h"
#include "Util.h"

Bohachevsky::Bohachevsky(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Bohachevsky::~Bohachevsky(){
	
}

void Bohachevsky::setFnEvals(int c){
	cont = c;	
}

int Bohachevsky::getFnEvals(){
	return cont;	
}

double Bohachevsky::getGap(){
	return (bestValue - minValue);
}

bool Bohachevsky::isNearOptimum(double fBest){
	cont++;
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Bohachevsky::calc(double *x){
	cont++;  
	long double value = 0, pi=3.14159265;

	value = pow(x[0],2) + 2.0*pow(x[1],2) - 0.3*cos(3*pi*x[0]) - 0.4*cos(4*pi*x[1]) + 0.7;

	return value;
}


double Bohachevsky::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Bohachevsky::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
