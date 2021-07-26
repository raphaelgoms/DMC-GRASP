#include <stdio.h>
#include <math.h>

#include "Booth.h"
#include "Util.h"

Booth::Booth(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Booth::~Booth(){

}

void Booth::setFnEvals(int c){
	cont = c;	
}

int Booth::getFnEvals(){
	return cont;	
}

double Booth::getGap(){
	return (bestValue - minValue);
}

bool Booth::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Booth::calc(double *x){
	cont++; 
	long double value = 0;

 	value = pow(x[0] + 2*x[1] - 7,2) + pow(2*x[0] + x[1] - 5,2);

	return value;
}


double Booth::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Booth::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
