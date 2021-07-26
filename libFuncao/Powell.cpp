#include <stdio.h>
#include <math.h>

#include "Powell.h"
#include "Util.h"

Powell::Powell(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Powell::~Powell(){

}

void Powell::setFnEvals(int c){
	cont = c;	
}

int Powell::getFnEvals(){
	return cont;	
}

double Powell::getGap(){
	return (bestValue - minValue);
}

bool Powell::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Powell::calc(double *x){
	cont++;
	long double value;

	value = pow(x[0] + 10*x[1],2) + 5*pow(x[2]-x[3],2) + pow(x[1] - 2*x[2],4) + 10*pow(x[0]-x[3],4);
	return value;
}


double Powell::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Powell::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
