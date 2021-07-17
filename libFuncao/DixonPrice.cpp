#include <stdio.h>
#include <math.h>

#include "DixonPrice.h"
#include "Util.h"

DixonPrice::DixonPrice(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

DixonPrice::~DixonPrice(){

}

void DixonPrice::setFnEvals(int c){
	cont = c;	
}

int DixonPrice::getFnEvals(){
	return cont;	
}

double DixonPrice::getGap(){
	return (bestValue - minValue);
}

bool DixonPrice::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double DixonPrice::calc(double *x){
	cont++;  
	long double value=0;
	int i;

	value = pow(1 - x[0],2) + pow(1 - x[n-1],2);

	for (i = 1; i < n-1; i++)
		value += pow(pow(x[i],2) - x[i+1],2);

	return value;
}


double DixonPrice::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void DixonPrice::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
