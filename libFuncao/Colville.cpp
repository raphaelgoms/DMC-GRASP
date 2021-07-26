#include <stdio.h>
#include <math.h>

#include "Colville.h"
#include "Util.h"

Colville::Colville(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Colville::~Colville(){

}

void Colville::setFnEvals(int c){
	cont = c;	
}

int Colville::getFnEvals(){
	return cont;	
}

double Colville::getGap(){
	return (bestValue - minValue);
}

bool Colville::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Colville::calc(double *x){
	cont++;	
	long double value;

	value = 100*pow(x[1] - pow(x[0],2),2) + pow(1-x[0],2) + 90*pow(x[3]-pow(x[2],2),2) +
          pow(1-x[2],2) + 10.1*(pow(x[1] -1,2) + pow(x[3] -1,2)) + 19.8*(x[1] -1)*(x[3] -1);

	return value;
}


double Colville::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Colville::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
