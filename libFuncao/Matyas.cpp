#include <stdio.h>
#include <math.h>

#include "Matyas.h"
#include "Util.h"

Matyas::Matyas(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Matyas::~Matyas(){

}

void Matyas::setFnEvals(int c){
	cont = c;	
}

int Matyas::getFnEvals(){
	return cont;	
}

double Matyas::getGap(){
	return (bestValue - minValue);
}

bool Matyas::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Matyas::calc(double *x){
	cont++;	
	long double value = 0;

	value = 0.26*(pow(x[0],2) + pow(x[1],2)) - 0.48*x[0]*x[1];

	return value;
}


double Matyas::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Matyas::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
