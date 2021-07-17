#include <stdio.h>
#include <math.h>

#include "Perm.h"
#include "Util.h"

Perm::Perm(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Perm::~Perm(){

}

void Perm::setFnEvals(int c){
	cont = c;	
}

int Perm::getFnEvals(){
	return cont;	
}

double Perm::getGap(){
	return (bestValue - minValue);
}

bool Perm::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Perm::calc(double *x){
	cont++;	
	long double value = 0, beta = 0.5, sum = 0;
	int k, i;

	for (k=0;k<4;k++){
		sum = 0;
		for (i=0;i<4;i++)
			sum += (pow((long double) i+1,(long double) k+1) + beta)*(pow(x[i]/(i+1),k+1) - 1);

	      value += pow(sum,2);
	}

	return value;
}


double Perm::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Perm::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
