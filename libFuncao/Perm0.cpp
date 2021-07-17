#include <stdio.h>
#include <math.h>

#include "Perm0.h"
#include "Util.h"

Perm_0::Perm_0(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Perm_0::~Perm_0(){

}

void Perm_0::setFnEvals(int c){
	cont = c;	
}

int Perm_0::getFnEvals(){
	return cont;	
}

double Perm_0::getGap(){
	return (bestValue - minValue);
}

bool Perm_0::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Perm_0::calc(double *x){
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


double Perm_0::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Perm_0::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
