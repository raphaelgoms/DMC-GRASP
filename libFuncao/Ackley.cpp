#include <stdio.h>
#include <math.h>

#include "Ackley.h"
#include "Util.h"

Ackley::Ackley(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Ackley::~Ackley(){

}

void Ackley::setFnEvals(int c){
	cont = c;	
}

int Ackley::getFnEvals(){
	return cont;	
}

double Ackley::getGap(){
	return (bestValue - minValue);
}

bool Ackley::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Ackley::calc(double *x){
	cont++;
	long double value, sum1=0, sum2=0, pi=3.14159265;
	int i;
	long double value1, value2, value3;

	for (i = 0; i < n; i++)
		sum1 += pow(x[i], (long double) 2);

	for (i = 0; i < n; i++)
    	sum2 += cos(2 * pi * x[i]);
	
	value1 = -20*exp((long double) (-0.2* (long double) sqrt((long double) sum1/n)));
	value2 = - exp((long double) (sum2/n));
	value3 =  20 + exp((long double) 1);
	
	value = value1 + value2 + value3;

	return value;
}


double Ackley::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Ackley::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
