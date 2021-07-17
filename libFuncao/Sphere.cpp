#include <stdio.h>
#include <math.h>

#include "Sphere.h"
#include "Util.h"

Sphere::Sphere(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Sphere::~Sphere(){

}

void Sphere::setFnEvals(int c){
	cont = c;	
}

int Sphere::getFnEvals(){
	return cont;	
}

double Sphere::getGap(){
	return (bestValue - minValue);
}

bool Sphere::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Sphere::calc(double *x){
	cont++;  
	long double value,sum=0;
	int i;

	for (i = 0; i < n; i++)
		sum += pow(x[i],(long double) 2);

	value = sum;

	return value;
}


double Sphere::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Sphere::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
