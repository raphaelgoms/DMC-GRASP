#include <stdio.h>
#include <math.h>

#include "Trid.h"
#include "Util.h"

Trid::Trid(int n){
	cont = 0;
	this->n = n;
	switch (n){
		case 6: 	minValue = -50.0;
					break;
		case 10:	minValue = -210.0;
					break;
	}
}

Trid::~Trid(){

}

void Trid::setFnEvals(int c){
	cont = c;	
}

int Trid::getFnEvals(){
	return cont;	
}

double Trid::getGap(){
	return (bestValue - minValue);
}

bool Trid::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = minValue*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Trid::calc(double *x){
	cont++;
	long double value = 0;
	int i;

	for (i = 0; i < n; i++)
		value += pow(x[i]-1,2);

	for (i = 1; i < n; i++)
		value += -1*x[i]*x[i-1];

	return value;
}


double Trid::calc2(ap::real_1d_array x){
	double y = 0.0;	
	return y;

}


void Trid::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){

}
