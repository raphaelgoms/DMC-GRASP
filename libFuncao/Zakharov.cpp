#include <stdio.h>
#include <math.h>
#include "Zakharov.h"
#include "Util.h"


Zakharov::Zakharov(int n){
	contGrad = 0;	
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Zakharov::~Zakharov(){

}

void Zakharov::setFnEvals(int c){
	cont = c;	
	contGrad = 0;
}

int Zakharov::getFnEvals(){
	return cont;	
}

int Zakharov::getGradEvals(){
	return contGrad;	
}

double Zakharov::getGap(){
	return (bestValue - minValue);
}

bool Zakharov::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = fabs(bestValue)*0.0001 + 0.000001;	
//	equation = fabs(bestValue)*0.0001 + 0.01;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Zakharov::calc(double *x){
	double som1 = 0.0, som2 = 0.0;
	double fx = 0.0; 
	cont++;

	for (int i = 0; i < n; i++){		som1 += x[i]*x[i];
		som2 += 0.5*x[i]*(double(i+1));
	}

	//printf("Som1 = %lf - Som2 = %lf \n", som1, som2);
	fx = som1 + pow(som2, 2) + pow(som2,4);
	return fx;
}

double Zakharov::calc2(ap::real_1d_array x){
	double som1 = 0.0, som2 = 0.0;
	double fx = 0.0; 
	cont++;

	for (int i = 1; i <= n; i++){		som1 += x(i)*x(i);
		som2 += 0.5*x(i)*(double)i;
	}

	//printf("Som1 = %lf - Som2 = %lf \n", som1, som2);
	fx = som1 + pow(som2, 2) + pow(som2,4);
	return fx;
}


void Zakharov::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	double som = 0.0;
	contGrad++;
	
	for (int i = 1; i <= n; i++){		som += 0.5*x(i)*(double)i;
	}

	for (int i = 1; i < n; i++){
		//	g(i)= 2*x(i) + 2*som*(0.5*(double)i) + 4*pow(som,3)*(0.5*(double)i);
		g(i) = 2*x(i) + som*(double)i + 2*(double)i*pow(som,3);
	}
}
