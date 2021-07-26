#include <stdio.h>
#include <math.h>
#include "SumSquares.h"

SumSquares::SumSquares(int n){
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

SumSquares::~SumSquares(){

}

void SumSquares::setFnEvals(int c){
	cont = c;	
}

int SumSquares::getFnEvals(){
	return cont;	
}

double SumSquares::getGap(){
	return (bestValue - minValue);
}


bool SumSquares::isNearOptimum(double fBest){
/*	int max = 1;
	int min = -1;
	
	if (fBest > 1.0){
		return false;
	}
	
	int fX = (int)(fBest*1000000.0);
	//printf("Fbest = %.10lf - %lf,  Greater = %d \n", fBest, fBest, fX);
	if ((fX >= min) && (fX <= max)){
		return true;
	}
	 
	return false;
*/
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = fabs(minValue)*0.0001 + 0.000001;	
	//equation = fabs(bestValue)*0.0001 + 0.01;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double SumSquares::calc(double *x){
	double som1 = 0.0;
	double fx = 0.0; 
	cont++;
	
	for (int i = 0; i < n; i++){		som1 += (i+1)*(x[i]*x[i]);
	}

	//printf("Som1 = %lf - Som2 = %lf \n", som1, som2);
	fx = som1;
	return fx;
}

double SumSquares::calc2(double *x){
	double som1 = 0.0;
	double fx = 0.0; 
	cont++;
	
	for (int i = 1; i <= n; i++){		som1 += i*(x[i]*x[i]);
	}

	//printf("Som1 = %lf - Som2 = %lf \n", som1, som2);
	fx = som1;
	return fx;
}

void SumSquares::calcGrad(double *x, double *g){
	
//	g[0] = 2*x[0] + 2*x[1]*x[1];
//	g[1] = 1*x[0]*x[0] + 4*x[1];
	g[1] = 2*x[1] + 2*x[2]*x[2];
	g[2] = 1*x[1]*x[1] + 4*x[2];
		
}
