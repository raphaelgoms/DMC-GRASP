#include <stdio.h>
#include <math.h>
#include "Rosenbrock2.h"
#include "Util.h"

Rosenbrock2::Rosenbrock2(int n){
	contGrad = 0;	
	cont = 0;
	this->n = n;
	minValue = 0.0;
}

Rosenbrock2::~Rosenbrock2(){

}

void Rosenbrock2::setFnEvals(int c){
	cont = c;	
	contGrad = 0;
}

int Rosenbrock2::getFnEvals(){
	return cont;	
}

int Rosenbrock2::getGradEvals(){
	return contGrad;	
}

double Rosenbrock2::getGap(){
	return (bestValue - minValue);
}

bool Rosenbrock2::isNearOptimum(double fBest){
	//double minValue = 0.0;
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = fabs(minValue)*0.0001 + 0.000001;	
	//equation = fabs(bestValue)*0.0001 + 0.01;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;

}

double Rosenbrock2::calc(double *x){
	double parent1, parent2;
	double fx = 0.0; 
	cont++;
		
	for (int i = 0; i < n-1; i++){		parent1 = pow(((x[i]*x[i]) - x[i+1]), 2); //(x[0]^2 - x[1])^2
		parent2 = pow((x[i] - 1), 2);           //(x[0] - 1)^2
	
		fx += (100*parent1 + parent2); 
	}

	//printf("Fx = %lf \n", fx);

	return fx;
}

double Rosenbrock2::calc2(ap::real_1d_array x){
	double parent1, parent2;
	double fx = 0.0; 
	cont++;
	for (int i = 1; i < n; i++){		parent1 = pow(((x(i)*x(i))-x(i+1)), 2); //(x[0]^2 - x[1])^2
		parent2 = pow((x(i) - 1), 2);           //(x[0] - 1)^2
	
		fx += (100*parent1 + parent2); 
	}
	
	return fx;
}

void Rosenbrock2::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	contGrad++;
	for (int i = 1; i < n; i++){
		g(i) = -400*x(i)*(x(i+1)-(x(i)*x(i))) -2*(1-x(i)); 
	 }
	
	g(n) = 200*(x(n)-(x(n-1)*x(n-1)));
}



