#include "Funcao.h"
#include <limits>

Funcao::Funcao(){
	cont = 0;
	contGrad = 0;	
	bestValue = std::numeric_limits<double>::max();
}

Funcao::~Funcao(){
}

int Funcao::getGradEvals(){
	return contGrad;	
}

int Funcao::getFnEvals(){
	return cont;	
}

double Funcao::getGap(){
	return std::numeric_limits<double>::max();
}

void Funcao::setBestValue(double fX){
	bestValue = fX;
}
		
bool Funcao::isNearOptimum(double fBest){
	return true;
}
		
double Funcao::calc(double *x){
	return 0.0;
}

double Funcao::calc2(ap::real_1d_array x){
	return 0.0;
}

void Funcao::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	return;
}


