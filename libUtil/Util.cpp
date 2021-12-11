#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Util.h"

 Util::Util(){
 }

 Util::~Util(){
 }

 //void Util::setMtRand(MTRand *mtRand){
 //	this->mtRand = mtRand;			
 //}
 
 double Util::dRand(){
	return ((double)rand() / ((double)(RAND_MAX)) );
	//return mtRand->rand();
 }

 void Util::copy(double *xAux, double *x, int n){
  	for (int i = 0; i < n; i++){
  		xAux[i] = x[i];
  	}
  }

 // Se x está no intervalo entre 'l' e 'u' retorne true, senao retorne false. 
 bool Util::feasible(double *x, double *l, double *u, int n){
	 for (int i = 0; i < n; i++){
		 if ((x[i] < l[i]) || (x[i] > u[i])){
			return false;
		 }
	 }
	 
	 return true;
 }

 bool Util::equals(double x1, double x2){
	if (fabs(x1 - x2) > 0.000001){
		return false;
	}

	return true;
 }	

 void Util::printX(double *x, int n){
	printf("(");
	for (int i = 0; i < n; i++){
		printf(" %.8lf", x[i]);
		if (i < n-1){
			printf(",");
		}
	}
	printf(")");
 }	 	

 double Util::calcNorma(double *x, double *xGrid, int n){
 	double aux = 0.0, norma;
 	
 	for (int i = 0; i < n; i++){
 		aux += pow((xGrid[i]-x[i]), 2);
 	}
 
 	norma = sqrt(aux);
 	return norma;

 }


 double Util::calcNorma(double *x, int n){
 	double aux = 0.0, norma;
 	
 	for (int i = 0; i < n; i++){
 		aux += pow(x[i], 2);
 	}
 
 	norma = sqrt(aux);
 	return norma;
 }

 double Util::dist(double *x1, double *x2, int n){
	double aux = 0.0;
		
	for (int i = 0; i < n; i++){
		aux += pow(x2[i]-x1[i], 2);	
	}

	return	sqrt(aux);
 }

 double Util::maxreal(double x1, double x2){
	return	((x1 > x2) ? x1: x2);
 }

 double Util::minreal(double x1, double x2){
	return	((x1 < x2) ? x1: x2);
 }

 double Util::dotproduct(double *x1, double *x2, int n){
	double aux = 0.0;
	for (int i = 0; i < n; i++){
		aux += x1[i]*x2[i];
	}
	return aux;			
 }
 
 void Util::addvector(double *x1, double *x2, int n){
	for (int i = 0; i < n; i++){
		x1[i] += x2[i];
	}
 } 

 void Util::subvector(double *x1, double *x2, int n){
	for (int i = 0; i < n; i++){
		x1[i] -= x2[i];
	}
 }

 double Util::sqr(double x){
	return x*x;
 }



