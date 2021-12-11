#include <stdio.h>
#include <math.h>

#include "Shekel.h"
#include "Util.h"

Shekel::Shekel(int n, int m){
	cont = 0;
	this->n = n;
	this->m = m;
	minValue = -1.0;
	switch(m){
		case 5: minValue = -10.15319538;
				break;
		case 7: minValue = -10.40281868;
				break;
		case 10: minValue = -10.53628349;
				break;
	}
}

Shekel::~Shekel(){

}

void Shekel::setFnEvals(int c){
	cont = c;	
}

int Shekel::getFnEvals(){
	return cont;	
}

double Shekel::getGap(){
	return (bestValue - minValue);
}

bool Shekel::isNearOptimum(double fBest){
	/*double bestValue = -1.0;
	switch(m){
		case 5: bestValue = -10.15319538;
				break;
		case 7: bestValue = -10.40281868;
				break;
		case 10: bestValue = -10.53628349;
				break;
	}*/

	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = fabs(minValue)*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation))){
		return true;
	} 

	return false;
}

double Shekel::calc(double *x){
	cont++;
	long double value,sum,totsum,temp;
  	long double a[10][4]={{4,4,4,4},{1,1,1,1},{8,8,8,8},{6,6,6,6},{7,3,7,3},{2,9,2,9},
                   {5,5,3,3},{8,1,8,1},{6,2,6,2},{7,3.6,7,3.6}};
  	long double c[10]={0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5};
  	int i,j;

  	totsum = 0;
  	for (i = 0; i < m; i++){
    	sum = 0;
    	for (j = 0; j < 4; j++){
    	  temp = (x[j] - a[i][j]) * (x[j] - a[i][j]);
    	  sum += temp;
    	}
    	sum += c[i];
    	totsum += pow(sum,-1.0);
  	}	

  	value = -totsum;
  	return value;

 }


 double Shekel::calc2(ap::real_1d_array x){
 	cont++;
	long double value,sum,totsum,temp;
  	long double a[10][4]={{4,4,4,4},{1,1,1,1},{8,8,8,8},{6,6,6,6},{3,7,3,7},{2,9,2,9},
                   {5,5,3,3},{8,1,8,1},{6,2,6,2},{7,3.6,7,3.6}};
  	long double c[10]={0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5};
  	int i,j;

  	totsum = 0;
  	for (i = 0; i < m; i++){
    	sum = 0;
    	for (j = 0; j < 4; j++){
    	  temp = (x(j+1) - a[i][j]) * (x(j+1) - a[i][j]);
    	  sum += temp;
    	}
    	sum += c[i];
    	totsum += pow(sum,-1);
  	}	

  	value = -totsum;
  	return value;
 }

 // Calculando o gradiente
 void Shekel::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	contGrad++;
	long double sum,totsum,temp;
  	long double a[10][4]={{4,4,4,4},{1,1,1,1},{8,8,8,8},{6,6,6,6},{3,7,3,7},{2,9,2,9},
                   {5,5,3,3},{8,1,8,1},{6,2,6,2},{7,2.6,7,3.6}};
  	long double c[10]={0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5};
  	int i,j;

	
	for (j = 0; j < 4; j++){
    	sum = 0;
    	for (i = 0; i < m; i++){
    		temp = -((2*(x(j+1) - a[i][j])) / (pow((x(j+1) - a[i][j]), 2) + c[i]));
    	    sum += temp;
    	}
    	g(j+1) = sum;
  	}	
 }
