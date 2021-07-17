#include <stdio.h>
#include <math.h>

#include "Hartmann.h"
#include "Util.h"

Hartmann::Hartmann(int n, int m){
	cont = 0;
	this->n = n;
	this->m = m;
	switch(n){
		case 3: minValue = -3.86278;
				break;
		case 6: minValue = -3.32237;
				break;
	}
}

Hartmann::~Hartmann(){

}

void Hartmann::setFnEvals(int c){
	cont = c;	
}

int Hartmann::getFnEvals(){
	return cont;	
}

double Hartmann::getGap(){
	return (bestValue - minValue);
}

bool Hartmann::isNearOptimum(double fBest){
	double deltaValue =	fabs(fBest - minValue);
	double equation;
	
	equation = fabs(minValue)*0.0001 + 0.000001;	
	if ((deltaValue < equation) || (Util::equals(deltaValue, equation)) /* || fBest < minValue */){
		return true;
	} 

	return false;
}

 double Hartmann::calc(double *x){
	cont++;
	long double value;
  	
	switch (n){
		case 3: value = func43(x);
				break;
		case 6: value = func46(x);
				break;	
	}

	return value;
 }


 // Hartmann function (3,4)
 double Hartmann::func43(double *x){
 	 long double value = 0;
 	 int i, j;
 	 long double qi;

 	 long double A[4][3] = {{3.0,10.0,30.0},{0.1,10.0,35.0},{3.0,10.0,30.0},{0.1,10.0,35.0}};
 	 long double P[4][3] = {{0.6890,0.1170,0.2673},{0.4699,0.4387,0.7470},
				 {0.1091,0.8732,0.5547},{0.0381,0.5743,0.8828}};
 	 long double alpha[4] = {1.0,1.2,3.0,3.2};

 	 for (i=0;i<4;i++){
 	     qi = 0;
 	     for(j=0;j<3;j++)
			qi += A[i][j]*pow(x[j] - P[i][j],2);

    	 value += -1*alpha[i]*exp(-1*qi);
    }

  	return value; 
 }


 // Hartmann function (6,4)
 double Hartmann::func46(double *x){
	 double value = 0;
	 int i, j;
	 double qi;
 
	 double B[4][6] = {{10,3,17,3.5,1.7,8},
				 {0.05,10,17,0.1,8,14},
				 {3,3.5,1.7,10,17,8},
				 {17,8,0.05,10,0.1,14}};
	  double Q[4][6] = {{0.1312,0.1696,0.5569,0.0124,0.8283,0.5886},
				 {0.2329,0.4135,0.8307,0.3736,0.1004,0.9991},
				 {0.2348,0.1451,0.3522,0.2883,0.3047,0.6650},
				 {0.4047,0.8828,0.8732,0.5743,0.1091,0.0381}};
	  double alpha[4] = {1,1.2,3,3.2};
	

	  	for (i=0;i<4;i++){
	      	qi = 0.0;
    	  	for(j=0;j<6;j++)
				qi += B[i][j]*pow(x[j] - Q[i][j],2);
	
    	  	value += -1*alpha[i]*exp(-1*qi);
    	}

		//if (value < -3.32237)	
		//	printf("Value = %lf \n", value);
  		return value;
 	}


 double Hartmann::calc2(ap::real_1d_array x){
 	 long double value = 0;
 	 int i, j;
 	 long double qi;

	 cont++;
 	 long double A[4][3] = {{3.0,10.0,30.0},{0.1,10.0,35.0},{3.0,10.0,30.0},{0.1,10.0,35.0}};
 	 long double P[4][3] = {{0.6890,0.1170,0.2673},{0.4699,0.4387,0.7470},
				 {0.1091,0.8732,0.5547},{0.0381,0.5743,0.8828}};

  	long double B[4][6] = {{10.0,3.0,17.0,3.5,1.7,8.0},
				 {0.05,10.0,17.0,0.1,8.0,14.0},
				 {3.0,3.5,1.7,10.0,17.0,8.0},
				 {17.0,8.0,0.05,10.0,0.1,14.0}};
  	long double Q[4][6] = {{0.1312,0.1696,0.5569,0.0124,0.8283,0.5886},
				 {0.2329,0.4135,0.8307,0.3736,0.1004,0.9991},
				 {0.2348,0.1451,0.3522,0.2883,0.3047,0.6650},
				 {0.4047,0.8828,0.8732,0.5743,0.1091,0.0381}};
 	 long double alpha[4] = {1.0,1.2,3.0,3.2};

	 if (n == 3){
 		
		 for (i=0;i<4;i++){
 		     qi = 0;
 		     for(j=0;j<3;j++)
				qi += A[i][j]*pow(x(j+1) - P[i][j],2);
	
    		 value += -1*alpha[i]*exp(-1*qi);
    	}
	}
	else if (n==6){
	  	for (i=0;i<4;i++){
	      	qi = 0;
	      	for(j=0;j<6;j++)
				qi += B[i][j]*pow(x(j+1) - Q[i][j],2);

	      	value += -1*alpha[i]*exp(-1*qi);
	    }
	}
  	

	return value;
 }
 
 void Hartmann::calcGrad(ap::real_1d_array &x, ap::real_1d_array &g){
	 long double value = 0;
 	 int i, j, k;
 	 long double qi, qiAux;

 	 long double A[4][3] = {{3.0,10.0,30.0},{0.1,10.0,35.0},{3.0,10.0,30.0},{0.1,10.0,35.0}};
 	 long double P[4][3] = {{0.6890,0.1170,0.2673},{0.4699,0.4387,0.7470},
				 {0.1091,0.8732,0.5547},{0.0381,0.5743,0.8828}};
  	long double B[4][6] = {{10.0,3.0,17.0,3.5,1.7,8.0},
				 {0.05,10.0,17.0,0.1,8.0,14.0},
				 {3.0,3.5,1.7,10.0,17.0,8.0},
				 {17.0,8.0,0.05,10.0,0.1,14.0}};
  	long double Q[4][6] = {{0.1312,0.1696,0.5569,0.0124,0.8283,0.5886},
				 {0.2329,0.4135,0.8307,0.3736,0.1004,0.9991},
				 {0.2348,0.1451,0.3522,0.2883,0.3047,0.6650},
				 {0.4047,0.8828,0.8732,0.5743,0.1091,0.0381}};


 	 long double alpha[4] = {1.0,1.2,3.0,3.2};
	
	contGrad++;
	
	// Hartmann(3,4).	
	if (n == 3){

		for (k = 1; k <= 3; k++){
			value = 0.0;		
			for (i=0;i<4;i++){
 			     qi = 0;
				 
				for(j=0;j<3;j++)
					qi += A[i][j]*pow(x(j+1) - P[i][j],2);
	
				 qiAux =  2.0*A[i][k]*x(k) - P[i][k];		     
    			 value += -1*alpha[i]*qiAux*exp(-1*qi);
			}
			g(k) = value;
		}
	}
	// Hartmann(6,4).	
	else if (n == 6){

		for (k = 1; k <= 6; k++){
			value = 0.0;		
			for (i=0;i<4;i++){
 			     qi = 0;
				 
				for(j=0;j<6;j++)
					qi += B[i][j]*pow(x(j+1) - Q[i][j],2);
	
				 qiAux =  2.0*B[i][k]*x(k) - Q[i][k];		     
    			 value += -1*alpha[i]*qiAux*exp(-1*qi);
			}
			g(k) = value;
		}
	}

 }

