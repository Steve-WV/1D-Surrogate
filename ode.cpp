// Comment: 

// The following code uses the open source C++-package on cubic spline interpolation included in the class
// spline.h which is made available on: https://kluge.in-chemnitz.de/opensource/spline/

/*
Copyright (C) 2022 Steve Wolff-Vorbeck (steve.wolff-vorbeck at mathematik.uni-freibrg.de)
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <ctime>        // std::time
#include <cstdlib>     
#include <algorithm> 
#include <cmath> 
#include <cstdio>
#include <vector>
#include<math.h>
#include<stdio.h>
#include<cctype>
using namespace std;

#include "ode.hpp"
#include "spline.h"
#include "boost/tokenizer.hpp"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Generators for the distributions 
std::default_random_engine generator_norm1;
std::default_random_engine generator_norm2;

// generator for standard normal distribution
 std::normal_distribution<double> distribution_norm1(0.0,1.0);

const double E=80.05863;
const double Schub=30.8;
const double pi = 3.141592653589;

const double a_sum = -1.5;
const double k0=0.0;
//const double kl=0.5;
const double kl=0.0;

inline double sqr(double s){ 

	return s*s;
}

inline double cub(double s){ 

	return s*s*s;
}
	
void ODESystem::CholeskyLoop () {

std::vector<T> vals_ktl1, vals_ktl2;
double x1,y1,z1;
double x2,y2,z2;
KTL_new1=Mat(L,L);
KTL_new2=Mat(L,L);
double b1= L;

	for (int i=0; i < L; ++i) {

			double b=i;	
			z1=(b/b1);
					
				for (int j=0;j < L; ++j) {
				
					double c=j;
					z2= c/b1;
					
					double ktl1 = exp(-(abs(z1-z2)/l_c1));
					double ktl2 = exp(-(abs(z1-z2)/l_c2));
				
					KTL_new1(i,j)=ktl1;
					KTL_new2(i,j)=ktl2;

					}											

				}
		
}

// Perform the geometric perturbations
void ODESystem::Perform_KTL() {

		std::vector<double> X,Y,Knots;
		//std::cerr<< "Performing Geometric perturbations..";
		Vec KTL_Chol1;
		Vec KTL_Chol2;
		
		KTL_Chol1= sigma_new1*KTL_Cholnew1*C_norm1; 
		KTL_Chol2 = sigma_new2*KTL_Cholnew1*C_norm2; 
	 	
		double b2=(L-1);
		double b12=(L-2);
		double b22= L;
		double b1=(L-1);

		for (int i=0; i<L; ++i){
				double b3=i+1;
				X.push_back(KTL_Chol1[i]);
				Y.push_back(KTL_Chol2[i]);
				Knots.push_back((b3/b22));
				
			}
	
		X1 = Vec::Zero(L-2);
	 	X2 = Vec::Zero(L-2);
	 	
		dX1 = Vec::Zero(L-2);
	 	dX2 = Vec::Zero(L-2);
	 	dXM1 = Vec::Zero(L-1);
	 	dXM2 = Vec::Zero(L-1);
	 	dMU = Vec::Zero(L-1);
	 	
  // Perform cubic spline interpolation (also deriviative); currently it is required that X is already sorted	 	
	tk::spline s1(Knots,X);
	tk::spline s2(Knots,Y);

 	double msz=m_size;
 	double l =L;
 	double scal = l/msz; 
 	double scal2 = 1/msz;
    
 	dXM1[L-2]+=scal2*s1.deriv(1,((b1-0.5)/b1));
 	dXM2[L-2]+=scal2*s2.deriv(1,((b1-0.5)/b1));
	dMU[L-2]+= 1.0;
 		
	for ( int i = 0; i < L-2; ++i){
			double b3 = i;
			int count= int (b3/scal);
			
			dXM1[i]+=scal2*s1.deriv(1,((b3+0.5)/b1));
    			dXM2[i]+=scal2*s2.deriv(1,((b3+0.5)/b1));
    			dMU[i]+=1.0;
    	
 		}

}
	
// Generate random numbers V^i for computing density variation and geometric perturbation
// via L*V^i using Cholesk decomposition	
	void ODESystem::Prep_Norm() {
	
			double trans1, trans2,trans3;	
			
			C_norm1=Vec::Zero(L);
			C_norm2=Vec::Zero(L);
			C_norm3=Vec::Zero(L);

			for (int i=0; i < L; ++i) {

					// numbers for geometric perturbation
					trans1 = distribution_norm1(generator_norm1);
					trans2 = distribution_norm1(generator_norm2);
					
					C_norm1[i]+= trans1; 
					C_norm2[i]+= trans2;
				
		}
}

void ODESystem::KLoop_det ( SpMat* th_k1) {
	
	std::vector<T> vals_k1;
	double b1= (L-1);

	for (int i=0; i<L-2; ++i) {
										
					vals_k1.push_back( T(i, i+1, 2*b1) );
					vals_k1.push_back( T(i, i, -b1) );
					vals_k1.push_back( T(i, i+2, -b1) );
										
				   }

(*th_k1).setFromTriplets( vals_k1.begin(), vals_k1.end() );

}


void ODESystem::PrepK_det() {
	std::vector<SpMat> th_k1;
	
	for ( int j = 0; j < numThreads; ++j ) {
		th_k1.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		
		//th_k1.push_back( SpMat(3*(L-2), 3*(L-2)) );
		//th_k2.push_back( SpMat(3*(L-2), 3*(L-2)) );
		//th_k3.push_back( SpMat(3*(L-2), 3*(L-2)) );
	}
  
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) {
		threads.push_back( std::thread( &ODESystem::KLoop_det, this, 
	&th_k1[j])  );
	}

	// join threads
	for (auto &thread : threads) thread.join();
	
	
	K6 = SpMat(4*(L-2)+3,4*(L-2)+3);
	
	//K1 = SpMat(3*(L-2),3*(L-2));
	//K2 = SpMat(3*(L-2),3*(L-2));
	//K3 = SpMat(3*(L-2),3*(L-2));
	for (int j = 0; j < numThreads; ++j) {
		
		K6 += th_k1[j];
	}

}

void ODESystem::KLoop ( SpMat* th_k1, SpMat* th_k2, SpMat* th_k3,SpMat* th_k4,SpMat* th_k5) {
	
	std::vector<T> vals_k1,vals_k2,vals_k3,vals_k4,vals_k5,vals_k6;
	double b1= (L-1);

	for (int i=1; i<L-3; ++i) {
				double k1=0.0;
				double k2=0.0;
				double k3=0.0;
				double k11=0.0;
				double k12=0.0;
				double k13=0.0;
	
				k1+= (dMU[i]*b1+dMU[i+1]*b1);
				k2-=dMU[i]*b1;
				k3-=dMU[i+1]*b1;
				
				k11+= (dMU[i]*b1+dMU[i+1]*b1);
				k12-=dMU[i]*b1;
				k13-=dMU[i+1]*b1;
			
				vals_k1.push_back( T(i, i, k1) );
				vals_k1.push_back( T(i, i-1, k2) );
				vals_k1.push_back( T(i, i+1, k3) );
				
				vals_k2.push_back( T(i+(L-2), i+(L-2), k1) );
				vals_k2.push_back( T(i+(L-2), i+1+(L-2), k3) );
			        vals_k2.push_back( T(i+(L-2), i-1+(L-2), k2) );
			        
				vals_k3.push_back( T(i+2*(L-2), i+2*(L-2), k1) );
				vals_k3.push_back( T(i+2*(L-2), i+1+2*(L-2), k3) );
				vals_k3.push_back( T(i+2*(L-2), i-1+2*(L-2), k2) );
				
				// mass constraint part
				vals_k4.push_back(T(4*(L-2), i+(L-2), (1/b1)));
				vals_k4.push_back(T(4*(L-2)+1, i+2*(L-2), (1/b1)));
				vals_k4.push_back(T(4*(L-2)+2, i+3*(L-2), (1/b1)));
				vals_k4.push_back(T(i+(L-2),4*(L-2), (1/b1)));
				vals_k4.push_back(T(i+2*(L-2),4*(L-2)+1, (1/b1)));
				vals_k4.push_back(T(i+3*(L-2),4*(L-2)+2, (1/b1)));
				
				vals_k5.push_back( T(i+3*(L-2), i+3*(L-2), k1) );
				vals_k5.push_back( T(i+3*(L-2), i+1+3*(L-2), k3) );
				vals_k5.push_back( T(i+3*(L-2), i-1+3*(L-2), k2) );
							
				vals_k6.push_back( T(i+3*(L-2), i+3*(L-2), 2*b1) );
				vals_k6.push_back( T(i+3*(L-2), i+1+3*(L-2), -b1) );
				vals_k6.push_back( T(i+3*(L-2), i-1+3*(L-2), -b1) );
						
						}
				
				
		
	vals_k1.push_back( T(0, 0, dMU[0]*b1+dMU[1]*b1 ) );
	vals_k1.push_back( T(0, 1, -dMU[1]*b1) );
	vals_k1.push_back( T(L-3, L-3, dMU[L-3]*b1+dMU[L-2]*b1) );
	vals_k1.push_back( T(L-3, L-4, -dMU[L-3]*b1) );
	
	vals_k2.push_back( T((L-2),(L-2), dMU[0]*b1+dMU[1]*b1) );
	vals_k2.push_back( T((L-2),1+(L-2), -dMU[1]*b1) );
	vals_k2.push_back( T((L-3)+(L-2), L-3+(L-2), dMU[L-3]*b1+dMU[L-2]*b1) );
	vals_k2.push_back( T((L-3)+(L-2), L-4+(L-2), -dMU[L-3]*b1) );
	
	vals_k3.push_back( T(2*(L-2), 2*(L-2), dMU[0]*b1+dMU[1]*b1) );
	vals_k3.push_back( T(2*(L-2),1+2*(L-2), -dMU[1]*b1) );
	vals_k3.push_back( T((L-3)+2*(L-2), (L-3)+2*(L-2), dMU[L-3]*b1+dMU[L-2]*b1) );
	vals_k3.push_back( T((L-3)+2*(L-2), L-4+2*(L-2), -dMU[L-3]*b1) );
	
	vals_k5.push_back( T(3*(L-2), 3*(L-2), dMU[0]*b1+dMU[1]*b1) );
	vals_k5.push_back( T(3*(L-2),1+3*(L-2),  -dMU[1]*b1) );
	vals_k5.push_back( T((L-3)+3*(L-2), (L-3)+3*(L-2),  dMU[L-3]*b1+dMU[L-2]*b1) );
	vals_k5.push_back( T((L-3)+3*(L-2), L-4+3*(L-2),  -dMU[L-3]*b1) );
	
	//mass constraint 
	vals_k4.push_back(T(4*(L-2), (L-2), (1/b1)));
	vals_k4.push_back(T(4*(L-2), 2*(L-2)-1, (1/b1)));
	vals_k4.push_back(T(4*(L-2)+1, 2*(L-2), (1/b1)));
	vals_k4.push_back(T(4*(L-2)+1, 3*(L-2)-1, (1/b1)));
	vals_k4.push_back(T(4*(L-2)+2, 3*(L-2), (1/b1)));
	vals_k4.push_back(T(4*(L-2)+2, 4*(L-2)-1, (1/b1)));
	
	vals_k4.push_back(T((L-2),4*(L-2), (1/b1)));
	vals_k4.push_back(T(2*(L-2)-1,4*(L-2), (1/b1)));
	vals_k4.push_back(T(2*(L-2),4*(L-2)+1, (1/b1)));
	vals_k4.push_back(T(3*(L-2)-1,4*(L-2)+1, (1/b1)));
	vals_k4.push_back(T(3*(L-2),4*(L-2)+2, (1/b1)));
	vals_k4.push_back(T(4*(L-2)-1,4*(L-2)+2, (1/b1)));	
	
	(*th_k1).setFromTriplets( vals_k1.begin(), vals_k1.end() );
	(*th_k2).setFromTriplets( vals_k2.begin(), vals_k2.end() );
	(*th_k3).setFromTriplets( vals_k3.begin(), vals_k3.end() );
	(*th_k4).setFromTriplets( vals_k4.begin(), vals_k4.end() );
	(*th_k5).setFromTriplets( vals_k5.begin(), vals_k5.end() );
	
}

void ODESystem::PrepK() {
	std::vector<SpMat> th_k1,th_k2,th_k3,th_k4,th_k5;
	
	for ( int j = 0; j < numThreads; ++j ) {
		th_k1.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		th_k2.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		th_k3.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		th_k4.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		th_k5.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		
		
		//th_k1.push_back( SpMat(3*(L-2), 3*(L-2)) );
		//th_k2.push_back( SpMat(3*(L-2), 3*(L-2)) );
		//th_k3.push_back( SpMat(3*(L-2), 3*(L-2)) );
	}
  
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) {
		threads.push_back( std::thread( &ODESystem::KLoop, this, 
	&th_k1[j],&th_k2[j],&th_k3[j],&th_k4[j],&th_k5[j])  );
	}

	// join threads
	for (auto &thread : threads) thread.join();
	
	K1 = SpMat(4*(L-2)+3,4*(L-2)+3);
	K2 = SpMat(4*(L-2)+3,4*(L-2)+3);
	K3 = SpMat(4*(L-2)+3,4*(L-2)+3);
	K5 = SpMat(4*(L-2)+3,4*(L-2)+3);
	L1 = SpMat(4*(L-2)+3,4*(L-2)+3);
	
	//K1 = SpMat(3*(L-2),3*(L-2));
	//K2 = SpMat(3*(L-2),3*(L-2));
	//K3 = SpMat(3*(L-2),3*(L-2));
	for (int j = 0; j < numThreads; ++j) {
		K1 += th_k1[j];
		K2 += th_k2[j];
		K3 += th_k3[j];
		K5 += th_k5[j];
		L1 += th_k4[j];
	}
}

void ODESystem::PerturbeLoop(SpMat* th_k, Vec* th_b) { 

std::vector<T> vals_p;
*th_b = Vec::Zero(4*(L-2)+3);
double b1= (L-1);

	for (int i=1; i<L-3; ++i) {
				
				double k11x= 0.0;
				double k12x= 0.0;
				double k13x= 0.0;
				
				double k11y= 0.0;
				double k12y= 0.0;
				double k13y= 0.0;
				double kx1=0.0;
				double kx2=0.0;
				double kx3=0.0;
				double ky1=0.0;
				double ky2=0.0;
				double ky3=0.0;
				
				double kxy=0.0;
				double kxy1=0.0;
				double kxy2=0.0;
				double kxx=0.0;
				double kxx1=0.0;
				double kxx2=0.0;
				double kyy=0.0;
				double kyy1=0.0;
				double kyy2=0.0;
				double kx=0.0;
				double ky =0.0;
				double kp=0.0;
				
				// u' in perturbe matrix
				k11x+=0.5*(dMU[i]*dXM2[i]-dMU[i+1]*dXM2[i+1]);
				k12x+=0.5*dMU[i]*dXM2[i];
				k13x-=0.5*dMU[i+1]*dXM2[i+1];
				
				k11y+=0.5*(dMU[i]*dXM1[i]-dMU[i+1]*dXM1[i+1]);
				k12y+=0.5*dMU[i]*dXM1[i];
				k13y-=0.5*dMU[i+1]*dXM1[i+1];

				vals_p.push_back( T(i+(L-2), i, -k11x) );
				vals_p.push_back( T(i+(L-2), i-1, -k12x) );
				vals_p.push_back( T(i+(L-2), i+1, -k13x) );
							
				vals_p.push_back( T(i+2*(L-2), i, k11y) );
				vals_p.push_back( T(i+2*(L-2), i-1, k12y) );
				vals_p.push_back( T(i+2*(L-2), i+1, k13y) );

				// mass matrix part
				kx1= 0.5*(dMU[i]*dXM2[i]-dMU[i+1]*dXM2[i+1]);
				kx2= 0.5*(-dMU[i]*dXM2[i]);
				kx3= 0.5*(dMU[i+1]*dXM2[i+1]);
				
				ky1= 0.5*(dMU[i]*dXM1[i]-dMU[i+1]*dXM1[i+1]);
				ky2= 0.5*(-dMU[i]*dXM1[i]);
				ky3= 0.5*(dMU[i+1]*dXM1[i+1]);
											
				vals_p.push_back( T(i, i+(L-2), -kx1) );
				vals_p.push_back( T(i, i+1+(L-2), -kx3) );
				vals_p.push_back( T(i, i-1+(L-2), -kx2) );
							
				vals_p.push_back( T(i, i+2*(L-2), ky1) );
				vals_p.push_back( T(i, i+1+2*(L-2), ky3) );
				vals_p.push_back( T(i, i-1+2*(L-2), ky2) );
							
				/*kxy+=(1/6.0)*(1/b1)*(dXM1[i]*dXM2[i]+dX1[i]*dX2[i]+dXM1[i+1]*dXM2[i+1]+dX1[i+1]*dX2[i+1]);
				kxy1+=(1/6.0)*(1/b1)*(dXM1[i]*dXM2[i]);
				kxy2+=(1/6.0)*(1/b1)*(dXM1[i+1]*dXM2[i+1]);
				
				kxx+=(1/6.0)*(1/b1)*(dXM1[i]*dXM1[i]+dX1[i]*dX1[i]+dXM1[i+1]*dXM1[i+1]+dX1[i+1]*dX1[i+1]);
				kxx1+=(1/6.0)*(1/b1)*(dXM1[i]*dXM1[i]);
				kxx2+=(1/6.0)*(1/b1)*(dXM1[i+1]*dXM1[i+1]);
				
				kyy+=(1/6.0)*(1/b1)*(dXM2[i]*dXM2[i]+dX2[i]*dX2[i]+dXM2[i+1]*dXM2[i+1]+dX2[i+1]*dX2[i+1]);
				kyy1+=(1/6.0)*(1/b1)*(dXM2[i]*dXM2[i]);
				kyy2+=(1/6.0)*(1/b1)*(dXM2[i+1]*dXM2[i+1]);*/

				kxy+=(1/3.0)*(1/b1)*(dMU[i]*dXM1[i]*dXM2[i]+dMU[i+1]*dXM1[i+1]*dXM2[i+1]);
				kxy1+=(1/6.0)*(1/b1)*(dMU[i]*dXM1[i]*dXM2[i]);
				kxy2+=(1/6.0)*(1/b1)*(dMU[i+1]*dXM1[i+1]*dXM2[i+1]);
				
				kxx+=(1/3.0)*(1/b1)*(dMU[i]*dXM1[i]*dXM1[i]+dMU[i+1]*dXM1[i+1]*dXM1[i+1]);
				kxx1+=(1/6.0)*(1/b1)*(dMU[i]*dXM1[i]*dXM1[i]);
				kxx2+=(1/6.0)*(1/b1)*(dMU[i+1]*dXM1[i+1]*dXM1[i+1]);
				
				kyy+=(1/3.0)*(1/b1)*(dMU[i]*dXM2[i]*dXM2[i]+dMU[i+1]*dXM2[i+1]*dXM2[i+1]);
				kyy1+=(1/6.0)*(1/b1)*(dMU[i]*dXM2[i]*dXM2[i]);
				kyy2+=(1/6.0)*(1/b1)*(dMU[i+1]*dXM2[i+1]*dXM2[i+1]);
				
				vals_p.push_back( T(i+(L-2), i+(L-2), kyy) );
				vals_p.push_back( T(i+(L-2), i+1+(L-2), kyy2) );
				vals_p.push_back( T(i+(L-2), i-1+(L-2), kyy1 ));
							
				vals_p.push_back( T(i+(L-2), i+2*(L-2), -kxy) );
				vals_p.push_back( T(i+(L-2), i+1+2*(L-2), -kxy2) );
				vals_p.push_back( T(i+(L-2), i-1+2*(L-2), -kxy1) );
							
				vals_p.push_back( T(i+2*(L-2), i+(L-2), -kxy) );
				vals_p.push_back( T(i+2*(L-2), i+1+(L-2), -kxy2) );
				vals_p.push_back( T(i+2*(L-2), i-1+(L-2), -kxy1 ));
									
				vals_p.push_back( T(i+2*(L-2), i+2*(L-2), kxx) );
				vals_p.push_back( T(i+2*(L-2), i+1+2*(L-2), kxx2) );
				vals_p.push_back( T(i+2*(L-2), i-1+2*(L-2), kxx1) );
							
				// Right hand side
				kx+=0.5*(dMU[i]*dXM1[i]+dMU[i+1]*dXM1[i+1]);
				ky+=0.5*(dMU[i]*dXM2[i]+dMU[i+1]*dXM2[i+1]);			
				kp+=(dMU[i]-dMU[i+1]);
				//(*th_b)[(L-2)+i]+= (ky*(1/b1));
				//(*th_b)[2*(L-2)+i]-= (kx*(1/b1));
				
				(*th_b)[(L-2)+i]+= (1.0+a_sum)*(ky*(1/b1));
				(*th_b)[2*(L-2)+i]-= (1.0+a_sum)*(kx*(1/b1));		
				(*th_b)[3*(L-2)+i]-= (kl-k0)*kp;
									
					}

				// u' in perturbe matrix
				vals_p.push_back( T(L-2, 0, -0.5*(dMU[0]*dXM2[0]-dMU[1]*dXM2[1])) );
				vals_p.push_back( T(L-2, 1, -0.5*(dMU[1]*dXM2[1])) );
				vals_p.push_back( T(L-2+L-3, L-3, -0.5*(dMU[L-3]*dXM2[L-3]-dMU[L-2]*dXM2[L-2])) );
				vals_p.push_back( T(L-2+L-3, L-4, 0.5*(dMU[L-3]*dXM2[L-3])) );
	
				vals_p.push_back( T(2*(L-2), 0, 0.5*(dMU[0]*dXM1[0]-dMU[1]*dXM1[1])) );
				vals_p.push_back( T(2*(L-2), 1, 0.5*(dMU[1]*dXM1[1])) );
				vals_p.push_back( T(2*(L-2)+(L-3), L-3, 0.5*(dMU[L-3]*dXM1[L-3]-dMU[L-2]*dXM1[L-2])) );
				vals_p.push_back( T(2*(L-2)+(L-3), L-4, -0.5*(dMU[L-3]*dXM1[L-3])) );

				// mass matrix part
				vals_p.push_back( T(0, L-2, -0.5*(dMU[0]*dXM2[0]-dMU[1]*dXM2[1])) );
				vals_p.push_back( T(0, 1+L-2, 0.5*(dMU[1]*dXM2[1])) );
				vals_p.push_back( T(L-3, L-3+L-2, -0.5*(dMU[L-3]*dXM2[L-3]-dMU[L-2]*dXM2[L-2])) );
				vals_p.push_back( T(L-3, L-4+L-2, -0.5*(dMU[L-3]*dXM2[L-3])) );
				vals_p.push_back( T(0, 2*(L-2), 0.5*(dMU[0]*dXM1[0]-dMU[1]*dXM1[1])) );
				vals_p.push_back( T(0, 1+2*(L-2), -0.5*(dMU[1]*dXM1[1])) );
				vals_p.push_back( T(L-3, L-3+2*(L-2), 0.5*(dMU[L-3]*dXM1[L-3]-dMU[L-2]*dXM1[L-2])) );
				vals_p.push_back( T(L-3, L-4+2*(L-2), 0.5*(dMU[L-3]*dXM1[L-3])) );
					
				/*double kxy_left=(1/6.0)*(1/b1)*(dXM1[0]*dXM2[0]+dXM1[1]*dXM2[1]+dX1[0]*dX2[0]+dX1[1]*dX2[1]);
				double kxy_right=(1/6.0)*(1/b1)*(dXM1[L-3]*dXM2[L-3]+dXM1[L-2]*dXM2[L-2]+dX1[L-3]*dX2[L-3]+dX1[L-4]*dX2[L-4]);
				double kxy1_left=(1/6.0)*(1/b1)*(dXM1[1]*dXM2[1]);
				double kxy1_right=(1/6.0)*(1/b1)*(dXM1[L-3]*dXM2[L-3]);
				double kxx_left=(1/6.0)*(1/b1)*(dXM1[0]*dXM1[0]+dXM1[1]*dXM1[1]+dX1[0]*dX1[0]+dX1[1]*dX1[1]);
				double kxx_right=(1/6.0)*(1/b1)*(dXM1[L-3]*dXM1[L-3]+dXM1[L-2]*dXM1[L-2]+dX1[L-3]*dX1[L-3]+dX1[L-4]*dX1[L-4]);
				double kxx1_left=(1/6.0)*(1/b1)*(dXM1[1]*dXM1[1]);
				double kxx1_right=(1/6.0)*(1/b1)*(dXM1[L-3]*dXM1[L-3]);
				double kyy_left=(1/6.0)*(1/b1)*(dXM2[0]*dXM2[0]+dXM2[1]*dXM2[1]+dX2[0]*dX2[0]+dX2[1]*dX2[1]);
				double kyy_right=(1/6.0)*(1/b1)*(dXM2[L-3]*dXM2[L-3]+dXM2[L-2]*dXM2[L-2]+dX2[L-3]*dX2[L-3]+dX2[L-4]*dX2[L-4]);
				double kyy1_left=(1/6.0)*(1/b1)*(dXM2[1]*dXM2[1]);
				double kyy1_right=(1/6.0)*(1/b1)*(dXM2[L-3]*dXM2[L-3]);*/

				double kxy_left=(1/3.0)*(1/b1)*(dMU[0]*dXM1[0]*dXM2[0]+dMU[1]*dXM1[1]*dXM2[1]);
				double kxy_right=(1/3.0)*(1/b1)*(dMU[L-3]*dXM1[L-3]*dXM2[L-3]+dMU[L-2]*dXM1[L-2]*dXM2[L-2]);
				double kxy1_left=(1/6.0)*(1/b1)*(dMU[1]*dXM1[1]*dXM2[1]);
				double kxy1_right=(1/6.0)*(1/b1)*(dMU[L-3]*dXM1[L-3]*dXM2[L-3]);
				
				double kxx_left=(1/3.0)*(1/b1)*(dMU[0]*dXM1[0]*dXM1[0]+dMU[1]*dXM1[1]*dXM1[1]);
				double kxx_right=(1/3.0)*(1/b1)*(dMU[L-3]*dXM1[L-3]*dXM1[L-3]+dMU[L-2]*dXM1[L-2]*dXM1[L-2]);
				double kxx1_left=(1/6.0)*(1/b1)*(dMU[1]*dXM1[1]*dXM1[1]);
				double kxx1_right=(1/6.0)*(1/b1)*(dMU[L-3]*dXM1[L-3]*dXM1[L-3]);
				
				double kyy_left=(1/3.0)*(1/b1)*(dMU[0]*dXM2[0]*dXM2[0]+dMU[1]*dXM2[1]*dXM2[1]);
				double kyy_right=(1/3.0)*(1/b1)*(dMU[L-3]*dXM2[L-3]*dXM2[L-3]+dMU[L-2]*dXM2[L-2]*dXM2[L-2]);
				double kyy1_left=(1/6.0)*(1/b1)*(dMU[1]*dXM2[1]*dXM2[1]);
				double kyy1_right=(1/6.0)*(1/b1)*(dMU[L-3]*dXM2[L-3]*dXM2[L-3]);

        			vals_p.push_back( T(L-2, L-2, kyy_left ));
				vals_p.push_back( T(L-2, 1+(L-2), kyy1_left ));
				vals_p.push_back( T((L-2)+L-3, L-4+(L-2), kyy1_right ) );
				vals_p.push_back( T((L-2)+L-3, L-3+(L-2), kyy_right) );
	
				vals_p.push_back( T(L-2, 2*(L-2), -kxy_left ));
				vals_p.push_back( T(L-2, 1+2*(L-2), -kxy1_left ));
				vals_p.push_back( T((L-2)+L-3, L-3+2*(L-2), -kxy_right) );
				vals_p.push_back( T((L-2)+L-3, L-4+2*(L-2), -kxy1_right) );
	
				vals_p.push_back( T(2*(L-2), L-2, -kxy_left ));
				vals_p.push_back( T(2*(L-2), 1+(L-2), -kxy1_left ));
				vals_p.push_back( T(2*(L-2)+L-3, L-3+(L-2), -kxy_right ) );
				vals_p.push_back( T(2*(L-2)+L-3, L-4+(L-2), -kxy1_right) );
	
				vals_p.push_back( T(2*(L-2), 2*(L-2), kxx_left ));
				vals_p.push_back( T(2*(L-2), 1+2*(L-2), kxx1_left ));
				vals_p.push_back( T(2*(L-2)+L-3, L-3+2*(L-2), kxx_right) );
				vals_p.push_back( T(2*(L-2)+L-3, L-4+2*(L-2), kxx1_right) );	
				
				// right hand side
				double kx1=0.5*(dMU[0]*dXM1[0]+dMU[1]*dXM1[1]);
				double kx2=0.5*(dMU[L-2]*dXM1[L-2]+dMU[L-3]*dXM1[L-3]);
				double ky1=0.5*(dMU[0]*dXM2[0]+dMU[1]*dXM2[1]);
				double ky2=0.5*(dMU[L-2]*dXM2[L-2]+dMU[L-3]*dXM2[L-3]);
				double kp1= (dMU[0]-dMU[1]);
				double kp2= (dMU[L-3]-dMU[L-2]);
	
			(*th_b)[(L-2)]+= (ky1*(1/b1));
			(*th_b)[(L-2)+(L-3)]+= (ky2*(1/b1));
			(*th_b)[2*(L-2)]-= (kx1*(1/b1));
			(*th_b)[2*(L-2)+(L-3)]-= (kx2*(1/b1));
			(*th_b)[3*(L-2)]-= ((kl-k0)*kp1);
			(*th_b)[3*(L-2)+(L-3)]-= ((kl-k0)*kp2);
   
	(*th_k).setFromTriplets( vals_p.begin(), vals_p.end() );	

}




void ODESystem::PrepPerturbe() {
	std::vector<SpMat> th_k;
	std::vector<Vec> th_b;

	for ( int j = 0; j < numThreads; ++j ) {
		
		th_k.push_back( SpMat(4*(L-2)+3, 4*(L-2)+3) );
		//th_k.push_back( SpMat(3*(L-2), 3*(L-2)) );
		th_b.push_back( Vec() );
		
	}

	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) {
		threads.push_back( std::thread( &ODESystem::PerturbeLoop, this, 
	&th_k[j],&th_b[j])  );
	}

	// join threads
	//for (auto &thread : threads) thread.join();
	for ( int j= 0; j< numThreads; ++j)
	{	
			threads[j].join();	
	}

	P = SpMat(4*(L-2)+3,4*(L-2)+3);
	B = Vec::Zero(4*(L-2)+3);
	
	//P = SpMat(3*(L-2),3*(L-2));
	//B = Vec::Zero(3*(L-2));
	for (int j = 0; j < numThreads; ++j) {
		P += th_k[j];
		B += th_b[j];
	}

}


void ODESystem::Solve(){	

m_size=20;
L=5*m_size;
double m_size2=m_size;
l_c1=1/m_size2;
l_c2=1/m_size2;
sigma_new1=0.3;
sigma_new2=0.3;

// Performing Cholesky Decomposition for geometric perturbations; the decomposition needs to be performed only once
CholeskyLoop();

for (int i=1; i<2; ++i){
	
	Eigen::LLT<Mat> lltOfA(KTL_new1);
	KTL_Cholnew1= lltOfA.matrixL(); 
	 	}

		//Eigen::LLT<Mat> lltOfA(KTL_new2);
	  	//KTL_Cholnew2= lltOfA.matrixL();


//  Number sample of deformations
 int Numb = 500;	
 
for (int i =1; i<=Numb; ++i){
	
Prep_Norm();
Perform_KTL();
PrepK();
PrepK_det();
PrepPerturbe();

// right hand side
B*=(2*pi*E);
K_new= 2.0*E*pi*K1+0.5*E*pi*K2+0.5*E*pi*K3+2.0*E*pi*P+2*pi*E*L1+Schub*pi*K5;

bst.compute(K_new);

std::cerr << "Solving... ";
Z=bst.solve(B);
std::cerr << "done." << std::endl;

CalcE();
std::cerr<< "Energie:" << 0.5*energy1 << std::endl;

	}	
}

// compute elastic energy
void ODESystem::CalcE() {

	energy1=0.0;
	energy2=0.0;
	double b1= (1/(L-1));
	double b2=(L-1);
	double th1=0.0;
	double th2 =0.0;
	double th3=0.0;
	
	th1+=dMU[0]*sqr((1+a_sum)+b2*Z[0]+ 0.5*(Z[2*(L-2)]*dXM1[0]-Z[(L-2)]*dXM2[0]))*(1/b2);	
	th1+=dMU[L-2]*sqr((1+a_sum)-b2*Z[L-3]+ 0.5*(Z[2*(L-2)+(L-3)]*dXM1[(L-2)]-Z[(L-2)+(L-3)]*dXM2[(L-2)]))*(1/b2);

	th2+=(2.0*(kl-k0)*dMU[0]*b2*Z[3*(L-2)]*(1/b2)+sqr(kl-k0)*dMU[0]*(1/b2));
	th2+=(2.0*(kl-k0)*dMU[L-2]*b2*Z[3*(L-2)+(L-3)]*(1/b2)+sqr(kl-k0)*dMU[L-2]*(1/b2));
	
	for (int i=0; i<L-3; ++i) {
	
			th1+=dMU[i+1]*sqr(0.5*dXM1[i+1]*(Z[2*(L-2)+i]+Z[2*(L-2)+i+1])-0.5*dXM2[i+1]*(Z[(L-2)+i]+Z[(L-2)+i+1])+((1+a_sum)+dMU[i+1]*b2*(Z[i+1]-Z[i])))*(1/b2);
			th2+= (2.0*(kl-k0)*dMU[i+1]*(Z[3*(L-2)+i+1]-Z[3*(L-2)+i])*b2*(1/b2)+sqr(kl-k0)*dMU[i+1]*(1/b2));
		
		}
				  
	energy1+= E*pi*th1+(E/(4.0))*pi*(Z.dot(K2*Z)+Z.dot(K3*Z));
}





