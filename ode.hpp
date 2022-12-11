#ifndef __PDE_HPP__
#define __PDE_HPP__

#include <string>

#include <thread>
#include <mutex>

#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>
//#include "spline.h"

using namespace std;
#include <random>

class ODESystem {

public:
	typedef Eigen::VectorXd Vec;
        //typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
	typedef Eigen::MatrixXd Mat;
        typedef Eigen::Triplet<double> T;
        typedef Eigen::Spline<double,1> Spline1d;
        typedef Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> >  BCGST;
   
  
  	void KLoop(SpMat*,SpMat*,SpMat*,SpMat*,SpMat*);
 	 void KLoop_det(SpMat*);
  	void Gen_Norm();
  	void CholeskyLoop();
  	void PerturbeLoop(SpMat*, Vec*);
  	void PrepK();
  	void PrepK_det();
  
  	void Solve(); 
  	void Prep_Norm();
  	void Perform_KTL();
  	void PrepPerturbe();
  	void CalcE();
  
private:

SpMat K1,K2,K3,K5,K6,P,K_new,L1;
Mat KTL_new1, KTL_new2,KTL_new3,KTL_Cholnew1,KTL_Cholnew2,KTL_Cholnew3;
Vec B, Z,C_norm1,C_norm2,C_norm3,X2,X1,dX1,dX2,dXM1,dXM2,dMU;

int L, m_size; 
double l_c1,l_c2,l_c3,sigma_new1,sigma_new2,sigma_new3,kernel_b,energy1,energy2,mu1,mu2;

int numThreads=1; // number of threads
//int numThreads=4; // number of threads
 
 BCGST bst;


};


#endif
