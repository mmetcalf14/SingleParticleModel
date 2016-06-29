//
//  SingleParticle_Hamiltonian_MLM.h
//  Single Particle Model
//
//  Created by mekena McGrew on 7/14/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <vector>
#include "/usr/local/include/Eigen/Eigen"
#include "/usr/local/include/Eigen/Sparse"
#include "/usr/local/include/Eigen/Eigenvalues"
using namespace Eigen;
using namespace std;

#ifndef Single_Particle_Model_SingleParticle_Hamiltonian_MLM_h
#define Single_Particle_Model_SingleParticle_Hamiltonian_MLM_h

class Hamiltonian
{
    friend class CorrMat;
private:
    
    double J1, J2, Delta;
    int Ro;
    double y;
    MatrixXd Friendly_Ham;
    
    
protected:
    
    int L;
    MatrixXd Ham_Mat;
    
public:
    
    VectorXd EVal;
    MatrixXd EVec;
    
//    VectorXcd EVal_c;
//    MatrixXcd EVec_c;
    
    
    Hamiltonian();
    Hamiltonian(double _J1, double _J2, int _L, double _Delta, int _Ro, double _y);
    
    void Set_Mat_Dim();
    void Build_Hamiltonian();
    MatrixXd Thouless_Hamiltonian(double t1, double t2, double h);
    void Diagonalize_Hamiltonian(MatrixXd H);
    //void Diagonalize_ThoulessHamiltonian();
    void Write_File_EV();
    inline void FriendHam(){Friendly_Ham = Ham_Mat;}
    
    
    
};

class SlaterDet :public Hamiltonian
{
private:
    int N;
    int count;
    
    vector<int> basis;
    vector<double> G_state;
    MatrixXd SlaterMat;
    vector<double> n_state;
    
    
    
protected:
    
public:
    SlaterDet(double _J1, double _J2, int _L, double _Delta, int _Ro, double _y):Hamiltonian( _J1,  _J2, _L, _Delta, _Ro, _y){};
    void SetNumber(int Npart);
    
    void FockBasis();
    void Get_Gstate_FB();
    void Onsite_Gstate();
    void Write_File(ofstream &OutputFile);
    
    
};

inline int MY_bittest(int m, int n)// m -> basis integer, n -> site
{
    int Eval;//if Eval is size_t I get a totally wrong number compared to int
    //seg fault occurring regardless of whether return value is correct or incorrect
    //std::cout << m << " " << n << std::endl;
    Eval = (m & (1 << n));//I haven't changed anything why not working all of a sudden?
    return Eval;
}


#endif
