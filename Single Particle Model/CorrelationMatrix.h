//
//  CorrelationMatrix.h
//  Single Particle Model
//
//  Created by mekena McGrew on 5/10/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#include "SingleParticle_Hamiltonian_MLM.h"
#include <cmath>
#ifndef CorrelationMatrix_h
#define CorrelationMatrix_h
using namespace Eigen;

class CorrMat
{
private:
    
    MatrixXd B_mat;
    MatrixXcd C_mat;
    
    int Nsite;
    int Np;
    
    complex<double> I;
    complex<double> tbar;
    
protected:
    
public:
    
    CorrMat(const Hamiltonian&){};
    void GetNumberVariables(int _Np, int _Nsite);
    void SetEnergyMat();
    void BuildCorrMat(const Hamiltonian&);
    void RungeKuttaOnCMat(const Hamiltonian&, ofstream &fout, double dt, double N_it, double Phi_max, double t_p);
    //void OnsiteDensity(ofstream &fout);
    
};

#endif /* CorrelationMatrix_h */
