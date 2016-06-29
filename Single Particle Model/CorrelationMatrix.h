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
    
    MatrixXcd B_mat;
    MatrixXcd C_mat;
    
    int Nsite;
    int Np;
    
    complex<double> I;
    complex<double> tbar;
    
protected:
    
public:
    
    CorrMat(const Hamiltonian&){};
    void GetNumberVariables(int _Np, int _Nsite);
    void SetEnergyMat(const Hamiltonian& h);
    void BuildCorrMat(const Hamiltonian&);
    void RKPeierlsOnCMat(const Hamiltonian&, ofstream &fout, double dt, double N_it, double Phi_max, double t_p);
    void RKAdiabatic(Hamiltonian& ham, ofstream &fout, ofstream &Fout, double dt, double N_it, double T, double J1, double J2, double h_0);
    void ClearHam();
    
    //void RKThouless(const Hamiltonian&, ofstream &fout, double dt, double N_it, double Phi_max, double t_p);
    //void OnsiteDensity(ofstream &fout);
    
};

#endif /* CorrelationMatrix_h */
