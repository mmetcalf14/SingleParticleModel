//
//  main.cpp
//  Single Particle Model
//
//  Created by mekena McGrew on 7/14/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>

#include "/usr/local/include/Eigen/Eigen"
#include "/usr/local/include/Eigen/Sparse"
#include "SingleParticle_Hamiltonian_MLM.h"
#include "CorrelationMatrix.h"

int main(int argc, const char * argv[])
{
    using namespace std;
    
    int Nsite;//keep odd number sites to pair correctly
    int Npart;
    double J1;
    double J2;
    double Delta = 0.0;
    int Ro;//(Nsite +1)/2;
    double y = 0.;//y = 1/2*m*\omega^2
    double h;
    
    double dt;
    int Tfinal;

    double Pi = 3.14159265;
    double Phi_m = Pi/2.;
    double t_phi;
    char output[100];
    char AD_OP[100];
    
    if(Nsite % 2 == 0)
    {
        Ro = (Nsite)/2;
    }
    else
    {
        Ro = (Nsite +1)/2;
    }
    
    ofstream fout;
    fout.open("AlternatingChain_J1-2_J2-1_Delta0_y0.01_L21_N10_OnsiteDensity_MLM_032316.dat");
    
    
    ifstream read_file("AlternatingChain_DensityProfile_RK4_DataInput.cfg");
    assert(read_file.is_open());
    
    read_file >> Nsite;
    read_file >> Npart;
    read_file >> J1;
    read_file >> J2;
    read_file >> Tfinal;
    read_file >> dt;
    read_file >> t_phi;
    read_file >> h;
    read_file >> output;
    read_file >> AD_OP;
    
    cout << Nsite << endl;
    cout << h << endl;
    cout << output << endl;
    read_file.close();
    
    ofstream Fout(output);
    assert(Fout.is_open());
    Fout.setf(ios::scientific);
    Fout.precision(11);
    
    ofstream FFout(AD_OP);
    assert(FFout.is_open());
    FFout.setf(ios::scientific);
    FFout.precision(11);


    
    int Time_it = Tfinal/dt;
    
    cout << "1\n";
    Hamiltonian ham(J1, J2, Nsite, Delta, Ro, y);
    //SlaterDet SD(J1, J2, Nsite, Delta, Ro, y);
    CorrMat CM(ham);

//    //ham.Write_File(fout);
//    //fout << Hamiltonian::EVal << endl;
    
//Running Dynamics using microcanonical formalism
    CM.GetNumberVariables(Npart, Nsite);
    //CM.SetEnergyMat();//building <b^{\dagger}b> = \theta(N_p-k)this is unnecessary
    cout << "2\n";
    CM.BuildCorrMat(ham);
    cout << "3\n";
    //CM.RungeKuttaOnCMat(ham, Fout, dt, Time_it, Phi_m, t_phi);
    CM.RKAdiabatic(ham, Fout, FFout, dt, Time_it, J1, J2, h);
    

//    SD.SetNumber(Npart);
//    SD.Get_Gstate_FB();
//    SD.Onsite_Gstate();
//    SD.Write_File(fout);
    
    cout << "Done \n";
    fout.close();
    
    
    return 0;
}
