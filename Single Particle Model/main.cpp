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

int main(int argc, const char * argv[])
{
    using namespace std;
    
    int Nsite = 10;//keep odd number sites to pair correctly
    int Npart = 6;
    double J1 = 1.0;
    double J2 = 1.0;
    double Delta = 0.0;
    int Ro;//(Nsite +1)/2;
    double y = 0.;//y = 1/2*m*\omega^2
    
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


        
    
    Hamiltonian ham(J1, J2, Nsite, Delta, Ro, y);
    SlaterDet SD(J1, J2, Nsite, Delta, Ro, y);
    
    ham.Set_Mat_Dim();
    ham.Build_Hamiltonian();
    ham.Diagonalize_Hamiltonian();
    //ham.Write_File(fout);
    //fout << Hamiltonian::EVal << endl;
    
    SD.SetNumber(Npart);
    SD.Get_Gstate_FB();
    SD.Onsite_Gstate();
    SD.Write_File(fout);
    
    cout << "Done \n";
    fout.close();
    
    return 0;
}
