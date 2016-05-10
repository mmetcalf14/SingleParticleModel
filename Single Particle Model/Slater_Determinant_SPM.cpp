//
//  Slater_Determinant_SPM.cpp
//  Single Particle Model
//
//  Created by mekena McGrew on 8/12/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "SingleParticle_Hamiltonian_MLM.h"


void SlaterDet::SetNumber(int Npart)
{
    N = Npart;
    FockBasis();
}

void SlaterDet::FockBasis()
{
    cout << "Creating basis \n";
    int minrange = 0, maxrange = 0;
    for (int i = 1; i <= N; i++)
    {
        
        minrange += pow(2,(i-1));
        maxrange += pow(2,(L-i));
        
    }
    
    for (int i = minrange; i <= maxrange; i++) //create spin up basis and vectors
    {
        int nbit = 0;
        for(int j = 0; j < L; j++)
        {
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == N)
        {
            count++;
            basis.push_back(i);
        }
    }
    cout << "Basis complete \n";
//    for(int i = 0; i < basis.size(); i++)
//    {
//        cout << basis[i] << endl;
//    }
    //end function
}

void SlaterDet::Get_Gstate_FB()
{
    SlaterMat.resize(N,N);
    
    cout << "Determining Ground state \n";
    
    for(int ct = 0; ct < count; ct++)
    {
        int bas = basis[ct];
        for(int n = 0; n < N; n++)
        {
            int i = 0;
            for(int l = 0; l < L; l++)
            {
                if(MY_bittest(bas, l))
                {
                    SlaterMat(n,i) = EVec(l,n);
                    i++;
                }
            }
            
            
        }
        //cout << "SlaterMat: \n" << SlaterMat << endl;
        G_state.push_back(SlaterMat.determinant());//only works up to 4x4 mat
        //should I set this to zero
    }
//       cout << "The ground state is: \n";
//    for(int i = 0; i < G_state.size(); i++)
//    {
//        cout << G_state[i] << endl;
//    }
    
        
}

void SlaterDet::Onsite_Gstate()
{
    n_state.resize(L,0.0);
    cout << "Onsite ground state\n";
    for(int ct = 0; ct < count; ct++)
    {
        double cf = G_state[ct]*G_state[ct];
        for(int n =0; n < L; n++)
        {
            int bas = basis[ct];
            if(MY_bittest(bas, n))
            {
                n_state.at(n) +=cf;
            }
        }
    }
    
    double Part_tot = 0.0;
    cout << "The onsite ground state is: \n";
    for(int i = 0; i < n_state.size(); i++)
    {
        Part_tot += n_state[i];
        cout << n_state[i] << endl;
    }
    cout << "Total Particle number? " << Part_tot << endl;

}
void SlaterDet::Write_File(ofstream &OutputFile)
{
    for(int i = 0; i < n_state.size(); i++)
    {
        OutputFile << i << " " << n_state[i] << endl;
    }
}


