//
//  SingleParticle_Hamiltonian_MLM.cpp
//  Single Particle Model
//
//  Created by mekena McGrew on 7/14/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>

#include "SingleParticle_Hamiltonian_MLM.h"
using namespace std;



Hamiltonian::Hamiltonian(double _J1, double _J2, int _L, double _Delta, int _Ro, double _y)
{
    J1 = _J1;
    J2 = _J2;
    L = _L;
    Delta = _Delta;
    Ro = _Ro;
    y = _y;
    
    Set_Mat_Dim();
}

void Hamiltonian::Set_Mat_Dim()
{
   
    Ham_Mat.resize(L,L);
    for(int i = 0; i < L; i++)
    {
        for(int j =0 ; j < L; j++)
            Ham_Mat(i,j) = 0;
    }
    
    Build_Hamiltonian();
    
}

void Hamiltonian::Build_Hamiltonian()
{
    
    for(int i = 0; i < L; i++)
    {
        if( (i%2) == 0)
        {
            
            if( (i+1) < L)
            {Ham_Mat(i+1, i) = -J1;}
            //Ham_Mat(i,i) = Delta;
            Ham_Mat(i,i) = y*((i+1)-Ro)*((i+1)-Ro);
        }
        else
        {
            if( (i+1) < L)
            {
               
                Ham_Mat(i+1, i) = -J2;
            }
            //Ham_Mat(i,i) = -Delta;
            Ham_Mat(i,i) = y*((i+1)-Ro)*((i+1)-Ro);
        }
    }
    
//    for(int i = 0; i < L; i++)
//    {
//        cout << Ham_Mat(i,i) << endl;
//    }
    MatrixXd tmp = 0.50 * (Ham_Mat + Ham_Mat.transpose());
    //it is necessary to construct the full Hermitian matrix in order to use compute()
    //which computes tri-diagonal
    //if computing eigenvalues using tridiag method create full matrix!
    Ham_Mat = tmp;
    cout << "This is the Hamiltonian: \n" << Ham_Mat << endl;
    Diagonalize_Hamiltonian();
}

void Hamiltonian::Diagonalize_Hamiltonian()
{
 SelfAdjointEigenSolver<MatrixXd> Diag(Ham_Mat);
    EVal = Diag.eigenvalues();
    EVec = Diag.eigenvectors();
    
//    
//    cout << "Here are the eigenvalues: \n" << EVal << endl;
//    ofstream output;
//    output.open("AlternatingChain_J1-1_J2-2_Delta0_y0.4_L10_EnergyEV_MLM_032316.dat");
//    cout << "Eval size: " << EVal.size() << endl;
//    VectorXd numbers(L);
//    for(int i = 0; i < L; i++)
//    {
//        numbers(i) = i;
//        //cout << "numbers: " << numbers(i) << endl;
//    }
//    output << EVal << endl;
//    output.close();
//
//    cout << "Here are the eigenvectors: \n" << EVec << endl;
   // Write_File_EV();
    
    
}

void Hamiltonian::Write_File_EV()
{
    ofstream fout;
    fout.open("AlternatingChain_J1isJ2_Delta0_y1_L7_N5_Prob_MLM_090115.dat");
    
    VectorXd n_vec(L);
    for(int i = 0; i < L; i++)
    {
        n_vec(i) = i;
    }
    
    if((L%2) != 0)
    {   
        for(int i = 0; i < L ; i++)
        {
            
            fout << n_vec(i)  << " "<< EVec.col(((L+1)/2)-1).row(i)*EVec.col(((L+1)/2)-1).row(i)  <<" " << EVec.col(L-1).row(i)*EVec.col(L-1).row(i)<< " " << EVec.col(L-2).row(i)*EVec.col(L-2).row(i)<<  endl;//" " << EVec.col(22).row(i)*EVec.col(22).row(i)<<
        }
        cout << "Highest Energy Eigenvecotr: \n"<< EVec.col(L-3) << endl;
        cout << "Second Highest Energy Eigenvecotr: \n"<< EVec.col(L-4) << endl;
    }
    else
    {
        for(int i = 0; i < L ; i++)
        {
         //fout << n_vec(i) << " " << EVec.col(L/2).row(i).adjoint()*EVec.col(L/2).row(i) <<endl;//<< EVec.col(19).row(i) << " "
            cout << "Highest Energy Eigenvecotr: \n"<< EVec.col(L-1) << endl;
        }
    }
    
    fout.close();
}







