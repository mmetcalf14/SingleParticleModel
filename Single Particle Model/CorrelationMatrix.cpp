//
//  CorrelationMatrix.cpp
//  Single Particle Model
//
//  Created by mekena McGrew on 5/10/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "CorrelationMatrix.h"

using namespace Eigen;
void CorrMat::GetNumberVariables(int _Np, int _Nsite)
{
    Np = _Np;
    Nsite = _Nsite;
    I.real(0.0);
    I.imag(1.0);
}



void CorrMat::BuildCorrMat(const Hamiltonian & H)
{
    C_mat = MatrixXcd::Zero(Nsite, Nsite);
    B_mat = MatrixXcd::Zero(Nsite, Nsite);
    
    //cout << "Pre Corr Assign: \n" << C_mat << endl;
    complex<double> sum;
    for( int i = 0; i < Nsite; i++)
    {
        for(int j = 0; j < Nsite; j++)
        {
            sum = 0.;
            for(int k = 0; k < Np; k++)
            {
                sum += H.EVec(i,k)*H.EVec(j,k);
                
            }
            if ( sum.real() <= .0000000001)
            {
                sum = 0.;
            }
            //cout << "sum: " << sum.real() << endl;
            C_mat(i,j) = sum;
        }
    }
    
//    cout << "C_mat initial: \n" << C_mat(0,0).real() << endl;
//    cout << "C_mat initial: \n" << C_mat(12,12).real() << endl;
    
}

void CorrMat::RKPeierlsOnCMat(const Hamiltonian& h, ofstream &fout, double dt, double N_it, double Phi_max, double t_p)
{
    MatrixXcd k1(Nsite, Nsite);
    MatrixXcd k2(Nsite, Nsite);
    MatrixXcd k3(Nsite, Nsite);
    MatrixXcd k4(Nsite, Nsite);
    MatrixXcd Ck_mat(Nsite, Nsite);
    MatrixXcd Cmat_new(Nsite, Nsite);
    
    double t = 0;
    double Phi_t = 0;
    
    int flag = 0;
    
    Matrix2cd T_mat;
    complex<double> t1;//t_{j-1}
    complex<double> t2;//t_j
    complex<double> t3;//t_i
    complex<double> t4;//t_{i-1}
    
    
    for(int it = 0; it <= N_it; it++)
    {
        t = it*dt;
        if(t <= t_p)
        {
        Phi_t = -1*(t/t_p)*Phi_max;
        }
        else
        {
            Phi_t = Phi_max;
        }
        //cout << "t: " << t << " Phi: " << Phi_t << endl;
        
        T_mat(0,0) = h.J1*exp(I*Phi_t);//tbar_{j-1},tbar_i only need this information for if loop
        T_mat(0,1) = h.J2*exp(I*Phi_t);
        T_mat(1,0) = h.J1*exp(-I*Phi_t);//tbar_j, tbar_i-1
        T_mat(1,1) = h.J2*exp(-I*Phi_t);
        //cout << "Problem with Tmat\n";
        
        
        
        for(int i = 0; i < Nsite; i++)//so begins evaluation of k1
            //double loop for each k vector
        {
            

            for(int j = 0; j < Nsite; j++)//only doing upper triangle
            {
                if(((j-1) >= 0) && ((j-1) < (Nsite - 1)))//t_{j-1}
                {
                    if( (j-1)%2 == 0)
                    {
                        t1 = T_mat(0,0)*C_mat(i,j-1);//all information is held in t#
                        //cout << "t1 assigned for even j-1 "<<j-1 <<endl;
                    }
                    else
                    {
                        t1 = T_mat(0,1)*C_mat(i,j-1);
                        //cout << "t1 assigned for odd j-1 "<<j-1 <<endl;
                    }
                }
                else{
                    t1 = 0.;
                    //cout << "t1 is zero "<< j-1 <<endl;
                }
                
                if((j >= 0) && (j < (Nsite - 1)))//t_j
                {
                    if( j%2 == 0)
                    {
                        t2 = T_mat(1,0)*C_mat(i,j+1);
                        //cout << "t2 assigned for even j "<<j <<endl;
                    }
                    else
                    {
                        t2 = T_mat(1,1)*C_mat(i,j+1);
                        //cout << "t2 assigned for odd j "<<j <<endl;
                    }
                }
                else{
                    t2 = 0.;
                    //cout << "t2 is zero "<< j <<endl;
                }
                
                if((i >= 0) && (i < (Nsite - 1)))//t_i
                {
                    if( i%2 == 0)
                    {
                        t3 = T_mat(0,0)*C_mat(i+1,j);
                        //cout << "t3 assigned for even i "<<i <<endl;
                    }
                    else
                    {
                        t3 = T_mat(0,1)*C_mat(i+1,j);
                        //cout << "t3 assigned for odd i "<<i <<endl;
                    }
                }
                else{
                    t3 = 0.;
                }
                if(((i-1) >= 0) && ((i-1) < (Nsite - 1)))//t_{i-1}
                {
                    if( (i-1)%2 == 0)
                    {
                        t4 = T_mat(1,0)*C_mat(i-1,j);
                        //cout << "t4 assigned for even i-1 "<<i-1 <<endl;
                    }
                    else
                    {
                        t4 = T_mat(1,1)*C_mat(i-1,j);
                        //cout << "t4 assigned for odd i-1 "<<i-1 <<endl;
                    }
                }
                else{
                    t4 = 0.;
                   // cout << "t4 is zero "<<i-1 <<endl;
                }
                k1(i,j) = dt*I*(t1+t2-t3-t4);
                Ck_mat(i,j) = C_mat(i,j) + (1./2.)*k1(i,j);
                
                //cout << "k1 is set.\n";
            }
        }//end k1 loop onto k2
        
        for(int i = 0; i < Nsite; i++)//so begins evaluation of k2
            //double loop for each k vector
        {
            
            for(int j = 0; j < Nsite; j++)//only doing upper triangle
            {
                if(((j-1) >= 0) && ((j-1) < (Nsite - 1)))//t_{j-1}
                {
                    if( (j-1)%2 == 0)
                    {
                        t1 = T_mat(0,0)*Ck_mat(i,j-1);//all information is held in t#
                        //cout << "t1 assigned for even j-1 "<<j-1 <<endl;
                    }
                    else
                    {
                        t1 = T_mat(0,1)*Ck_mat(i,j-1);
                        //cout << "t1 assigned for odd j-1 "<<j-1 <<endl;
                    }
                }
                else{
                    t1 = 0.;
                    //cout << "t1 is zero "<< j-1 <<endl;
                }
                
                if((j >= 0) && (j < (Nsite - 1)))//t_j
                {
                    if( j%2 == 0)
                    {
                        t2 = T_mat(1,0)*Ck_mat(i,j+1);
                        //cout << "t2 assigned for even j "<<j <<endl;
                    }
                    else
                    {
                        t2 = T_mat(1,1)*Ck_mat(i,j+1);
                        //cout << "t2 assigned for odd j "<<j <<endl;
                    }
                }
                else{
                    t2 = 0.;
                    //cout << "t2 is zero "<< j <<endl;
                }
                
                if((i >= 0) && (i < (Nsite - 1)))//t_i
                {
                    if( i%2 == 0)
                    {
                        t3 = T_mat(0,0)*Ck_mat(i+1,j);
                        //cout << "t3 assigned for even i "<<i <<endl;
                    }
                    else
                    {
                        t3 = T_mat(0,1)*Ck_mat(i+1,j);
                        //cout << "t3 assigned for odd i "<<i <<endl;
                    }
                }
                else{
                    t3 = 0.;
                }
                if(((i-1) >= 0) && ((i-1) < (Nsite - 1)))//t_{i-1}
                {
                    if( (i-1)%2 == 0)
                    {
                        t4 = T_mat(1,0)*Ck_mat(i-1,j);
                        //cout << "t4 assigned for even i-1 "<<i-1 <<endl;
                    }
                    else
                    {
                        t4 = T_mat(1,1)*Ck_mat(i-1,j);
                        //cout << "t4 assigned for odd i-1 "<<i-1 <<endl;
                    }
                }
                else{
                    t4 = 0.;
                    //cout << "t4 is zero "<<i-1 <<endl;
                }
                k2(i,j) = dt*I*(t1+t2-t3-t4);
                Ck_mat(i,j) = C_mat(i,j) + (1./2.)*k2(i,j);
                //cout << "k2 is set.\n";
            }
        }//end k2 loop
        
        for(int i = 0; i < Nsite; i++)//so begins evaluation of k3
            //double loop for each k vector
        {
            
            for(int j = 0; j < Nsite; j++)//only doing upper triangle
            {
                if(((j-1) >= 0) && ((j-1) < (Nsite - 1)))//t_{j-1}
                {
                    if( (j-1)%2 == 0)
                    {
                        t1 = T_mat(0,0)*Ck_mat(i,j-1);//all information is held in t#
                        //cout << "t1 assigned for even j-1 "<<j-1 <<endl;
                    }
                    else
                    {
                        t1 = T_mat(0,1)*Ck_mat(i,j-1);
                        //cout << "t1 assigned for odd j-1 "<<j-1 <<endl;
                    }
                }
                else{
                    t1 = 0.;
                    //cout << "t1 is zero "<< j-1 <<endl;
                }
                
                if((j >= 0) && (j < (Nsite - 1)))//t_j
                {
                    if( j%2 == 0)
                    {
                        t2 = T_mat(1,0)*Ck_mat(i,j+1);
                        //cout << "t2 assigned for even j "<<j <<endl;
                    }
                    else
                    {
                        t2 = T_mat(1,1)*Ck_mat(i,j+1);
                        //cout << "t2 assigned for odd j "<<j <<endl;
                    }
                }
                else{
                    t2 = 0.;
                    //cout << "t2 is zero "<< j <<endl;
                }
                
                if((i >= 0) && (i < (Nsite - 1)))//t_i
                {
                    if( i%2 == 0)
                    {
                        t3 = T_mat(0,0)*Ck_mat(i+1,j);
                        //cout << "t3 assigned for even i "<<i <<endl;
                    }
                    else
                    {
                        t3 = T_mat(0,1)*Ck_mat(i+1,j);
                        //cout << "t3 assigned for odd i "<<i <<endl;
                    }
                }
                else{
                    t3 = 0.;
                }
                if(((i-1) >= 0) && ((i-1) < (Nsite - 1)))//t_{i-1}
                {
                    if( (i-1)%2 == 0)
                    {
                        t4 = T_mat(1,0)*Ck_mat(i-1,j);
                        //cout << "t4 assigned for even i-1 "<<i-1 <<endl;
                    }
                    else
                    {
                        t4 = T_mat(1,1)*Ck_mat(i-1,j);
                        //cout << "t4 assigned for odd i-1 "<<i-1 <<endl;
                    }
                }
                else{
                    t4 = 0.;
                    //cout << "t4 is zero "<<i-1 <<endl;
                }
                k3(i,j) = dt*I*(t1+t2-t3-t4);
                Ck_mat(i,j) = C_mat(i,j) + k3(i,j);
                
                //cout << "k3 is set.\n";
            }
        }
        
        for(int i = 0; i < Nsite; i++)//so begins evaluation of k4
            //double loop for each k vector
        {
            
            for(int j = 0; j < Nsite; j++)//only doing upper triangle
            {
                if(((j-1) >= 0) && ((j-1) < (Nsite - 1)))//t_{j-1}
                {
                    if( (j-1)%2 == 0)
                    {
                        t1 = T_mat(0,0)*Ck_mat(i,j-1);//all information is held in t#
                       // cout << "t1 assigned for even j-1 "<<j-1 <<endl;
                    }
                    else
                    {
                        t1 = T_mat(0,1)*Ck_mat(i,j-1);
                        //cout << "t1 assigned for odd j-1 "<<j-1 <<endl;
                    }
                }
                else{
                    t1 = 0.;
                    //cout << "t1 is zero "<< j-1 <<endl;
                }
                
                if((j >= 0) && (j < (Nsite - 1)))//t_j
                {
                    if( j%2 == 0)
                    {
                        t2 = T_mat(1,0)*Ck_mat(i,j+1);
                        //cout << "t2 assigned for even j "<<j <<endl;
                    }
                    else
                    {
                        t2 = T_mat(1,1)*Ck_mat(i,j+1);
                        //cout << "t2 assigned for odd j "<<j <<endl;
                    }
                }
                else{
                    t2 = 0.;
                    //cout << "t2 is zero "<< j <<endl;
                }
                
                if((i >= 0) && (i < (Nsite - 1)))//t_i
                {
                    if( i%2 == 0)
                    {
                        t3 = T_mat(0,0)*Ck_mat(i+1,j);
                        //cout << "t3 assigned for even i "<<i <<endl;
                    }
                    else
                    {
                        t3 = T_mat(0,1)*Ck_mat(i+1,j);
                        //cout << "t3 assigned for odd i "<<i <<endl;
                    }
                }
                else{
                    t3 = 0.;
                }
                if(((i-1) >= 0) && ((i-1) < (Nsite - 1)))//t_{i-1}
                {
                    if( (i-1)%2 == 0)
                    {
                        t4 = T_mat(1,0)*Ck_mat(i-1,j);
                        //cout << "t4 assigned for even i-1 "<<i-1 <<endl;
                    }
                    else
                    {
                        t4 = T_mat(1,1)*Ck_mat(i-1,j);
                        //cout << "t4 assigned for odd i-1 "<<i-1 <<endl;
                    }
                }
                else{
                    t4 = 0.;
                    //cout << "t4 is zero "<<i-1 <<endl;
                }
                k4(i,j) = dt*I*(t1+t2-t3-t4);
                //cout << "k4 is set.\n";
            }
        }
        Cmat_new = C_mat+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4;
        C_mat = Cmat_new;
        //cout << "new Corr mat is set.\n";
        
        if(flag == 1000)
        {
            for(int i = 0; i < Nsite; i++)
            {
                fout << i << " " << C_mat(i,i).real() << endl;
            }
            
            //cout <<"Phi: " << Phi_t<< "Cmat evolve: \n"<< C_mat.real() << endl;
            fout << endl;
            flag = 0;
        }
        
        flag++;
        //prepare output formula tomorrow
        T_mat.setZero();

    }
    
    fout.close();
    //for 1D chain there are no nearest neigbors. THis has been taken care of in the equation
}

void CorrMat::RKAdiabatic(Hamiltonian & ham, ofstream &fout, ofstream &Fout, double dt, double N_it, double J1, double J2, double h_0)
{
    
    MatrixXcd k1(Nsite, Nsite);
    MatrixXcd k2(Nsite, Nsite);
    MatrixXcd k3(Nsite, Nsite);
    MatrixXcd k4(Nsite, Nsite);
    MatrixXcd Ck_mat(Nsite, Nsite);
    MatrixXcd Cmat_new(Nsite, Nsite);
    
    double delta_0 = 0.5;
    double delta;
    double J_0 = 1.;
    //double h_0 = 1.;
    double h;
    double omega = 2* 3.14159265;

    //double h_0;
    
    double t = 0;
    //ham.FriendHam();
    
    for(int it = 0; it <= N_it; it++)
    {
        t = it*dt;

            delta = delta_0*cos(omega*t);
            J1 = J_0 + delta;
            J2 = J_0 - delta;
            h = h_0*sin(omega*t);
            ham.Thouless_Hamiltonian(J1, J2, h);//to get H(t)

        
        k1 = dt* I*(C_mat*ham.Friendly_Ham - ham.Friendly_Ham*C_mat);//k1 = f(t,y)->H(t)
        Ck_mat = C_mat + (1./2.)*k1;
        
            delta = delta_0*cos(omega*(t+(0.5*dt)));
            J1 = J_0 + delta;
            J2 = J_0 - delta;
            h = h_0*sin(omega*(t+(0.5*dt)));
            ham.Thouless_Hamiltonian(J1, J2, h); //to get H(t+dt)
        
        k2 = dt* I*(Ck_mat*ham.Friendly_Ham - ham.Friendly_Ham*Ck_mat);//k2 = f(t+1/2dt, y+1/2k1)
        Ck_mat = C_mat + (1./2.)*k2;
            
        k3 = dt* I*(Ck_mat*ham.Friendly_Ham - ham.Friendly_Ham*Ck_mat);//k3 = f(t+1/2dt, y+1/2k2)
        Ck_mat = C_mat + k3;
        
            delta = delta_0*cos(omega*(t+dt));
            J1 = J_0 + delta;
            J2 = J_0 - delta;
            h = h_0*sin(omega*(t+dt));
            ham.Thouless_Hamiltonian(J1, J2, h);//to get H(t)
        
        k4 = dt*I*(Ck_mat*ham.Friendly_Ham - ham.Friendly_Ham*Ck_mat);//k3 = f(t+dt, y+k3)-> H(t+dt)
        
        Cmat_new = C_mat+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4;
        C_mat = Cmat_new;
        

        if(it == 0 || it == N_it/2. || it == N_it)
        {
            ham.Diagonalize_Hamiltonian(ham.Friendly_Ham);//Should I diagonalize H(t) or H(t+dt)?
            //I chose H(t+dt) because Cij(t+dt)
            cout << "Setting Energy \n";
            SetEnergyMat(ham);
            int q = (Nsite/2)+1;
            //complex<double> E = sqrt(conj(B_mat(q,q))*B_mat(q,q));
            
            Fout << J1 << " " << B_mat(q,q).real() << endl;
            
            for(int i = 0; i < Nsite; i++)
                {
                 fout << i << " " << C_mat(i,i).real() << endl;
                }
            fout << endl;
            
            cout << ham.Friendly_Ham << endl;

        }
        
    }
    //end RK4
    
    fout.close();
}

void CorrMat::SetEnergyMat(const Hamiltonian& h)
{
    MatrixXcd EVec_c;
    MatrixXcd EVec_adj;
    //cout << "DOn't xHAve complex Eiegen Mat!\n";
    EVec_c = h.EVec.cast<std::complex<double>>();
    
    EVec_adj = EVec_c.adjoint();
    
    B_mat = EVec_c*C_mat*EVec_adj;
    
    //cout << "B_mat: \n" << B_mat << endl;
    
    C_mat = EVec_adj*B_mat*EVec_c;
    
    //cout << "C matrix: \n" << C_mat << endl;
    
}

//void CorrMat::OnsiteDensity(ofstream &fout)
//{
//    for(int i = 0; i < Nsite; i++)
//    {
//        fout << i << " " << C_mat(i,i) << endl;
//    }
//}




