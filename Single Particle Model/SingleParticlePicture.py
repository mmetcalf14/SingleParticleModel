#!/usr/bin/env python3
# coding=utf-8

import numpy as np
import scipy
import scipy.linalg
import copy
import os
import matplotlib.pyplot as plt

class Hamiltonian(object):
    """docstring for Hamiltonian"""
    def __init__(self, L, dtype=np.complex128):
        super(Hamiltonian, self).__init__()
        self.L = L
        self.dtype = dtype

    def UpdateHopping(self, J1, J2):
        if self.L % 2 == 0:
            num = np.int(self.L / 2 - 1)
            hop = [-J1, -J2] * num
            hop.append(-J1)
        else:
            num = np.int((self.L - 1) / 2)
            hop = [-J1, -J2] * num
        self.hop = np.array(hop, dtype=self.dtype)

    def UpdateOnsiteV(self, V1, V2):
        if self.L % 2 == 0:
            num = np.int(self.L / 2)
            onsite = [V1, V2] * num
        else:
            num = np.int((self.L - 1) / 2)
            onsite = [V1, V2] * num
            onsite.append(V1)
        self.onsite = np.array(onsite, dtype=self.dtype)

    def UpdateOnsite(self, Vlist):
        assert len(Vlist) == self.L, "The inout potential list has wrong dimension!"
        self.onsite = np.array(Vlist, dtype=self.dtype)

    def BuildHamiltonian(self):
        assert len(self.onsite) == len(self.hop) + 1, "Dimension is not fit!"
        self.Ham = np.diag(self.onsite, 0) + np.diag(self.hop, +1) + np.diag(self.hop, -1)
        return self.Ham

    def eigh(self):
        if not hasattr(self, "Ham"):
            self.BuildHamiltonian()
        self.EigVal, self.EigVec = scipy.linalg.eigh(self.Ham)

    def getEigVal(self):
        return self.EigVal

    def getEigVec(self):
        return self.EigVec

    def BuildCorrelationMatrix(self, Ef):
        if not hasattr(self, "EigVal"):
            self.BuildHamiltonian()
            self.eight()
        EMask = self.EigVal < Ef
        NumStates = np.sum(EMask)
        ESpace = np.diag(EMask.astype(self.dtype), k=0)
        Rho = np.dot(self.EigVec.dot(ESpace), np.transpose(np.conj(self.EigVec)))
        return Rho, NumStates
# Thouless Pump

def Equation_of_Motion(A, H):
    Commutator = -1.0j * (A.dot(H) - H.dot(A))
    return Commutator

def RK4(CMt, H1, H2, H4, dt):
    """
    My f() is Equation of Motion.

        d A = -i[A, H(t)]

    # k1 = dt * f( x[i], t )
    # k2 = dt * f( x[i] + 0.5 * k1, t + 0.5 * dt )
    # k3 = dt * f( x[i] + 0.5 * k2, t + 0.5 * dt )
    # k4 = dt * f( x[i] + k3, t + dt )
    # x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
    """
    """ Construct H(t) = H1 """
    """ k1 = dt * f( x[t], t ) """
    C1 = dt * Equation_of_Motion(CMt, H1)

    """ Construct H(t+0.5*dt) = H2 """
    """ Construct x[t] + 0.5 * k1 """
    Cw = CMt + 0.50 * C1
    """ k2 = dt * f( x[i] + 0.5 * k1, t[i] + 0.5 * dt ) """
    C2 = dt * Equation_of_Motion(Cw, H2)

    """
    H3 = H2
    Construct x[t] + 0.5 * k2
    """
    Cw = CMt + 0.50 * C2
    """ k3 = dt * f( x[i] + 0.5 * k2, t[i] + 0.5 * dt ) """
    C3 = dt * Equation_of_Motion(Cw, H2)

    """ Construct H(t+dt) = H4 """
    """ Construct x[i] + k3 """
    Cw = CMt + C3
    """ k4 = dt * f( x[i] + k3, t + dt ) """
    C4 = dt * Equation_of_Motion(Cw, H4)

    """ x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0 """
    return CMt + (C1 + 2.0 * (C2 + C3) + C4) / 6.0


# Initial Condition

L = 111
Ham = Hamiltonian(L, dtype=np.complex128)

J = 1.0e0
dJ = 0.40e0
J1 = J - dJ
J2 = J + dJ
Ham.UpdateHopping(J1, J2)

V1 = 0.0e0
V2 = 0.0e0
Ham.UpdateOnsiteV(V1, V2)

Ham.eigh()

w = Ham.getEigVal()
v = Ham.getEigVec()

Ef = 0.5
CM, Ns = Ham.BuildCorrelationMatrix(Ef)

Tr = 40.0
omega = 2.0* np.pi / Tr
dtime = 1.0e-2
Tsteps = 100 * (np.int(Tr) + 1 )
Tf = Tsteps * dtime
V = 0.40e0

CMt = copy.copy(CM)
Density = [np.diag(CM).real,]
ts = 0
while ts < Tsteps:
    time = ts * dtime
    t1 = time
    J1 = J - dJ * np.cos(omega * t1)
    J2 = J + dJ * np.cos(omega * t1)
    Ham.UpdateHopping(J1, J2)
    V1 = -V * np.sin(omega * t1)
    V2 = V * np.sin(omega * t1)
    Ham.UpdateOnsiteV(V1, V2)
    H1 = copy.copy(Ham.BuildHamiltonian())
    if ts == 0 or ts == 4000 or ts == 2000:
        print(time, J1, J2, V1, V2)

    t2 = time + 0.50 * dtime
    J1 = J - dJ * np.cos(omega * t2)
    J2 = J + dJ * np.cos(omega * t2)
    Ham.UpdateHopping(J1, J2)
    V1 = -V * np.sin(omega * t2)
    V2 = V * np.sin(omega * t2)
    Ham.UpdateOnsiteV(V1, V2)
    H2 = copy.copy(Ham.BuildHamiltonian())

    t4 = time + dtime
    J1 = J - dJ * np.cos(omega * t4)
    J2 = J + dJ * np.cos(omega * t4)
    Ham.UpdateHopping(J1, J2)
    V1 = -V * np.sin(omega * time)
    V2 = V * np.sin(omega * time)
    Ham.UpdateOnsiteV(V1, V2)
    H4 = copy.copy(Ham.BuildHamiltonian())
    CMt = RK4(CMt, H1, H2, H4, dtime)
    Density.append(np.diag(CMt).real)

    ts += 1

fig = plt.figure(figsize=(16,6))
ax1 = fig.add_subplot(121)
ax1.imshow(Density, extent=(0, L-1, 0, Tf), origin='lower', aspect='auto', cmap=plt.get_cmap('YlGnBu'))
ax1.set_xlabel(r"Site Index $i$", fontsize=30)
ax1.set_ylabel(r"$t/t_0$", fontsize=30)
# ticks_size(ax1, 25)
# frame_width(ax1, 2)

ax2 = fig.add_subplot(122)
ax2.plot(Density[0], "o-", color='b', mec='b', label=r"$t=0$", ms=8)
ax2.plot(Density[np.int(0.5*Tr/dtime)], "^-", color='r', mec='r', label=r"$t=20t_0$", ms=8)
ax2.legend(loc="upper center", fontsize=30, frameon=0, numpoints=1)
ax2.set_xlim(-1, L)
ax2.set_xlabel(r"Site Index $i$", fontsize=30)
ax2.set_ylim(0.35, 1.05)
ax2.set_ylabel(r"$N_i$", fontsize=30)
# ticks_size(ax2, 25)
# frame_width(ax2, 2)

plt.savefig("Pump.pdf", format='pdf')
