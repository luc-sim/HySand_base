# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 11:41:33 2025

@author: Luc Simonin

Copy of HySand, taking of consolidation

"""
# import autograd.numpy as np
import numpy as np
from HyperDrive import Utils as hu

check_eps = np.array([0.01,0.04])
check_sig = np.array([60.1,13.1])
check_alp = np.array([[0.05,0.13], [0.09,0.18], [0.02,0.05],
                      [0.01,0.04], [0.005,0.36], [0.02,0.1], [0.1,2.1]])
check_chi = np.array([[9.0,0.1], [10,0.2], [11.0,0.3],
                      [14.0,0.1],[15.0,0.2], [14.0,0.3],[5.0,3.0]])

file = "Dilation and anisotropy"
name = "Dilation and anisotropy file"
mode = 1
ndim = 2
n_y = 3
n_int = 7
n_inp = 1
const =      [  3 ,1.8 ,50000,40000,0.7, 34   ,  0.8 , 2.057, 1.98  , 1.677 ,0.01,0.002, 60 ,0.25,10 , 5 , 0.0]
name_const = ['NN','vi', "K" , "G" ,'m',"phic","bmax",'Beta','Gamma','Delta','lb','ld' ,'A0','r' ,'h','b','a0']


def deriv():
    global NN,vi,K,G,m,phic,bmax,Be,Ga,De,lb,ld,A0,r,h0,b,a0
    NN = int(const[0])
    [vi,K,G,m,phic,bmax,Be,Ga,De,lb,ld,A0,r,h0,b,a0] =  [float(i) for i in const[1:]]
    global rN, n_int, n_y
    rN = 1/NN
    n_int = 2*NN+1
    n_y = NN
    global pref, kref, gref
    pref = 100
    kref = K/pref
    gref = G/pref
    global mu
    mu = np.tan(phic*np.pi/180)
    global H0
    H0 = np.zeros([NN])    
    for i in range(NN):
        H0[i]=h0*(1-(i+1)*rN)**b
    global eps0
    eps0 = 0
    global count
    count = 0
deriv()
def setalp(alp):
    pr = 0
    if pr==1:
        print("\nIN SETALP:")
    pi = 10**-6
    pc0 = 0.1
    # Calculating the initial elastic strain
    eps_p0 = -dgds(np.array([pi,0]),np.zeros([n_int,ndim]))[0]
    if pr==1:
        print("eps_p0:",eps_p0)
    # Calculating the normal consolidation line slope at the test point of initiation (pi=10**-6,0), vi
    pp_i = ( (pi/pref)**(1-m) - 1 ) / (1-m)
    lnvlratio = (np.log(Be)-(np.log(vi))-lb*pp_i) / (np.log(Be/De)-(lb-ld)*pp_i)
    lam_i = lb - (lb-ld)*lnvlratio
    alp_pci = (lam_i-1/kref) / (1-m) * (pc0/pref)**(1-m)
    if pr==1:
        print("lam test:",lam_i)
        print("alp_pci:",alp_pci)
    for i in range(NN):
        alp[NN+i,0] = alp_pci
    # Volumetric strain at initialisation at (p0,0) for v0
    global eps0
    eps0 = eps_p0 + alp_pci
    # eps0 = 0
    if pr==1:
        print("eps0:",eps0)
        print("\n")
    return(alp)

def update(t,eps,sig,alp,chi,dt,deps,dsig,dalp,dchi):
    global count
    if count<1:
        alp[NN,1] = a0      # Initialising anisotropy
        count = 1000
    if alp[-1,0] > 0:
        if sig[1]<=0:
            alp[-1,1] +=1   # Used to count number of cycles for 2-way cycles
    alp[-1,0] = sig[1]
    return alp

def pm(sig):
    pm =  ((sig[0])**2 + (kref/(3*gref)) * (1-m) * sig[1]**2 )**(1/2)
    return (pm)

def g(sig,alp):
    temp= - pref / ( kref*(1-m)*(2-m) ) * (pm(sig)/pref)**(2-m) \
          - rN * sig[0] * sum(alp[i,0] for i in range(NN)) \
          - rN * sig[0] * sum(alp[NN+i,0] for i in range(NN)) \
          - rN * sig[1] * sum(alp[i,1] for i in range(NN))
    return (temp)
def dgds(sig,alp):
    temp=-np.array([ ((sig[0]) / (pref*kref*(1-m)) * ( ((sig[0])**2 + kref*(1-m)*sig[1]**2/(3*gref)) / pref**2)**(-m/2))  + rN * sum(alp[i,0] for i in range(2*NN)),
                      (sig[1]/(3*gref*pref) * ( ((sig[0])**2 + kref*(1-m)*sig[1]**2/(3*gref)) / pref**2)**(-m/2))            + rN * sum(alp[i,1] for i in range(NN)) ])
    return (temp)    
def dgda(sig,alp):
    temp = np.zeros([n_int,ndim])
    for i in range(NN):
        temp[i,0] = -rN * sig[0]
        temp[NN+i,0] = -rN * sig[0]
        temp[i,1] = - rN * sig[1]
    return temp
def d2gdsds(sig,alp):
    brack=( ((sig[0])**2 + kref*(1-m)*sig[1]**2 / (3*gref)) / pref**2)
    temp=-np.array([[ (1 / (kref*(1-m)*pref) * brack**(-m/2) - m*(sig[0])**2 / ((1-m)*kref*pref**3) * brack**(-m/2-1))  
                            ,  -m*(sig[0])*sig[1] / (3*gref*pref**3) * brack**(-m/2-1)] ,
                      [-m*(sig[0])*sig[1] / (3*gref*pref**3) * brack**(-m/2-1)
                            ,  1 / (3*gref*pref) * brack**(-m/2) - m*kref*(1-m)*sig[1]**2 / (((3*gref)**2)*pref**3)  * brack**(-m/2-1) ]]) 
    return (temp)
def d2gdsda(sig,alp):
    temp = np.zeros([ndim,n_int,ndim])
    for j in range(NN):
        temp[0,j,0] = - rN
        temp[0,NN+j,0] = - rN
        temp[1,j,1] = - rN
    return temp
def d2gdads(sig,alp):
    temp = np.zeros([n_int,ndim,ndim])
    for j in range(NN):
        temp[j,0,0] = - rN
        temp[NN+j,0,0] = - rN
        temp[j,1,1] = - rN
    return temp
def d2gdada(sig,alp):
    temp = np.zeros([n_int,ndim,n_int,ndim])
    return temp

y_exclude = True

# denominator of Matsuoka-Nakai gauge function
def f_MN_pq(sig):
    return ( 4 * mu**2 * (sig[0]+2/3*sig[1])*(sig[0]-sig[1]/3) )
def d_f_MN_pq_ds(sig):
    temp = np.zeros([ndim])
    temp[0] = 4 * mu**2 * ( (sig[0]+2/3*sig[1]) + (sig[0]-sig[1]/3) )
    temp[1] = 4 * mu**2 * ( -1/3 * (sig[0]+2/3*sig[1]) + 2/3 * (sig[0]-sig[1]/3) )
    return temp

# Density ratio comparing the distance between equivalent density from loosest state at reference pressure 
# to distance between densest and loosest states at reference pressure
def f_lnvlratio(eps,sig):
    pp = ( (sig[0]/pref)**(1-m) - 1 ) / (1-m)
    temp = (np.log(Be)-(np.log(vi)-(eps[0]-eps0))-lb*pp) / (np.log(Be/De)-(lb-ld)*pp)
    return temp
def d_f_lnvlratio_de(eps,sig):
    temp = np.zeros([ndim])
    pp = ( (sig[0]/pref)**(1-m) - 1 ) / (1-m)
    temp[0] = 1 / (np.log(Be/De)-(lb-ld)*pp)
    return temp
def d_f_lnvlratio_ds(eps,sig):
    temp = np.zeros([ndim])
    pp = ( (sig[0]/pref)**(1-m) - 1 ) / (1-m)
    dpp_dp = (1-m) * sig[0]**(-m) * (pref)**(m-1)  / (1-m)
    u = (np.log(Be)-(np.log(vi)-(eps[0]-eps0))-lb*pp)
    du = -lb*dpp_dp
    v = (np.log(Be/De)-(lb-ld)*pp)
    dv = -(lb-ld)*dpp_dp
    temp[0] = (du*v-u*dv)/v**2
    return temp

# Consolidation pressure(s)
def f_pc(eps,sig,alp):
    lnvlratio = f_lnvlratio(eps,sig)
    temp = np.zeros([NN])
    lam = lb - (lb-ld)*lnvlratio
    for i in range(NN):
        temp[i] = pref * ( (1-m)/(lam-1/kref) * alp[NN+i,0] )**(1/(1-m))
    return temp
def d_f_pc_de(eps,sig,alp):
    lnvlratio = f_lnvlratio(eps,sig)
    dlnvlratio_de = d_f_lnvlratio_de(eps, sig)
    temp = np.zeros([NN,ndim])
    lam = lb - (lb-ld)*lnvlratio
    dlam = - (lb-ld)*dlnvlratio_de
    for i in range(NN):
        temp[i] = pref * ( (1-m) * alp[NN+i,0] )**(1/(1-m)) /(m-1) * dlam * (lam-1/kref)**(1/(m-1)-1)
    return temp
def d_f_pc_ds(eps,sig,alp):
    lnvlratio = f_lnvlratio(eps,sig)
    dlnvlratio_ds = d_f_lnvlratio_ds(eps, sig)
    temp = np.zeros([NN,ndim])
    lam = lb - (lb-ld)*lnvlratio
    dlam = - (lb-ld)*dlnvlratio_ds
    for i in range(NN):
        temp[i] = pref * ( (1-m)  * alp[NN+i,0] )**(1/(1-m)) / (m-1) * dlam*(lam-1/kref)**(1/(m-1)-1)
    return temp
def d_f_pc_da(eps,sig,alp):
    lnvlratio = f_lnvlratio(eps,sig)
    temp = np.zeros([NN,n_int,ndim])
    lam = lb - (lb-ld)*lnvlratio
    for i in range(NN):
        temp[i,NN+i,0] = pref * ( (1-m)/(lam-1/kref) )**(1/(1-m)) / (1-m) * alp[NN+i,0]**(1/(1-m)-1)
    return temp

# Modified deviatoric generalised stress
def f_CHI(eps,sig,alp,chi):
    lnvlratio = f_lnvlratio(eps,sig)
    temp = np.zeros([NN])
    H = H0 * lnvlratio
    beta_c = bmax * ( np.log(Be/Ga) - lnvlratio*np.log(Be/De) ) / np.log(De/Ga)
    A = A0 * lnvlratio
    for i in range(NN):
        temp[i] = (chi[i,1]*NN \
                  - 3*H[i]*sig[0]*alp[i,1] \
                  - alp[NN,1] * ((i+1)*rN)*beta_c*chi[i,0]*NN \
                + A * ( 1 - alp[NN,1] * np.sign(alp[NN,1]) ) * chi[NN,1] )
    return temp
def d_f_CHI_de(eps,sig,alp,chi):
    lnvlratio = f_lnvlratio(eps,sig)
    dlnvlratio_de = d_f_lnvlratio_de(eps, sig)
    temp = np.zeros([NN,ndim])
    H = H0 * lnvlratio
    dH = H0 * dlnvlratio_de[0]
    beta_c = bmax * ( np.log(Be/Ga) - lnvlratio*np.log(Be/De) ) / np.log(De/Ga)
    dbeta_c = bmax * ( - dlnvlratio_de[0]*np.log(Be/De) ) / np.log(De/Ga)
    A = A0 * lnvlratio
    dA = A0 * dlnvlratio_de[0]
    for i in range(NN):
        temp[i,0] = - 3*dH[i]*sig[0]*alp[i,1] \
                  - alp[NN,1] * ((i+1)*rN)*dbeta_c*chi[i,0]*NN \
                + dA * ( 1 - alp[NN,1] * np.sign(alp[NN,1]) ) * chi[NN,1]
    return temp
def d_f_CHI_ds(eps,sig,alp,chi):
    lnvlratio = f_lnvlratio(eps,sig)
    dlnvlratio_ds = d_f_lnvlratio_ds(eps, sig)
    temp = np.zeros([NN,ndim])
    H = H0 * lnvlratio
    dH = H0 * dlnvlratio_ds[0]
    beta_c = bmax * ( np.log(Be/Ga) - lnvlratio*np.log(Be/De) ) / np.log(De/Ga)
    dbeta_c = bmax * ( - dlnvlratio_ds[0]*np.log(Be/De) ) / np.log(De/Ga)
    A = A0 * lnvlratio
    dA = A0 * dlnvlratio_ds[0]
    for i in range(NN):
        temp[i,0] = - 3*dH[i]*sig[0]*alp[i,1] \
                  - 3*H[i]*alp[i,1] \
                  - alp[NN,1] * ((i+1)*rN)*dbeta_c*chi[i,0]*NN \
                + dA * ( 1 - alp[NN,1] * np.sign(alp[NN,1]) ) * chi[NN,1]
    return temp
def d_f_CHI_da(eps,sig,alp,chi):
    lnvlratio = f_lnvlratio(eps,sig)
    temp = np.zeros([NN,n_int,ndim])
    H = H0 * lnvlratio
    beta_c = bmax * ( np.log(Be/Ga) - lnvlratio*np.log(Be/De) ) / np.log(De/Ga)
    A = A0 * lnvlratio
    for i in range(NN):
        temp[i,i,1] = - 3*H[i]*sig[0]
        temp[i,NN,1] = - ((i+1)*rN)*beta_c*chi[i,0]*NN \
                - A * np.sign(alp[NN,1]) * chi[NN,1] 
    return temp
def d_f_CHI_dc(eps,sig,alp,chi):
    lnvlratio = f_lnvlratio(eps,sig)
    temp = np.zeros([NN,n_int,ndim])
    H = H0 * lnvlratio
    beta_c = bmax * ( np.log(Be/Ga) - lnvlratio*np.log(Be/De) ) / np.log(De/Ga)
    A = A0 * lnvlratio
    for i in range(NN):
        temp[i,i,0] = - alp[NN,1] * ((i+1)*rN)*beta_c*NN
        temp[i,i,1] = NN
        temp[i,NN,1] = A * ( 1 - alp[NN,1] * np.sign(alp[NN,1]) )
    return temp


# Definition of the yield surface(s)
def y(eps,sig,alp,chi):
    temp=np.zeros([n_y])
    MN_pq = f_MN_pq(sig)
    pc = f_pc(eps,sig,alp)
    CHI = f_CHI(eps,sig,alp,chi)
    for i in range(NN):
        temp[i] = CHI[i]**2 / ( ((i+1)*rN)**2 * MN_pq )     # Matsuoka-Nakai part and no consolidation part
        temp[i] -= 1
    return temp
def dyde(eps,sig,alp,chi):
    tempy=np.zeros([n_y,ndim])
    MN_pq = f_MN_pq(sig)
    pc = f_pc(eps,sig,alp)
    dpc = d_f_pc_de(eps,sig,alp)
    CHI = f_CHI(eps,sig,alp,chi)
    dCHI = d_f_CHI_de(eps,sig,alp,chi)
    for i in range(NN):
        tempy[i,0] = 2*dCHI[i,0]*CHI[i] / ( ((i+1)*rN)**2 * MN_pq )         # Matsuoka-Nakai part
    return tempy
def dyds(eps,sig,alp,chi):
    tempy=np.zeros([n_y,ndim])
    MN_pq = f_MN_pq(sig)
    dMN_pq = d_f_MN_pq_ds(sig)
    pc = f_pc(eps,sig,alp)
    dpc = d_f_pc_ds(eps,sig,alp)
    CHI = f_CHI(eps,sig,alp,chi)
    dCHI = d_f_CHI_ds(eps,sig,alp,chi)
    for i in range(NN):
        u = CHI[i]**2
        du_dp = 2*dCHI[i,0]* CHI[i]
        du_dq = 2*dCHI[i,1]* CHI[i]
        v = ((i+1)*rN)**2 * MN_pq
        dv_dp = ((i+1)*rN)**2 * dMN_pq[0]
        dv_dq = ((i+1)*rN)**2 * dMN_pq[1]
        tempy[i,0] = (du_dp*v-u*dv_dp) / v**2                               # Matsuoka-Nakai part
        tempy[i,1] = (du_dq*v-u*dv_dq) / v**2                               # Matsuoka-Nakai part
    return tempy
def dyda(eps,sig,alp,chi):
    tempy = np.zeros([n_y,n_int,ndim])
    MN_pq = f_MN_pq(sig)
    pc = f_pc(eps,sig,alp)
    dpc = d_f_pc_da(eps,sig,alp)
    CHI = f_CHI(eps,sig,alp,chi)
    dCHI = d_f_CHI_da(eps,sig,alp,chi)
    for i in range(NN):
        for j in range(n_int):
            tempy[i,j] = 2*dCHI[i,j] *CHI[i] / ( ((i+1)*rN)**2 * MN_pq )        # Matsuoka-Nakai part
    return tempy
def dydc(eps,sig,alp,chi):
    tempy = np.zeros([n_y,n_int,ndim])
    MN_pq = f_MN_pq(sig)
    pc = f_pc(eps,sig,alp)
    CHI = f_CHI(eps,sig,alp,chi)
    dCHI = d_f_CHI_dc(eps,sig,alp,chi)
    for i in range(NN):
        for j in range(n_int):
            tempy[i,j] = 2*dCHI[i,j] *CHI[i] / ( ((i+1)*rN)**2 * MN_pq )    # Matsuoka-Nakai part
    return tempy