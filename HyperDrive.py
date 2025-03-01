#
# HyperDrive
# Version 2.0.0
# (c) G.T. Houlsby, 2018-2021
#
# Routines to run the Hyperdrive system for implementing hyperplasticity
# models. Load this software then run "drive()" or "check()".
#

# set following line to auto = False if autograd not available
auto = True

if auto:
    import autograd as ag
    import autograd.numpy as np
else:
    import numpy as np

import copy
import importlib
import matplotlib.pyplot as plt
import os.path
import re
from scipy import optimize
import sys
import time

global rwt # reciprocal weights in internal variables * generalised stress product
global rec # record of test

global tfac, t5th, terr, use_RKDP
tfac = np.array([[   1.0/5.0,     0.0,        0.0,         0.0,        0.0      ],
                 [   3.0/40.0,    9.0/40.0,   0.0,         0.0,        0.0      ],
                 [   3.0/10.0,   -9.0/10.0,   6.0/5.0,     0.0,        0.0      ],
                 [ 226.0/729.0, -25.0/27.0, 880.0/729.0,  55.0/729.0,  0.0      ],
                 [-181.0/270.0,   5.0/2.0, -226.0/297.0, -91.0/27.0, 189.0/55.0]])
t5th = np.array([19.0/216.0, 0.0, 1000.0/2079.0, -125.0/216.0, 81.0/88.0,  5.0/56.0 ])
terr = np.array([11.0/360.0, 0.0,  -10.0/63.0,     55.0/72.0, -27.0/40.0, 11.0/280.0])
use_RKDP = False

S_txl_d = np.array([[0.0, -1.0,  1.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  1.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  1.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  1.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 1.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0]])
E_txl_d = np.array([[0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [1.0,  0.0,  0.0,  0.0,  0.0, 0.0]])
S_txl_u = np.array([[0.0, -1.0,  1.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  1.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  1.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 1.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0]])
E_txl_u = np.array([[0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [1.0,  1.0,  1.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [1.0,  0.0,  0.0,  0.0,  0.0, 0.0]])
S_dss   = np.array([[0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  1.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 1.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0]])
E_dss   = np.array([[1.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  1.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  1.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  0.0,  0.0, 0.0],
                    [0.0,  0.0,  0.0,  1.0,  0.0, 0.0]])

global quiet # option on printing at runtime
quiet = False
def qprint(*args):
    if not quiet: print(*args)
    
global voigt # option to use Voigt notation for input and output
             # internal calculations use Mandel vectors
voigt = False
def sig_from_voigt(sigi): 
    if voigt:
        temp = copy.deepcopy(sigi)
        temp[3:6] = temp[3:6]*Utils.root2
        return temp
    else:
        return sigi
def sig_to_voigt(sigo): 
    if voigt:
        temp = copy.deepcopy(sigo)
        temp[3:6] = temp[3:6]*Utils.rooth
        return temp
    else:
        return sigo
def eps_from_voigt(epsi): 
    if voigt:
        temp = copy.deepcopy(epsi)
        temp[3:6] = temp[3:6]*Utils.rooth
        return temp
    else:
        return epsi
def eps_to_voigt(epso): 
    if voigt:
        temp = copy.deepcopy(epso)
        temp[3:6] = temp[3:6]*Utils.root2
        return temp
    else:
        return epso

def read_def(text, default):
    temp = input(text + " [" + default + "]: ",)
    if len(temp) > 0: 
        return temp
    else:
        return default

global auto_t_step # option for setting automatic timestep
def substeps(nsub, tinc):
    if auto_t_step:
        return int(tinc / at_step) + 1
    else:
        return nsub

udef = "undefined"
def undef(arg1=[],arg2=[],arg3=[],arg4=[]): return udef

def error(text = "Unspecified error"):
    print(text)
    sys.exit()
def pause(args=[]):
    if len(args) > 0:
        message = " ".join(args)
        if len(message) > 0: print(message)
    text = input("Processing paused: hit ENTER to continue (x to exit)... ")
    if text == "x" or text == "X": sys.exit()

#printing routine for matrices
def pprint(x, label, n1=14, n2=8):
    def ftext(x,n1=14,n2=8):
        form = "{:"+str(n1)+"."+str(n2)+"g}"
        return form.format(x)
    if type(x) == str:
        print(label,x)
        return
    if hasattr(x,"shape"):
        if len(x.shape) == 0:
            print(label+ftext(x,n1,n2))
        elif len(x.shape) == 1:
            text = label+" ["
            for i in range(x.shape[0]):
                text = text+ftext(x[i],n1,n2)
            text = text+"]"
            print(text)
        elif len(x.shape) == 2:
            for i in range(x.shape[0]):
                if i == 0:
                    text = label+" ["
                    lenstart = len(text)
                else:
                    text = " " * lenstart
                text = text+"["
                for j in range(x.shape[1]):
                    text = text+ftext(x[i,j],n1,n2)
                text = text+"]"
                if i == x.shape[0]-1:
                    text = text+"]"                
                print(text)
        elif len(x.shape) == 3:
            for i in range(x.shape[0]):
                for j in range(x.shape[1]):
                    if i == 0 and j == 0:
                        text = label+" ["
                        lenstart = len(text)
                    else:
                        text = " " * lenstart
                    if j == 0:
                        text = text+"["
                    else:
                        text = text+" "
                    text = text+"["
                    for k in range(x.shape[2]):
                        text = text+ftext(x[i,j,k],n1,n2)
                    text = text+"]"
                    if j == x.shape[1]-1:
                        text = text+"]"                
                        if i == x.shape[0]-1:
                            text = text+"]"                
                    print(text)
        elif len(x.shape) == 4:
            for i in range(x.shape[0]):
                for j in range(x.shape[1]):
                    for k in range(x.shape[2]):
                        if i == 0 and j == 0 and k == 0:
                            text = label+" ["
                            lenstart = len(text)
                        else:
                            text = " " * lenstart
                        if j == 0 and k == 0:
                            text = text+"["
                        else:
                            text = text+" "
                        if k == 0:
                            text = text+"["
                        else:
                            text = text+" "
                        text = text+"["
                        for l in range(x.shape[3]):
                            text = text+ftext(x[i,j,k,l],n1,n2)
                        text = text+"]"
                        if k == x.shape[2]-1:
                            text = text+"]"      
                            if j == x.shape[1]-1:
                                text = text+"]"                
                                if i == x.shape[0]-1:
                                    text = text+"]"                
                        print(text)
        else:
            print(label)
            print(x)
    else:
        print(label+ftext(x,n1,n2))
    return

#Utility routines
class Utils:
    big   = 1.0e12
    small = 1.0e-12

    def mac(x): # Macaulay bracket
        if x <= 0.0: return 0.0
        else: return x
    def macm(x, delta = 0.0001): # Macaulay bracket, with possible rounding
        if x <= -delta: return 0.0
        if x >= delta: return x
        else: return ((x + delta)**2) / (4.0*delta)
    def S(x): # modified Signum function
        if abs(x) < Utils.small: return x / Utils.small
        elif x > 0.0: return 1.0  
        else: return -1.0
    def non_zero(x): # return a finite small value if argument close to zero
        if np.abs(x) < Utils.small: return Utils.small
        else: return x
    def Ineg(x): # Indicator function for set of negative reals
        if x <= 0.0: return 0.0 
        else: return Utils.big * 0.5*(x**2)
    def Nneg(x): # Normal Cone for set of negative reals
        if x <= 0.0: return 0.0 
        else: return Utils.big * x
    def w_rate_lin(y, mu): # convert canonical y to linear w function
        return (Utils.mac(y)**2) / (2.0*mu) 
    def w_rate_lind(y, mu): # differential of w_rate_lin
        return Utils.mac(y) / mu
    def w_rate_rpt(y, mu,r): # convert canonical y to Rate Process Theory w function
        return mu*(r**2) * (np.cosh(Utils.mac(y) / (mu*r)) - 1.0)
    def w_rate_rptd(y, mu,r): # differential of w_rate_rpt
        return r * np.sinh(Utils.mac(y) / (mu*r))

# Tensor utilities for mode 2 3x3 tensors
    # basic definitions
    delta = np.eye(3) # Kronecker delta, unit tensor
    def dev(t): # deviator
        return t - (Utils.tr1(t)/3.0)*Utils.delta
    def trans(t): # transpose
        return np.einsum("ij->ji",t)
    def sym(t): # symmetric part
        return (t + Utils.trans(t))/2.0
    def skew(t): # skew (antisymmetric) part
        return (t - Utils.trans(t))/2.0
    # shorthand for various products
    def cont(t1,t2): # double contraction of two tensors
        return np.einsum("ij,ij->",t1,t2)
    def iprod(t1,t2): # inner product of two tensors
        return np.einsum("ij,jk->ik",t1,t2)
    # 4th order products
    def pijkl(t1,t2): return np.einsum("ij,kl->ijkl",t1,t2)
    def pikjl(t1,t2): return np.einsum("ik,jl->ijkl",t1,t2)
    def piljk(t1,t2): return np.einsum("il,jk->ijkl",t1,t2)
    # 4th order unit and projection tensors
    II    = pikjl(delta, delta)  # 4th order unit tensor, contracts with t to give t
    IIb   = piljk(delta, delta)  # 4th order unit tensor, contracts with t to give t-transpose
    IIbb  = pijkl(delta, delta)  # 4th order unit tensor, contracts with t to give tr(t).delta
    IIsym = (II + IIb) / 2.0     # 4th order unit tensor (symmetric)
    PP    = II    - (IIbb / 3.0) # 4th order projection tensor, contracts with t to give dev(t)
                                 # also equal to differtential d_dev(t) / dt
    PPb   = IIb   - (IIbb / 3.0) # 4th order projection tensor, contracts with t to give dev(t)-transpose
    PPsym = IIsym - (IIbb / 3.0) # 4th order projection tensor (symmetric)
    # traces
    # _e variants are einsum versions for checking
    def tr1_e(t): return np.einsum("ii->",t)
    def tr2_e(t): return np.einsum("ij,ji->",t,t)
    def tr3_e(t): return np.einsum("ij,jk,ki->",t,t,t) 
    # mixed invariants of two tensors
    def trm_ab(a,b): return np.einsum("ij,ji->",a,b)
    def trm_a2b(a,b): return np.einsum("ij,jk,ki->",a,a,b)
    def trm_ab2(a,b): return np.einsum("ij,jk,ki->",a,b,b)
    def trm_a2b2(a,b): return np.einsum("ij,jk,kl,li->",a,a,b,b)
    # faster versions
    def tr1(t): return t[0,0] + t[1,1] + t[2,2] # trace
    def tr2(t): # trace of square
        return t[0,0]*t[0,0] + t[0,1]*t[1,0] + t[0,2]*t[2,0] + \
               t[1,0]*t[0,1] + t[1,1]*t[1,1] + t[1,2]*t[2,1] + \
               t[2,0]*t[0,2] + t[2,1]*t[1,2] + t[2,2]*t[2,2]
    def tr3(t): # trace of cube
        return t[0,0]*(t[0,0]*t[0,0] + t[0,1]*t[1,0] + t[0,2]*t[2,0]) + \
               t[0,1]*(t[1,0]*t[0,0] + t[1,1]*t[1,0] + t[1,2]*t[2,0]) + \
               t[0,2]*(t[2,0]*t[0,0] + t[2,1]*t[1,0] + t[2,2]*t[2,0]) + \
               t[1,0]*(t[0,0]*t[0,1] + t[0,1]*t[1,1] + t[0,2]*t[2,1]) + \
               t[1,1]*(t[1,0]*t[0,1] + t[1,1]*t[1,1] + t[1,2]*t[2,1]) + \
               t[1,2]*(t[2,0]*t[0,1] + t[2,1]*t[1,1] + t[2,2]*t[2,1]) + \
               t[2,0]*(t[0,0]*t[0,2] + t[0,1]*t[1,2] + t[0,2]*t[2,2]) + \
               t[2,1]*(t[1,0]*t[0,2] + t[1,1]*t[1,2] + t[1,2]*t[2,2]) + \
               t[2,2]*(t[2,0]*t[0,2] + t[2,1]*t[1,2] + t[2,2]*t[2,2]) 
    # invariants - use basic definitions
    def i1(t): return Utils.tr1(t) # 1st invariant
    def i2(t): # 2nd invariant (NB some sources define this with opposite sign)
        return (Utils.tr2(t) - Utils.tr1(t)**2)/2.0 
    def i3(t): # 3rd invariant
        return (2.0*Utils.tr3(t) - 3.0*Utils.tr2(t)*Utils.tr1(t) + Utils.tr1(t)**3)/6.0 
    def j2(t): # 2nd invariant of deviator
        return (3.0*Utils.tr2(t) - Utils.tr1(t)**2)/6.0 
    def j2_a(t): # 2nd invariant of deviator, alternative form for checking
        return Utils.i2(Utils.dev(t)) 
    def j3(t): # 3rd invariant of deviator
        return (9.0*Utils.tr3(t) - 9.0*Utils.tr2(t)*Utils.tr1(t) + 2.0*Utils.tr1(t)**3)/27.0 
    def j3_a(t): # 3rd invariant of deviator, alternative form for checking
        return Utils.i3(Utils.dev(t)) 
    def det(t): # determinant - should be same as 3rd invariant
        return t[0,0]*(t[1,1]*t[2,2] - t[1,2]*t[2,1]) + \
               t[0,1]*(t[1,2]*t[2,0] - t[1,0]*t[2,2]) + \
               t[0,2]*(t[1,0]*t[2,1] - t[1,1]*t[2,0])
    def dtr1(t): return Utils.delta # differential of trace
    def di1(t): return Utils.delta  # differential of 1st invariant

    # Voigt and Mandel notation conversion utilities
    # _ve for Voigt strain-like, _vs for Voigt stress-like,
    # _m for mandel (applicable to both)
    root2 = np.sqrt(2.0)
    rooth = np.sqrt(0.5)
    delta_m = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    def t_to_m(t): return np.array([t[0,0], t[1,1], t[2,2],
                                    Utils.root2*t[1,2], Utils.root2*t[2,0], Utils.root2*t[0,1]])
    def t_to_ve(t): return np.array([t[0,0], t[1,1], t[2,2],
                                     2.0*t[1,2], 2.0*t[2,0], 2.0*t[0,1]])
    def t_to_vs(t): return np.array([t[0,0], t[1,1], t[2,2],
                                     t[1,2], t[2,0], t[0,1]])
    def m_to_t(t): return np.array([[             t[0], Utils.rooth*t[5], Utils.rooth*t[4]], 
                                     [Utils.rooth*t[5],             t[1], Utils.rooth*t[3]], 
                                     [Utils.rooth*t[4], Utils.rooth*t[3],             t[2]]])
    def ve_to_t(t): return np.array([[    t[0], 0.5*t[5], 0.5*t[4]], 
                                     [0.5*t[5],     t[1], 0.5*t[3]], 
                                     [0.5*t[4], 0.5*t[3],     t[2]]])
    def vs_to_t(t): return np.array([[t[0], t[5], t[4]], 
                                     [t[5], t[1], t[3]], 
                                     [t[4], t[3], t[2]]])
    def m_to_ve(tm): return np.array([tm[0], tm[1], tm[2],
                                      Utils.root2*tm[3], Utils.root2*tm[4], Utils.root2*tm[5]])
    def m_to_vs(tm): return np.array([tm[0], tm[1], tm[2],
                                      Utils.rooth*tm[3], Utils.rooth*tm[4], Utils.rooth*tm[5]])
    def ve_to_m(ve): return np.array([ve[0], ve[1], ve[2],
                                      Utils.rooth*ve[3], Utils.rooth*ve[4], Utils.rooth*ve[5]])
    def vs_to_m(vs): return np.array([vs[0], vs[1], vs[2],
                                      Utils.root2*vs[3], Utils.root2*vs[4], Utils.root2*vs[5]])
    def dev_m(t): # deviator
        return t - (Utils.tr1_m(t)/3.0)*Utils.delta_m
    def cont_m(t1,t2): # contraction
        return np.einsum("i,i->",t1,t2) # dyadic product
    def pij_m(t1,t2): return np.einsum("i,j->ij",t1,t2)
    II_m = pij_m(delta_m,delta_m)
    PP_m = np.eye(6) - II_m/3.0 # projection tensor, also differential d_dev(t) / dt
    # traces
    def tr1_m(t): return t[0] + t[1] + t[2] # Mandel trace
    def tr2_m(t): return t[0]**2 + t[1]**2 + t[2]**2 + \
                         t[3]**2 + t[4]**2 + t[5]**2     # Mandel trace of square
    def tr3_m(t): # Mandel trace of cube
        return t[0]**3 + t[1]**3 + t[2]**3 + \
                1.5*(t[0]*(t[4]**2 + t[5]**2) + t[1]*(t[5]**2 + t[3]**2) + t[2]*(t[3]**2 + t[4]**2)) + \
                3.0*Utils.rooth*t[3]*t[4]*t[5] 
    # invariants
    def i1_m(t): # 1st invariant
        return Utils.tr1_m(t) 
    def i2_m(t): # 2nd invariant
        return (Utils.tr2_m(t) - Utils.tr1_m(t)**2)/2.0 
    def i3_m(t): # 3rd invariant
        return (2.0*Utils.tr3_m(t) - 3.0*Utils.tr2_m(t)*Utils.tr1_m(t) + Utils.tr1_m(t)**3)/6.0 
    def j2_m(t): # 2nd invariant of deviator
        return (3.0*Utils.tr2_m(t) - Utils.tr1_m(t)**2)/6.0 
    def j3_m(t): # 3rd invariant of deviator
        return (9.0*Utils.tr3_m(t) - 9.0*Utils.tr2_m(t)*Utils.tr1_m(t) + 2.0*Utils.tr1_m(t)**3)/27.0 
    def det_m(t): # determinant -should be equal to i3
        return t[0]*t[1]*t[2] + Utils.rooth*t[3]*t[4]*t[5] \
               - 0.5*(t[0]*t[3]*t[3] + t[1]*t[4]*t[4] + t[2]*t[5]*t[5])
                 
    #mixed invariants
    def trm_ab_m(a,b): 
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4] + a[5]*b[5]
    def trm_a2b_m(a,b): 
        return b[0]*(a[0]**2 + 0.5*(a[5]**2 + a[4]**2)) + \
               b[1]*(a[1]**2 + 0.5*(a[3]**2 + a[5]**2)) + \
               b[2]*(a[2]**2 + 0.5*(a[4]**2 + a[3]**2)) + \
               b[3]*(a[1]*a[3] + a[2]*a[3] + Utils.rooth*a[4]*a[5]) + \
               b[4]*(a[2]*a[4] + a[0]*a[4] + Utils.rooth*a[5]*a[3]) + \
               b[5]*(a[0]*a[5] + a[1]*a[5] + Utils.rooth*a[3]*a[4])
    def trm_ab2_m(a,b): 
        return a[0]*(b[0]**2 + 0.5*(b[5]**2 + b[4]**2)) + \
               a[1]*(b[1]**2 + 0.5*(b[3]**2 + b[5]**2)) + \
               a[2]*(b[2]**2 + 0.5*(b[4]**2 + b[3]**2)) + \
               a[3]*(b[1]*b[3] + b[2]*b[3] + Utils.rooth*b[4]*b[5]) + \
               a[4]*(b[2]*b[4] + b[0]*b[4] + Utils.rooth*b[5]*b[3]) + \
               a[5]*(b[0]*b[5] + b[1]*b[5] + Utils.rooth*b[3]*b[4])
#    def trm_a2b2_m(a,b): to be defined when needed
    # derivatives of invariants
    def di1_m(t): return Utils.delta_m
    def di2_m(t): return t - Utils.di1_m(t)*Utils.tr1_m(t)
#   def di3_m(t): to be defined when needed
    def di3_m(t):
        return np.array([t[1]*t[2] - 0.5*t[3]*t[3],
                         t[2]*t[0] - 0.5*t[4]*t[4],
                         t[0]*t[1] - 0.5*t[5]*t[5],
                         Utils.rooth*t[4]*t[5] - t[3]*t[0],
                         Utils.rooth*t[5]*t[3] - t[4]*t[1],
                         Utils.rooth*t[3]*t[4] - t[5]*t[2]])
    def dj2_m(t): return t - Utils.tr1_m(t)*Utils.di1_m(t)/3.0
#   def dj3_m(t): to be defined when needed
    #differentials of mixed invariants - add others when needed
    def dtrm_ab_a_m(a,b): # d_tr(ab) / da
        return b
    def dtrm_ab_b_m(a,b): # d_tr(ab) / db
        return a
    def dtrm_a2b_a_m(a,b): # d_tr(aab) / da
        return np.array([2.0*a[0]*b[0] + a[4]*b[4] + a[5]*b[5],
                         2.0*a[1]*b[1] + a[5]*b[5] + a[3]*b[3],
                         2.0*a[2]*b[2] + a[3]*b[3] + a[4]*b[4],
                         (a[1] + a[2])*b[3] + a[3]*(b[1] + b[2]) + Utils.rooth*(a[4]*b[5] + a[5]*b[4]),
                         (a[2] + a[0])*b[4] + a[4]*(b[2] + b[0]) + Utils.rooth*(a[5]*b[3] + a[3]*b[5]),
                         (a[0] + a[1])*b[5] + a[5]*(b[0] + b[1]) + Utils.rooth*(a[3]*b[4] + a[4]*b[3])])
    def dtrm_a2b_b_m(a,b): # d_tr(aab) / db
        return np.array([a[0]**2 + 0.5*(a[5]**2 + a[4]**2),          \
                         a[1]**2 + 0.5*(a[3]**2 + a[5]**2),          \
                         a[2]**2 + 0.5*(a[4]**2 + a[3]**2),          \
                         (a[1] + a[2])*a[3] + Utils.rooth*a[4]*a[5], \
                         (a[2] + a[0])*a[4] + Utils.rooth*a[5]*a[3], \
                         (a[0] + a[1])*a[5] + Utils.rooth*a[3]*a[4]])
    def dtrm_ab2_a_m(a,b): # d_tr(abb) / da
        return np.array([b[0]**2 + 0.5*(b[5]**2 + b[4]**2),          \
                         b[1]**2 + 0.5*(b[3]**2 + b[5]**2),          \
                         b[2]**2 + 0.5*(b[4]**2 + b[3]**2),          \
                         (b[1] + b[2])*b[3] + Utils.rooth*b[4]*b[5], \
                         (b[2] + b[0])*b[4] + Utils.rooth*b[5]*b[3], \
                         (b[0] + b[1])*b[5] + Utils.rooth*b[3]*b[4]])
    def dtrm_ab2_b_m(a,b): # d_tr(abb) / db
        return np.array([2.0*b[0]*a[0] + b[4]*a[4] + b[5]*a[5],
                         2.0*b[1]*a[1] + b[5]*a[5] + b[3]*a[3],
                         2.0*b[2]*a[2] + b[3]*a[3] + b[4]*a[4],
                         (a[1] + a[2])*b[3] + a[3]*(b[1] + b[2]) + Utils.rooth*(a[4]*b[5] + a[5]*b[4]),
                         (a[2] + a[0])*b[4] + a[4]*(b[2] + b[0]) + Utils.rooth*(a[5]*b[3] + a[3]*b[5]),
                         (a[0] + a[1])*b[5] + a[5]*(b[0] + b[1]) + Utils.rooth*(a[3]*b[4] + a[4]*b[3])])
    # second differential - define others when needed
    def d2j2_m(t): return Utils.PP_m
    # signumof deviator and its differential
    def S_j2_m(t):
        if Utils.j2_m(t) == 0.0: 
            temp = np.zeros([6])
        else:
            temp = Utils.dev_m(t) / np.sqrt(2.0*Utils.j2_m(t))
        return temp
    def dS_j2_m(t):
        if Utils.j2_m(t) == 0.0: 
            temp = np.zeros([6,6])
        else:
            temp  =  Utils.PP_m / np.sqrt(2.0*Utils.j2_m(t))
            temp += -Utils.pij_m(Utils.dev_m(t),Utils.dev_m(t)) / ((2.0*Utils.j2_m(t))**1.5)
        return temp

def process(source): # processes the Hyperdrive commands
    global line # necessary to allow recording of history using histrec
    for line in source:
        text = line.rstrip("\n\r")
        textsplit = re.split(r'[ ,;]',text)
        keyword = textsplit[0]
        if len(keyword) == 0:
            continue
        if keyword[:1] == "#":
            continue
        elif keyword[:1] == "*" and hasattr(Commands, keyword[1:]):
            getattr(Commands, keyword[1:])(textsplit[1:])
            if keyword == "*end": break
        # else:
        #     pause("\x1b[0;31mWARNING - keyword not recognised:",keyword,"\x1b[0m")

class Commands:
# each of these functions corresponds to a command in Hyperdrive
    def title(args): # read a job title
        global title
        title = " ".join(args)
        print("Title: ", title)
    def mode(args): # set mode (0 = scalar, 1 = vector, 2 = tensor)
        global mode, n_dim
        mode = int(args[0])
        if mode == 0: 
            n_dim = 1
        elif mode == 1: 
            n_dim = int(args[1])
        elif mode == 2: 
            n_dim = int(args[1])
        else:
            error("Unrecognised mode:"+mode)
        print("Mode:", mode, "n_dim:", n_dim)
    def voigt(args):
        print("Setting Voigt mode (mode = 1, n_dim = 6)")
        global mode, n_dim, voigt
        mode = 1
        n_dim = 6
        voigt = True
    def model(args): # specify model
        global mode, n_dim
        global hm, model
        model_temp = args[0]
        print("Current directory:",os.getcwd())
        print("Current path directory:",os.path.curdir)
        if os.path.isfile("./" + model_temp + ".py"):
            model = model_temp + ""
            print("Importing hyperplasticity model: ", model)
            hm = importlib.import_module(model)
            if hasattr(hm, "deriv"): hm.deriv()
            princon(hm)
            print("Description: ", hm.name)
            if hm.mode != mode:
                print("Model mode mismatch: hm.mode:",hm.mode,", mode:", mode)
                error()
            if mode == 0:
                n_dim = 1
            else:
                if hm.ndim != n_dim:
                    print("Model n_dim mismatch: hm.ndim:",hm.ndim,", n_dim:", n_dim)
                    error()
        else:
            error("Model not found:" + model_temp)
        for fun in ["f", "g", "y", "w"]:
            if not hasattr(hm, fun): print("Function",fun,"not present in",model)
        set_up_analytical()
        set_up_auto()
        set_up_num()
        choose_diffs()
        
# commands for setting options
    def prefs(args=[]): # set preferences for method of differentiation
        global pref
        print("Setting differential preferences:")
        pref[0:3] = args[0:3]
        for i in range(3):
            print("Preference", i+1, pref[i])
        choose_diffs()
    def f_form(args=[]): # use f-functions for preference
        global fform, gform
        print("Setting f-form")
        fform = True
        gform = False
    def g_form(args=[]): # use g-functions for preference
        global gform, fform
        print("Setting g-form")
        gform = True
        fform = False
    def rate(args=[]): # use rate-dependent algorithm
        global rate, hm
        print("Setting rate dependent analysis")
        rate = True
        if len(args) > 0: hm.mu = float(args[0])
    def rateind(args=[]): # use rate-independent algorithm
        global rate
        print("Setting rate independent analysis")
        rate = False
    def acc(args): # set acceleration factor for rate-independent algorithm
        global acc
        acc = float(args[0])
        print("Setting acceleration factor:", acc)
    def RKDP(args): # set RKDP for general increment
        global use_RKDP
        use_RKDP = True
        print("Setting RKDP mode")
    def no_RKDP(args): # unset RKDP for general increment (default)
        global use_RKDP
        use_RKDP = False
        print("Unsetting RKDP mode")
    def quiet(args=[]): # set quiet mode
        global quiet
        quiet = True
        print("Setting quiet mode")
    def unquiet(args=[]): # unset quiet mode (default)
        global quiet
        quiet = False
        print("Unsetting quiet mode")
    def auto_t_step(args=[]): # set auto_t_step mode
        global auto_t_step, at_step
        auto_t_step = True
        at_step = float(args[0])
        print("Setting auto_t_step mode")
    def no_auto_t_step(args=[]): # unset auto_t_step (default)
        global auto_t_step
        auto_t_step = False 
        print("Unsetting auto_t_step mode")
    def colour(args): # select plot colour
        global colour
        colour = args[0]
        print("Setting plot colour:", colour)

# commands for setting constants
    def const(args): # read model constants
        global model, hm
        if model != "undefined":
            hm.const = [float(tex) for tex in args]
            qprint("Constants:", hm.const)
            if hasattr(hm, "deriv"): hm.deriv()
            princon(hm)
        else:
            pause("Cannot read constants: model undefined")
    def tweak(args): # tweak a single model constant
        global model, hm
        if model != "undefined":
            const_tweak = args[0]
            val_tweak = float(args[1])
            for i in range(len(hm.const)):
                if const_tweak == hm.name_const[i]: 
                    hm.const[i] = val_tweak 
                    print("Tweaked constant: "+hm.name_const[i]+", value set to", hm.const[i])
            if hasattr(hm, "deriv"): hm.deriv()
        else:
            pause("Cannot tweak constant: model undefined")
    def const_from_points(args):
        global model, hm
        if model != "undefined":
            npt = int(args[0])
            epsd = np.zeros(npt)
            sigd = np.zeros(npt)
            for ipt in range(npt):
                epsd[ipt] = float(args[1+2*ipt])
                sigd[ipt] = float(args[2+2*ipt])
            Einf = float(args[1+2*npt])
            epsmax = float(args[2+2*npt])
            HARM_R = float(args[3+2*npt])
            hm.const = derive_from_points(model,epsd,sigd,Einf,epsmax,HARM_R)
            print("Constants from points:", hm.const)
            if hasattr(hm, "deriv"): hm.deriv()
            princon(hm)
        else:
            pause("Cannot obtain constants, model undefined")
    def const_from_curve(args):
        global model, hm
        modtype = args[0]
        curve = args[1]
        npt = int(args[2])
        sigmax = float(args[3])
        HARM_R = float(args[4])
        param = np.zeros(3)
        param[0:3] = [float(item) for item in args[5:8]]
        epsd = np.zeros(npt)
        sigd = np.zeros(npt)
        print("Calculated points from curve:")
        for ipt in range(npt):
            sigd[ipt] = sigmax*float(ipt+1)/float(npt)
            if curve == "power":
                Ei = param[0]
                epsmax = param[1]
                power = param[2]
                epsd[ipt] = sigd[ipt]/Ei + (epsmax-sigmax/Ei)*(sigd[ipt]/sigmax)**power
            if curve == "jeanjean":                    
                Ei = param[0]
                epsmax = param[1]
                A = param[2]
                epsd[ipt] = sigd[ipt]/Ei + (epsmax-sigmax/Ei)*((np.atanh(np.tanh(A)*sigd[ipt]/sigmax)/A)**2)
            if curve == "PISA":
                Ei = param[0]
                epsmax = param[1]
                n = param[2]
                epspmax = epsmax - sigmax/Ei
                A = n*sigmax/(Ei*epspmax)
                B = -2.0*A*sigd[ipt]/sigmax + (1.0-n)*((1.0 + sigmax/(Ei*epspmax))**2)*(sigd[ipt]/sigmax - 1.0)
                C = A*(sigd[ipt]/sigmax)**2
                D = np.max(B**2-4.0*A*C, 0.0)
                epsd[ipt] = sigd[ipt]/Ei + 2.0*epspmax*C/(-B + np.sqrt(D))
            print(epsd[ipt],sigd[ipt])
        Einf = 0.5*(sigd[npt-1]-sigd[npt-2]) / (epsd[npt-1]-epsd[npt-2])
        epsmax = epsd[npt-1]
        hm.const = derive_from_points(modtype, epsd, sigd, Einf, epsmax, HARM_R)
        print("Constants from curve:", hm.const)
        if hasattr(hm, "deriv"): hm.deriv()
        princon(hm)
    def const_from_data(args): #not yet implemented
        global model, hm
        if model != "undefined":
            dataname = args[0]
            if dataname[-4:] != ".csv": dataname = dataname + ".csv"
            print("Reading data from", dataname)
            data_file = open(dataname,"r")
            data_text = data_file.readlines()
            data_file.close()
            for i in range(len(data_text)):
                data_split = re.split(r'[ ,;]',data_text[i])
#                ttest = float(data_split[0])
                epstest = float(data_split[1])
                sigtest = float(data_split[2])
            npt = int(args[1])
            maxsig = float(args[2])
            HARM_R = float(args[3])
            epsd = np.zeros(npt)
            sigd = np.zeros(npt)
            print("Calculating const from data")
            for ipt in range(npt):
                sigd[ipt] = maxsig*float(ipt+1)/float(npt)
                for i in range(len(data_text)-1):
                    if sigd[ipt] >= sigtest[i] and sigd[ipt] <= sigtest[i+1]:
                        epsd[ipt] = epstest[i] + (sigtest[i+1]-sigd[ipt]) * \
                                    (epstest[i+1]-epstest[i])/(sigtest[i+1]-sigtest[i])
                print(epsd[ipt],sigd[ipt])
            Einf = 0.5*(sigd[npt-1]-sigd[npt-2])/(epsd[npt-1]-epsd[npt-2])
            epsmax = epsd[npt-1]
            hm.const = derive_from_points(model,epsd,sigd,Einf,epsmax,HARM_R)
            print("Constants from data:", hm.const)
            if hasattr(hm, "deriv"): hm.deriv()
            princon(hm)
        else:
            pause("Cannot obtain constants, model undefined")

# commands for starting and stopping processing
    def init(args=[]):
        startup()
    def start(args=[]):
        global hm
        global n_stage, n_cyc, n_print
        global rec, curr_test, test_rec, test_col, colour
        global start, last, high, test_high
        global t
        if hasattr(hm, "deriv"): hm.deriv()
        modestart()
        if hasattr(hm, "setalp"): hm.setalp(alp)
        n_stage = 0
        n_cyc = 0
        n_print = 0
        t = 0.0
        rec=[[]]
        curr_test = 0
        test_rec=[]
        test_col=[]
        record(eps, sig)
        test_col.append(colour)
        start = time.process_time()
        last = start + 0.0
        high = [[]]
        test_high = []

    def restart(args=[]):
        global hm
        global n_stage, n_cyc, n_print
        global rec, test_rec, curr_test, test_col, colour, high
        global start, last
        global t
        t = 0.0
        if hasattr(hm, "deriv"): hm.deriv()
        rec.append([])
#new
        #test_rec.append([])
#old
#end change
        high.append([])
        curr_test += 1
        modestart()
        if hasattr(hm, "setalp"): hm.setalp(alp)
        n_stage = 0
        n_cyc = 0
        n_print = 0 
        record(eps, sig)
        test_col.append(colour)
        now = time.process_time()
        print("Time:",now-start,now-last)
        last = now
    def rec(args=[]):
        global recording
        print("Recording data")
        recording = True
    def stoprec(args=[]):
        global mode, recording
        if mode == 0: 
            temp = np.nan
        else: 
            temp = np.full(n_dim, np.nan)
        record(temp, temp) # write a line of nan to split plots
        print("Stop recording data")
        recording = False
    def end(args=[]):
        global start, last
        now = time.process_time()
        print("Time:",now-start,now-last)
        qprint("End of test")
        print("")

# commands for initialisation
    def init_stress(args):
        global sig, eps, chi
        histrec()
        sigi = readvs(args)
        sig = sig_from_voigt(sigi) # convert if necessary
        eps = eps_g(sig,alp)
        chi = chi_g(sig,alp)
        print("Initial stress:", sig)
        print("Initial strain:", eps)
        delrec()
        record(eps, sig)
    def init_strain(args):
        global sig, eps, chi
        histrec()
        epsi = readvs(args)
        eps = eps_from_voigt(epsi) # convert if necessary
        sig = sig_f(eps,alp)
        chi = chi_f(eps,alp)
        print("Initial strain:", eps)
        print("Initial stress:", sig)
        delrec()
        record(eps, sig)
    def save_state(args=[]):
        global eps, sig, alp, chi, save_eps, save_sig, save_alp, save_chi
        save_eps = copy.deepcopy(eps)
        save_sig = copy.deepcopy(sig)
        save_alp = copy.deepcopy(alp)
        save_chi = copy.deepcopy(chi)
        print("State saved")
    def restore_state(args=[]):
        global eps, sig, alp, chi, save_eps, save_sig, save_alp, save_chi
        eps = copy.deepcopy(save_eps)
        sig = copy.deepcopy(save_sig)
        alp = copy.deepcopy(save_alp)
        chi = copy.deepcopy(save_chi)
        delrec()
        record(eps, sig)
        print("State restored")
        
# commands for applying stress and strain increments
    def v_txl_d(args):
        histrec()
        Smat = copy.deepcopy(S_txl_d)
        Emat = copy.deepcopy(E_txl_d)
        Tdt = np.zeros(6)
        dt     = float(args[0])
        Tdt[5] = float(args[1])
        nprint =   int(args[2])
        nsub   =   int(args[3])
        qprint("Drained triaxial test:")
        qprint("deps11 =", Tdt[5])
        run_general_inc(Smat, Emat, Tdt, dt, nprint, nsub)
    def v_txl_u(args):
        histrec()
        Smat = copy.deepcopy(S_txl_u)
        Emat = copy.deepcopy(E_txl_u)
        Tdt = np.zeros(6)
        dt     = float(args[0])
        Tdt[5] = float(args[1])
        nprint =   int(args[2])
        nsub   =   int(args[3])
        qprint("Undrained triaxial test:")
        qprint("deps11 =", Tdt[5])
        run_general_inc(Smat, Emat, Tdt, dt, nprint, nsub)
    def v_dss(args):
        histrec()
        Smat = copy.deepcopy(S_dss)
        Emat = copy.deepcopy(E_dss)
        Tdt = np.zeros(6)
        dt     = float(args[0])
        Tdt[5] = float(args[1])
        nprint =   int(args[2])
        nsub   =   int(args[3])
        qprint("Direct shear test:")
        qprint("dgam23  =", Tdt[5])
        run_general_inc(Smat, Emat, Tdt, dt, nprint, nsub)
        
        
# commands for applying stress and strain increments
###############################################################################
###############################################################################
############################ ADDITIONAL FEATURE XXX ###########################
    def und_stress_from_strain(args):
        global sig, eps
        global iprint, isub
        global start_inc
        start_inc=True
        isub=0
        ccc=0
        histrec()
        t_inc = float(args[0])
        q_target = float(args[1])
        eps_inc = [0,float(args[2])]
        nsub = int(args[3])
        deps = eps_inc
        dt = t_inc
        # print("Strain increment:", eps_inc, "deps =", deps, "dt =", dt)
        if deps[1]>0:
            record(eps, sig)
            while sig[1]<q_target:
                if sig[0]<30:
                    ccc = 1
                if sig[0]<20:
                    ccc = 2
                if sig[0]<10:
                    ccc = 3
                de = deps[1] / 1**ccc 
                strain_inc([0,de], dt)
                isub +=1
                if isub == nsub:
                    record(eps, sig)
                    isub=0
        else:
            record(eps, sig)
            while sig[1]>q_target:
                strain_inc(deps, dt)
                isub +=1
                if isub == nsub:
                    record(eps, sig)
                    isub=0

    def dr_stress_trgt(args):
        global sig, eps
        global iprint, isub
        global start_inc
        start_inc=True
        isub=0
        ccc=0
        histrec()
        t_inc = float(args[0])
        q_target = float(args[1])
        eps_inc = [float(args[2])]
        nsub = int(args[3])
        deps = eps_inc[0]
        dt = t_inc
        # print("Strain increment:", eps_inc, "deps =", deps, "dt =", dt)
        if deps>0:
            record(eps, sig)
            while sig[1]<q_target:
                general_inc([[3,-1],[0,0]],[[0,0],[0,1]],[0.0,deps],1)
                isub +=1
                if isub == nsub:
                    record(eps, sig)
                    isub=0
        else:
            record(eps, sig)
            while sig[1]>q_target:
                general_inc([[3,-1],[0,0]],[[0,0],[0,1]],[0.0,deps],1)
                isub +=1
                if isub == nsub:
                    record(eps, sig)
                    isub=0

############################ ADDITIONAL FEATURE XXX ###########################
###############################################################################
###############################################################################
        
        
    def general_inc(args):
        histrec()
        Smat = np.reshape(np.array([float(i) for i in args[:n_dim*n_dim]]), (n_dim, n_dim))
        Emat = np.reshape(np.array([float(i) for i in args[n_dim*n_dim:2*n_dim*n_dim]]), (n_dim, n_dim))
        Tdt = np.array([float(i) for i in args[2*n_dim*n_dim:2*n_dim*n_dim+n_dim]])
        dt     = float(args[2*n_dim*n_dim+n_dim])
        nprint =   int(args[2*n_dim*n_dim+n_dim+1])
        nsub   =   int(args[2*n_dim*n_dim+n_dim+2])
        qprint("General control increment:")
        qprint("S   =", Smat)
        qprint("E   =", Emat)
        qprint("Tdt =", Tdt)
        run_general_inc(Smat, Emat, Tdt, dt, nprint, nsub)
    def general_cyc(args):
        global sig, eps
        global start_inc
        histrec()
        Smat = np.reshape(np.array([float(i) for i in args[:n_dim*n_dim]]), (n_dim, n_dim))
        Emat = np.reshape(np.array([float(i) for i in args[n_dim*n_dim:2*n_dim*n_dim]]), (n_dim, n_dim))
        Tdt = np.array([float(i) for i in args[2*n_dim*n_dim:2*n_dim*n_dim+n_dim]])
        tper   = float(args[2*n_dim*n_dim + n_dim])
        ctype  =       args[2*n_dim*n_dim + n_dim + 1]
        ncyc   =   int(args[2*n_dim*n_dim + n_dim + 2])
        nprint =   int(args[2*n_dim*n_dim + n_dim + 3])
        if nprint%2 == 1: nprint += 1
        nsub = int(args[2*n_dim*n_dim + n_dim + 4])
        qprint("General control cycles:")
        qprint("S   =", Smat)
        qprint("E   =", Emat)
        qprint("Tdt =", Tdt)
        start_inc = True
        if ctype == "saw":
            dTdt  = Tdt / float(nprint/2) / float(nsub)
            dtper = tper / float(nprint/2) / float(nsub)
            qprint("General cycles (saw): tper=", tper, "ncyc =", ncyc)
            for icyc in range(2*ncyc):
                for iprint in range(int(nprint/2)):
                    for isub in range(nsub): general_inc(Smat, Emat, dTdt, dtper)
                    record(eps, sig)
                dTdt = -dTdt
        if ctype == "sine":
            dtper = tper / float(nprint) / float(nsub)
            qprint("General cycles (sine): tper=", tper, "ncyc =", ncyc)
            for icyc in range(ncyc):
                for iprint in range(int(nprint)):
                    th1 = 2.0*np.pi*float(iprint)  /float(nprint)
                    th2 = 2.0*np.pi*float(iprint+1)/float(nprint)
                    dTdt = Tdt * (np.sin(th2) - np.sin(th1)) / float(nsub)
                    for isub in range(nsub): general_inc(Smat, Emat, dTdt, dtper)
                    record(eps, sig)
        if ctype == "haversine":
            dtper = tper / float(nprint) / float(nsub)
            qprint("General cycles (haversine): tper=", tper, "ncyc =", ncyc)
            for icyc in range(ncyc):
                for iprint in range(int(nprint)):
                    th1 = 2.0*np.pi*float(iprint)  /float(nprint)
                    th2 = 2.0*np.pi*float(iprint+1)/float(nprint)
                    dTdt = Tdt * (np.sin(th2/2.0)**2 - np.sin(th1/2.0)**2) / float(nsub)
                    for isub in range(nsub): general_inc(Smat, Emat, dTdt, dtper)
                    record(eps, sig)
        qprint("Cycles complete")
    def strain_inc(args):
        global sig, eps
        global start_inc
        histrec()
        t_inc = float(args[0])
        eps_vali = readv(args)
        eps_val = eps_from_voigt(eps_vali) # convert if necessary
        if mode == 2:
            nprint = int(args[n_dim**2+1])
            nsub   = int(args[n_dim**2+2])
        else:
            nprint = int(args[n_dim+1])
            nsub   = int(args[n_dim+2])
        nsub = substeps(nsub, t_inc/float(nprint))
        deps = eps_val / float(nprint*nsub)
        dt = t_inc / float(nprint*nsub)
        qprint("Strain increment:", nprint, "steps,", nsub, "substeps")
        qprint("eps_inc =", eps_val)
        qprint("deps =", deps)
        qprint("t_inc =", t_inc, ", dt =",dt)
        start_inc = True
        for iprint in range(nprint):
            for isub in range(nsub): 
                strain_inc(deps, dt)
            record(eps, sig)
        qprint("Increment complete")
    def strain_targ(args):
        global sig, eps
        global start_inc
        histrec()
        t_inc = float(args[0])
        eps_vali = readv(args)
        eps_val = eps_from_voigt(eps_vali) # convert if necessary
        if mode == 2:
            nprint = int(args[n_dim**2+1])
            nsub   = int(args[n_dim**2+2])
        else:
            nprint = int(args[n_dim+1])
            nsub   = int(args[n_dim+2])    
        deps = (eps_val - eps) / float(nprint*nsub)
        dt = t_inc / float(nprint*nsub)
        qprint("Strain target:", nprint, "steps,", nsub, "substeps")
        qprint("eps_targ =", eps_val)
        qprint("deps =", deps)
        qprint("t_inc =", t_inc, ", dt =",dt)
        start_inc = True
        for iprint in range(nprint):
            for isub in range(nsub): 
                strain_inc(deps, dt)
            record(eps, sig)
        qprint("Increment complete")
    def strain_cyc(args):
        global sig, eps
        global start_inc
        histrec()
        tper = float(args[0])
        eps_cyci = readv(args)
        eps_cyc = eps_from_voigt(eps_cyci) # convert if necessary
        ctype = args[n_dim+1]
        ncyc = int(args[n_dim+2])
        nprint = int(args[n_dim+3])
        if nprint%2 == 1: nprint += 1
        nsub = int(args[n_dim+4])
        dt = tper / float(nprint*nsub)
        start_inc = True
        if ctype == "saw":
            deps = eps_cyc / float(nprint/2) / float(nsub)
            qprint("Strain cycle (saw):", eps_cyc, "tper=", tper, "ncyc =", ncyc)
            for icyc in range(2*ncyc):
                for iprint in range(int(nprint/2)):
                    for isub in range(nsub): strain_inc(deps,dt)
                    record(eps, sig)
                deps = -deps
        if ctype == "sine":
            qprint("Strain cycle (sine):", eps_cyc, "tper=", tper, "ncyc =", ncyc)
            for icyc in range(ncyc):
                for iprint in range(int(nprint)):
                    th1 = 2.0*np.pi*float(iprint)/float(nprint)
                    th2 = 2.0*np.pi*float(iprint+1)/float(nprint)
                    deps = eps_cyc * (np.sin(th2) - np.sin(th1)) / float(nsub)
                    for isub in range(nsub): strain_inc(deps,dt)
                    record(eps, sig)
        if ctype == "haversine":
            qprint("Strain cycle (haversine):", eps_cyc, "tper=", tper, "ncyc =", ncyc)
            for icyc in range(ncyc):
                for iprint in range(int(nprint)):
                    th1 = 2.0*np.pi*float(iprint)/float(nprint)
                    th2 = 2.0*np.pi*float(iprint+1)/float(nprint)
                    deps = eps_cyc * (np.sin(th2/2.0)**2 - np.sin(th1/2.0)**2) / float(nsub)
                    for isub in range(nsub): strain_inc(deps,dt)
                    record(eps, sig)
        qprint("Cycles complete")
    def stress_inc(args):
        global sig, eps
        global start_inc
        histrec()
        sig_vali = readv(args)
        sig_val = sig_from_voigt(sig_vali) # convert if necessary
        if mode == 2:
            nprint = int(args[n_dim**2+1])
            nsub   = int(args[n_dim**2+2])
        else:
            nprint = int(args[n_dim+1])
            nsub   = int(args[n_dim+2])    
        t_inc = float(args[0])
        dsig = sig_val / float(nprint*nsub)
        dt = t_inc / float(nprint*nsub)
#        print("Stress increment:", "{:14.6f}".format(sig_val), 
#              "dsig =", "{:14.6f}".format(dsig), 
#              "dt =", "{:14.6f}".format(dt))
        qprint("Stress increment:", nprint, "steps,", nsub, "substeps")
        qprint("sig_inc =", sig_val)
        qprint("dsig =", dsig)
        qprint("t_inc =", t_inc, ", dt =",dt)
        start_inc = True
        for iprint in range(nprint):
            for isub in range(nsub): 
                stress_inc(dsig, dt)
            record(eps, sig)
        qprint("Increment complete")
    def stress_targ(args):
        global sig, eps
        global start_inc
        histrec()
        sig_vali = readv(args)
        sig_val = sig_from_voigt(sig_vali) # convert if necessary
        if mode == 2:
            nprint = int(args[n_dim**2+1])
            nsub   = int(args[n_dim**2+2])
        else:
            nprint = int(args[n_dim+1])
            nsub   = int(args[n_dim+2])    
        t_inc = float(args[0])
        dsig = (sig_val - sig) / float(nprint*nsub)
        dt = t_inc / float(nprint*nsub)
        qprint("Stress target:", nprint, "steps,", nsub, "substeps")
        qprint("sig_targ =", sig_val)
        qprint("dsig =", dsig)
        qprint("t_inc =", t_inc, ", dt =",dt)
#        print("Stress target:", "{:14.6f}".format(sig_val), 
#              "deps =", "{:14.6f}".format(dsig), 
#              "dt =", "{:14.6f}".format(dt))
        start_inc = True
        for iprint in range(nprint):
            for isub in range(nsub): 
                stress_inc(dsig, dt)
            record(eps, sig)
        qprint("Increment complete")
    def stress_cyc(args):
        global sig, eps
        global start_inc
        histrec()
        sig_cyci = readv(args)
        sig_cyc = sig_from_voigt(sig_cyci) # convert if necessary
        tper = float(args[0])
        ctype = args[n_dim+1]
        ncyc = int(args[n_dim+2])
        nprint = int(args[n_dim+3])
        if nprint%2 == 1: nprint += 1
        nsub = int(args[n_dim+4])
        dt = tper / float(nprint*nsub)
        start_inc = True
        if ctype == "saw":
            dsig = sig_cyc / float(nprint/2) / float(nsub)
            qprint("Stress cycle (saw):", sig_cyc, "tper=", tper, "ncyc =", ncyc)
            for icyc in range(2*ncyc):
                for iprint in range(nprint/2):
                    for isub in range(nsub): 
                        stress_inc(dsig,dt)
                    record(eps, sig)
                dsig = -dsig
        if ctype == "sine":
            qprint("Stress cycle (sine):", sig_cyc, "tper=", tper, "ncyc =", ncyc)
            for icyc in range(ncyc):
                for iprint in range(int(nprint)):
                    th1 = 2.0*np.pi*float(iprint)/float(nprint)
                    th2 = 2.0*np.pi*float(iprint+1)/float(nprint)
                    dsig = sig_cyc * (np.sin(th2) - np.sin(th1)) / float(nsub)
                    for isub in range(nsub): 
                        stress_inc(dsig,dt)
                    record(eps, sig)
        if ctype == "haversine":
            qprint("Stress cycle (haversine):", sig_cyc, "tper=", tper, "ncyc =", ncyc)
            for icyc in range(ncyc):
                for iprint in range(int(nprint)):
                    th1 = 2.0*np.pi*float(iprint)/float(nprint)
                    th2 = 2.0*np.pi*float(iprint+1)/float(nprint)
                    dsig = sig_cyc * (np.sin(th2/2.0)**2 - np.sin(th1/2.0)**2) / float(nsub)
                    for isub in range(nsub): 
                        stress_inc(dsig,dt)
                    record(eps, sig)
        qprint("Cycles complete")
    def test(args, ptype): # hidden command used by stress_test and strain_test
        global mode, n_dim, t, test, test_rec
        global sig, eps
        global start_inc
        histrec()
        test = True
        testname = args[0]
        if testname[-4:] != ".csv": testname = testname + ".csv"
        print("Reading from", testname)
        test_file = open(testname,"r")
        test_text = test_file.readlines()
        test_file.close()
        nsub = int(args[1])
#new
        #test_rec = []
#old
        test_rec = []
#end change
        start_inc = True
        for iprint in range(len(test_text)):
            ltext = test_text[iprint].replace(" ","")
            #print(ltext)
            test_split = re.split(r'[ ,;]',ltext)
            ttest = float(test_split[0])
            if mode==0:
                epstest = float(test_split[1])
                sigtest = float(test_split[2])
            else:                
                epstest = np.array([float(tex) for tex in test_split[1:n_dim+1]])
                sigtest = np.array([float(tex) for tex in test_split[n_dim+1:2*n_dim+1]])
            if start_inc:
                eps = eps_from_voigt(epstest)
                sig = sig_from_voigt(sigtest)
                delrec()
                record(eps, sig)
                recordt(epstest, sigtest)
                start_inc= False
            else:
                dt = (ttest - t) / float(nsub)
                if ptype == "strain":
                    deps = (eps_from_voigt(epstest) - eps) / float(nsub)
                    for isub in range(nsub): 
                        strain_inc(deps,dt)
                elif ptype == "stress":
                    dsig = (sig_from_voigt(sigtest) - sig) / float(nsub)
                    for isub in range(nsub): 
                        stress_inc(dsig,dt)                    
                recordt(epstest, sigtest)
                record(eps, sig)
    def strain_test(args):
        Commands.test(args, "strain")
    def stress_test(args):
        Commands.test(args, "stress")
    def path(args, ptype): # hidden command used by stress_path and strain_path
        global mode, n_dim, t
        global eps, sig
        global start_inc
        histrec()
        testname = args[0]
        if testname[-4:] != ".csv": testname = testname + ".csv"
        print("Reading from", testname)
        test_file = open(testname,"r")
        test_text = test_file.readlines()
        test_file.close()
        nsub = int(args[1])
        start_inc = True
        for iprint in range(len(test_text)):
            test_split = test_text[iprint].split(r'[ ,;]')
            ttest = float(test_split[0])
            if mode == 0:
                val = float(test_split[1])
            else:                
                val = np.array([float(tex) for tex in test_split[1:n_dim+1]])
            if start_inc:
                if ptype == "strain":
                    eps = eps_from_voigt(val)
                elif ptype == "stress":
                    sig = sig_from_voigt(val)
                delrec()
                record(eps, sig)
                start_inc= False
            else:
                dt = (ttest - t) / float(nsub)
                if ptype == "strain":
                    deps = (eps_from_voigt(val) - eps) / float(nsub)
                    for isub in range(nsub): 
                        strain_inc(deps,dt)
                elif ptype == "stress":
                    dsig = (sig_from_voigt(val) - sig) / float(nsub)
                    for isub in range(nsub): 
                        stress_inc(dsig,dt)
                record(eps, sig)
    def strain_path(args):
        Commands.path(args, "strain")
    def stress_path(args):
        Commands.path(args, "stress")
        
# commands for plotting and printing
    def printrec(args=[]):
        oname = "hout_" + model
        if len(args) > 0: oname = args[0]
        results_print(oname)
    def csv(args=[]):
        oname = "hcsv_" + model
        if len(args) > 0: oname = args[0]
        results_csv(oname)
    def specialprint(args=[]):
        oname = "hout_" + model
        if len(args) > 0: oname = args[0]
        if hasattr(hm,"specialprint"): hm.specialprint(oname, eps, sig, alp, chi)
    def printstate(args=[]):
        print("eps =", eps)
        print("sig =", sig)
        print("alp =", alp)
        print("chi =", chi)
    def plot(args=[]):
        pname = "hplot_" + model
        if len(args) > 0: pname = args[0]
        results_plot(pname)
    def plotCS(args=[]):
        pname = "hplot_" + model
        if len(args) > 0: pname = args[0]
        results_plotCS(pname)
    def graph(args):
        #print(args)
        pname = "hplot_" + model
        xsize = 6.0
        ysize = 6.0
        axes = args[0:2]
        if len(args) > 2: pname = args[2]
        if len(args) > 3: 
            xsize = args[3]
            ysize = args[4]
        results_graph(pname, axes, xsize, ysize)
    def specialplot(args=[]):
        pname = "hplot_" + model
        if len(args) > 0: pname = args[0]
        if hasattr(hm,"specialplot"): hm.specialplot(pname, title, rec, eps, sig, alp, chi)
    def high(args=[]):
        global hstart
        print("Start highlighting plot")
        histrec()
        hstart = len(rec[curr_test])-1
    def unhigh(args=[]):
        global hend, high
        print("Stop highlighting plot")
        histrec()
        hend = len(rec[curr_test])
        high[curr_test].append([hstart,hend])
    def pause(args=[]):
        if len(args) > 0:
            message = " ".join(args)
            if len(message) > 0: print(message)
        text = input("Processing paused: hit ENTER to continue (x to exit)... ")
        if text == "x" or text == "X": sys.exit()

# commands for recording and running history
    def start_history(args=[]):
        global history, history_rec
        print("Start recording history")
        history = True
        history_rec = []
    def end_history(args=[]):
        global history
        print("Stop recording history")
        history = False
    def run_history(args=[]):
        global history
        history = False
        if len(args) > 0: runs = int(args[0])
        else: runs = 1
        for irun in range(runs):
            print("Running history, cycle",irun+1,"of",runs)
            process(history_rec)
#XXX      
    def returnrec(args=[]):
        global rec
        ultest = rec
        qprint('Your file now contains last record')
        return(ultest)
    def returnstate():
        global rec
        return(rec[0][-1])
    def returnitern():
        global rec
        itern = len(rec[0])
        return(itern)

def set_up_analytical():
    global hdfde, hdfda, hd2fdede, hd2fdeda, hd2fdade, hd2fdada
    global hdgds, hdgda, hd2gdsds, hd2gdsda, hd2gdads, hd2gdada
    global hdyde, hdyds, hdyda, hdydc
    global hdwdc

    def set_h(name):
        qprint("  "+name+"...")
        if hasattr(hm,name): 
            hname = getattr(hm,name)
            qprint("        ...found")
        else:
            hname = undef
            qprint("        ...undefined")
        return hname      
    print("Setting up analytical differentials (if available)")
    hdfde =    set_h("dfde")
    hdfda =    set_h("dfda")
    hd2fdede = set_h("d2fdede")
    hd2fdeda = set_h("d2fdeda")
    hd2fdade = set_h("d2fdade")
    hd2fdada = set_h("d2fdada")
    hdgds =    set_h("dgds")
    hdgda =    set_h("dgda")
    hd2gdsds = set_h("d2gdsds")
    hd2gdsda = set_h("d2gdsda")
    hd2gdads = set_h("d2gdads")
    hd2gdada = set_h("d2gdada")
    hdyde =  set_h("dyde")
    hdyds =  set_h("dyds")
    hdyda =  set_h("dyda")
    hdydc =  set_h("dydc")
    hdwdc =  set_h("dwdc")
    
def set_up_auto():
    global adfde, adfda, ad2fdede, ad2fdeda, ad2fdade, ad2fdada
    global adgds, adgda, ad2gdsds, ad2gdsda, ad2gdads, ad2gdada
    global adyde, adyds, adyda, adydc
    global adwdc
    
    adfde =    undef
    adfda =    undef
    ad2fdede = undef
    ad2fdeda = undef
    ad2fdade = undef
    ad2fdada = undef
    adgds =    undef
    adgda =    undef
    ad2gdsds = undef
    ad2gdsda = undef
    ad2gdads = undef
    ad2gdada = undef
    adyde =  undef
    adyds =  undef
    adyda =  undef
    adydc =  undef
    adwdc =  undef
    if not auto:
        print("\nAutomatic differentials not available")
        return
    print("\nSetting up automatic differentials")
    if hasattr(hm,"f"):
        if hasattr(hm,"f_exclude"):
            qprint("f excluded from auto-differentials in", hm.file)
        else:
            qprint("Setting up auto-differentials of f")
            qprint("  dfde...")
            adfde = ag.jacobian(hm.f,0)
            qprint("  dfda...")
            adfda = ag.jacobian(hm.f,1)
            qprint("  d2fdede...")
            ad2fdede = ag.jacobian(adfde,0)
            qprint("  d2fdeda...")
            ad2fdeda = ag.jacobian(adfde,1)
            qprint("  d2fdade...")
            ad2fdade = ag.jacobian(adfda,0)
            qprint("  d2fdada...")
            ad2fdada = ag.jacobian(adfda,1)
    else:
        qprint("f not specified in", hm.file)
    if hasattr(hm,"g"):
        if hasattr(hm,"g_exclude"):
            qprint("g excluded from auto-differentials in", hm.file)
        else:
            qprint("Setting up auto-differentials of g")
            qprint("  dgds...")
            adgds = ag.jacobian(hm.g,0)
            qprint("  dgda...")
            adgda = ag.jacobian(hm.g,1)
            qprint("  d2gdsds...")
            ad2gdsds = ag.jacobian(adgds,0)
            qprint("  d2gdsda...")
            ad2gdsda = ag.jacobian(adgds,1)
            qprint("  d2gdads...")
            ad2gdads = ag.jacobian(adgda,0)
            qprint("  d2gdada...")
            ad2gdada = ag.jacobian(adgda,1)
    else:
        qprint("g not specified in", hm.file)   
    if hasattr(hm,"y"):
        if hasattr(hm,"y_exclude"):
            qprint("y excluded from auto-differentials in", hm.file)
        else:
            qprint("Setting up auto-differentials of y")
            qprint("  dyde...")
            adyde = ag.jacobian(hm.y,0)
            qprint("  dyds...")
            adyds = ag.jacobian(hm.y,1)
            qprint("  dyda...")
            adyda = ag.jacobian(hm.y,2)
            qprint("  dydc...")
            adydc = ag.jacobian(hm.y,3)
    else:
        qprint("y not specified in", hm.file)    
    if hasattr(hm,"w"):
        if hasattr(hm,"w_exclude"):
            qprint("w excluded from auto-differentials in", hm.file)
        else:
            qprint("Setting up auto-differential of w")
            qprint("  dwdc...")
            adwdc = ag.jacobian(hm.w,3)
    else:
        qprint("w not specified in", hm.file)    

def set_up_num():
    global ndfde, ndfda, nd2fdede, nd2fdeda, nd2fdade, nd2fdada
    global ndgds, ndgda, nd2gdsds, nd2gdsda, nd2gdads, nd2gdada
    global ndyde, ndyds, ndyda, ndydc
    global ndwdc

    ndfde =    undef
    ndfda =    undef
    nd2fdede = undef
    nd2fdeda = undef
    nd2fdade = undef
    nd2fdada = undef
    ndgds =    undef
    ndgda =    undef
    nd2gdsds = undef
    nd2gdsda = undef
    nd2gdads = undef
    nd2gdada = undef
    ndyde =  undef
    ndyds =  undef
    ndyda =  undef
    ndydc =  undef
    ndwdc =  undef
    print("\nSetting up numerical differentials")
    if hasattr(hm,"f"):
        qprint("Setting up numerical differentials of f")
        qprint("  dfde...")
        def ndfde(eps,alp): return numdiff_1(mode, n_dim, hm.f, eps, alp, epsi)
        qprint("  dfda...")
        def ndfda(eps,alp): return numdiff_2(mode, n_dim, n_int, hm.f, eps, alp, alpi)
        qprint("  d2fdede...")
        def nd2fdede(eps,alp): return numdiff2_1(mode, n_dim, hm.f, eps, alp, epsi)
        qprint("  d2fdeda...")
        def nd2fdeda(eps,alp): return numdiff2_2(mode, n_dim, n_int, hm.f, eps, alp, epsi, alpi)
        qprint("  d2fdade...")
        def nd2fdade(eps,alp): return numdiff2_3(mode, n_dim, n_int, hm.f, eps, alp, epsi, alpi)
        qprint("  d2fdada...")
        def nd2fdada(eps,alp): return numdiff2_4(mode, n_dim, n_int, hm.f, eps, alp, alpi)
    else:
        qprint("f not specified in", hm.file)
    if hasattr(hm,"g"):
        qprint("Setting up numerical differentials of g")
        qprint("  dgds...")
        def ndgds(sig,alp): return numdiff_1(mode, n_dim, hm.g, sig, alp, sigi)
        qprint("  dgda...")
        def ndgda(sig,alp): return numdiff_2(mode, n_dim, n_int, hm.g, sig, alp, alpi)
        qprint("  d2gdsds...")
        def nd2gdsds(sig,alp): return numdiff2_1(mode, n_dim, hm.g, sig, alp, sigi)
        qprint("  d2gdsda...")
        def nd2gdsda(sig,alp): return numdiff2_2(mode, n_dim, n_int, hm.g, sig, alp, sigi, alpi)
        qprint("  d2gdads...")
        def nd2gdads(sig,alp): return numdiff2_3(mode, n_dim, n_int, hm.g, sig, alp, sigi, alpi)
        qprint("  d2gdada...")
        def nd2gdada(sig,alp): return numdiff2_4(mode, n_dim, n_int, hm.g, sig, alp, alpi)
    else:
        qprint("g not specified in", hm.file)
    if hasattr(hm,"y"):
        qprint("Setting up numerical differentials of y")
        qprint("  dydc...")
        def ndydc(eps,sig,alp,chi): return numdiff_3(mode, n_dim, n_int, n_y, hm.y, eps,sig,alp,chi, chii)
        qprint("  dyde...")
        def ndyde(eps,sig,alp,chi): return numdiff_4e(mode, n_dim, n_int, n_y, hm.y, eps,sig,alp,chi, epsi)
        qprint("  dyds...")
        def ndyds(eps,sig,alp,chi): return numdiff_4s(mode, n_dim, n_int, n_y, hm.y, eps,sig,alp,chi, sigi)
        qprint("  dyda...")
        def ndyda(eps,sig,alp,chi): return numdiff_5(mode, n_dim, n_int, n_y, hm.y, eps,sig,alp,chi, alpi)
    else:
        qprint("y not specified in", hm.file)
    if hasattr(hm,"w"):
        qprint("Setting up numerical differential of w")
        qprint("  dwdc...")
        def ndwdc(eps,sig,alp,chi): return numdiff_6(mode, n_dim, n_int, hm.w, eps,sig,alp,chi, chii)
    else:
        qprint("w not specified in", hm.file)

def choose_diffs():
    global dfde, dfda, d2fdede, d2fdeda, d2fdade, d2fdada
    global dgds, dgda, d2gdsds, d2gdsda, d2gdads, d2gdada
    global dyde, dyds, dyda, dydc
    global dwdc
    
    def choose(name, hd, nd, ad):
        d = udef
        for i in range(3):
            if pref[i] == "analytical" and hd != udef: d = hd
            if pref[i] == "automatic"  and ad != udef: d = ad 
            if pref[i] == "numerical"  and nd != udef: d = nd
            if d != undef:
                qprint(name+":", pref[i])
                return d
        d = undef
        qprint(name+": undefined - will not run if this is required")
        return d        
    if not quiet:
        print("\nChoosing preferred differential methods")
    dfde    = choose("dfde",    hdfde,    ndfde,    adfde)
    dfda    = choose("dfda",    hdfda,    ndfda,    adfda)
    d2fdede = choose("d2fdede", hd2fdede, nd2fdede, ad2fdede)
    d2fdeda = choose("d2fdeda", hd2fdeda, nd2fdeda, ad2fdeda)
    d2fdade = choose("d2fdade", hd2fdade, nd2fdade, ad2fdade)
    d2fdada = choose("d2fdada", hd2fdada, nd2fdada, ad2fdada)
    dgds    = choose("dgds",    hdgds,    ndgds,    adgds)
    dgda    = choose("dgda",    hdgda,    ndgda,    adgda)
    d2gdsds = choose("d2gdsds", hd2gdsds, nd2gdsds, ad2gdsds)
    d2gdsda = choose("d2gdsda", hd2gdsda, nd2gdsda, ad2gdsda)
    d2gdads = choose("d2gdads", hd2gdads, nd2gdads, ad2gdads)
    d2gdada = choose("d2gdada", hd2gdada, nd2gdada, ad2gdada)
    dyde    = choose("dyde",    hdyde,    ndyde,    adyde)
    dyds    = choose("dyds",    hdyds,    ndyds,    adyds)
    dyda    = choose("dyda",    hdyda,    ndyda,    adyda)
    dydc    = choose("dydc",    hdydc,    ndydc,    adydc)
    dwdc    = choose("dwdc",    hdwdc,    ndwdc,    adwdc)
    print("")

def sig_f(e, a): 
#    print("eps in sig_f",e)
    return dfde(e, a)
def chi_f(e, a): # updated for weights
    return -np.einsum(ein_chi, rwt, dfda(e, a))
def eps_g(s, a): 
    return -dgds(s, a)
def chi_g(s, a): # updated for weights
    return -np.einsum(ein_chi, rwt, dgda(s, a))

def testch(text1, val1, text2, val2, failtext):
    global npass, nfail, fails
    print(text1, val1)
    print(text2, val2)
    if hasattr(val1, "shape") and hasattr(val2, "shape"):
        if val1.shape != val2.shape:
            print("\x1b[0;31mArrays different dimensions:",val1.shape,"and",val2.shape,"\x1b[0m")
            nfail += 1
            fails.append(failtext)
            pause()
            return
    maxv = np.maximum(np.max(val1),np.max(val2))
    minv = np.minimum(np.min(val1),np.min(val2))
    mv = np.maximum(maxv,-minv)
    #print("mv =",mv)
    testmat = np.isclose(val1, val2, rtol=0.0001, atol=0.000001*mv)
    if all(testmat.reshape(-1)): 
        print("\x1b[1;32m***PASSED***\n\x1b[0m")
        npass += 1
    else:
        print("\x1b[1;31m",testmat,"\x1b[0m")
        print("\x1b[1;31m***FAILED***\n\x1b[0m")
        nfail += 1
        fails.append(failtext)
        pause()

def testch2(ch, v1, v2, failtext):
    global nfail, npass, nmiss, fails, nzero, zeros
    try:
        _ = (e for e in v1)
        v1_iterable = True
    except TypeError:
        v1_iterable = False
    try:
        _ = (e for e in v2)
        v2_iterable = True
    except TypeError:
        v2_iterable = False
    if not v1_iterable:
        if v1 == "undefined":
            print(ch,"first variable undefined, comparison not possible")
            nmiss += 1
            return
    if not v2_iterable:
        if v2 == "undefined":
            print(ch,"second variable undefined, comparison not possible")
            nmiss += 1
            return
    if v1_iterable != v2_iterable:
        print(ch,"variables of different types, comparison not possible")
        nfail += 1
        fails.append(failtext)
        pause()
        return
    if v1_iterable and v2_iterable:
        if hasattr(v1,"shape") and hasattr(v2,"shape"):
            if v1.shape != v2.shape:
                print("\x1b[0;31mArrays different dimensions:",v1.shape,"and",v2.shape,"\x1b[0m")
                nfail += 1
                fails.append(failtext)
                pause()
                return
    maxv = np.maximum(np.max(v1),np.max(v2))
    minv = np.minimum(np.min(v1),np.min(v2))
    mv = np.maximum(maxv,-minv)
#    print("mv =",mv)
    testmat = np.isclose(v1, v2, rtol=0.0001, atol=0.000001*mv)
    if all(testmat.reshape(-1)): 
        print(ch,"\x1b[1;32m***PASSED***\x1b[0m")
        npass += 1
        if np.max(v1) == 0.0 and np.min(v1) == 0.0:
            nzero +=1
            zeros.append(failtext)
    else:
        print(ch,"\x1b[1;31m***FAILED***\x1b[0m")
        print("test\x1b[1;31m",testmat,"\x1b[0m")
        pprint(v1,"va11",14,6)
        pprint(v2,"val2",14,6)
        nfail += 1
        fails.append(failtext)
        pause()
    
def testch1(text, hval, aval, nval):
    pprint(hval,"Analytical "+text, 14, 6)
    pprint(aval,"Automatic  "+text, 14, 6)
    pprint(nval,"Numerical  "+text, 14, 6)
    def ex(fun):
        if hasattr(fun,"shape"):
            exists = True
        else:
            exists = (fun != "undefined")
        return exists
    if ex(hval) and ex(aval):
        testch2("analytical v. automatic: ", hval, aval, "analytical v. automatic:  "+text)
    if ex(aval) and ex(nval):
        testch2("automatic  v. numerical: ", aval, nval, "automatic  v. numerical:  "+text)
    if ex(nval) and ex(hval):
        testch2("numerical  v. analytical:", nval, hval, "numerical  v. analytical: "+text)

def check(arg="unknown"):
    print("")
    print("+---------------------------------------------------------+")
    print("| HyperCheck: checking routine for hyperplasticity models |")
    print("| (c) G.T. Houlsby, 2019-2021                             |")
    print("+---------------------------------------------------------+")
    global epsi, sigi, alpi, chii
    global n_int, n_dim, n_y, mode, hm
    global npass, nfail, fails, nzero, zeros, nmiss
    global hm
    global pref
    
    print("Current directory: " + os.getcwd())
    if os.path.isfile("HyperDrive.cfg"):
        input_file = open("HyperDrive.cfg", 'r')
        last = input_file.readline().split(",")
        drive_file_name = last[0]
        check_file_name = last[1]
        input_file.close()
    else:
        drive_file_name = "hyper.dat"
        check_file_name = "unknown"
    if arg != "unknown": 
        check_file_name = arg
    input_found = False
    count = 0
    while not input_found:
        count += 1
        if count > 5: error("Too many tries")
        check_file_temp = read_def("Input model for checking", check_file_name)
        if os.path.isfile(check_file_temp + ".py"):
            check_file_name = check_file_temp
            print("Importing hyperplasticity model:", check_file_name)
            hm = importlib.import_module(check_file_name)
            if hasattr(hm, "deriv"): hm.deriv()
            #hm.setvals()
            print("Description:", hm.name)
            input_found = True
        else:
            print("File not found:", check_file_temp)
    last_file = open("HyperDrive.cfg", 'w')
    last_file.writelines([drive_file_name+","+check_file_name])
    last_file.close()

    npass = 0
    nfail = 0
    fails = []
    nzero = 0
    zeros = []
    nmiss = 0
    model = hm.file
    mode  = hm.mode
    if mode == 0:
        n_dim = 1
    else:
        n_dim = hm.ndim
    n_int = hm.n_int
    n_y   = hm.n_y
    
    epsi = 0.000001
    sigi = 0.0001
    alpi = 0.0000011
    chii = 0.00011
    if hasattr(hm, "epsi"): epsi = hm.epsi
    if hasattr(hm, "sigi"): sigi = hm.sigi
    if hasattr(hm, "alpi"): alpi = hm.alpi
    if hasattr(hm, "chii"): chii = hm.chii
    
    pref = ["analytical", "automatic", "numerical"]
    set_up_analytical()
    set_up_auto()
    set_up_num()
    choose_diffs()
    
    ein_inner, ein_contract, sig, eps, chi, alp = modestart_short()
    
    if hasattr(hm, "check_eps"): 
        eps = hm.check_eps
    else:
        text = input("Input test strain: ",)
        if mode == 0: 
            if len(text) == 0.0: text = "0.01"
            eps = float(text)
        elif mode == 1: 
            eps = np.array([float(item) for item in re.split(r'[ ,;]',text)])
    pprint(eps, "eps =", 12, 6)
    
    if hasattr(hm, "check_sig"): 
        sig = hm.check_sig
    else:
        text = input("Input test stress: ",)
        if mode == 0: 
            if len(text) == 0.0: text = "1.0"
            sig = float(text)
        elif mode == 1: 
            sig = np.array([float(item) for item in re.split(r'[ ,;]',text)])
    pprint(sig, "sig =", 12, 4)
    
    if hasattr(hm, "check_alp"): 
        alp = hm.check_alp
    pprint(alp, "alp =", 12, 6)
    
    if hasattr(hm, "check_chi"): 
        chi = hm.check_chi
    pprint(chi, "chi =", 12, 4)
    
    princon(hm)   
    input("Hit ENTER to start checks",)
    
    if hasattr(hm, "f") and hasattr(hm, "g"):
        print("Checking consistency of f and g formulations ...\n")
        print("Checking inverse relations...")
        sigt = sig_f(eps, alp)
        epst = eps_g(sig, alp)
        print("Stress from strain:", sigt)
        print("Strain from stress:", epst)
        testch("Strain:             ", eps, "-> stress -> strain:", eps_g(sigt, alp),"strain v stress -> strain")
        testch("Stress:             ", sig, "-> strain -> stress:", sig_f(epst, alp),"stress v strain -> stress")
    
        print("Checking chi from different routes...")
        testch("from eps and f:", chi_f(eps, alp), "from g:        ", chi_g(sigt, alp),"chi from eps")
        testch("from sig and g:", chi_g(sig, alp), "from f:        ", chi_f(epst, alp),"chi from sig")
    
        print("Checking Legendre transform...")
        W = np.einsum(ein_contract, sigt, eps)
        testch("from eps: f + (-g) =", hm.f(eps,alp) - hm.g(sigt,alp), "sig.eps =           ", W,"f + (-g): 1")
        W = np.einsum(ein_contract, sig, epst)
        testch("from sig: f + (-g) =", hm.f(epst,alp) - hm.g(sig,alp), "sig.eps =           ", W,"f + (-g): 2")
    
        print("Checking elastic stiffness and compliance at specified strain...")
        if mode == 0:
            unit = 1.0
        elif mode == 1:
            unit = np.eye(n_dim)
        else:
            unit = Utils.IIsym
        D = d2fdede(eps,alp)
        C = -d2gdsds(sigt,alp)
        DC = np.einsum(ein_inner,D,C)
        print("Stiffness matrix D  =\n",D)
        print("Compliance matrix C =\n",C)
        testch("Product DC = \n",DC,"unit matrix =\n", unit,"DC")
        print("and at specified stress...")
        C = -d2gdsds(sig,alp)
        D = d2fdede(epst,alp)
        CD = np.einsum(ein_inner,C,D)
        print("Compliance matrix D =\n",C)
        print("Stiffness matrix C  =\n",D)
        testch("Product CD = \n",CD,"unit matrix =\n", unit,"CD")
    
    if hasattr(hm, "f"):
        print("Checking derivatives of f...")
#        pprint(eps,"eps")
#        pprint(alp,"alp")
#        pprint(hdfde(eps,alp),"hdfde")
#        pprint(ndfde(eps,alp),"ndfde")
#        pprint(adfde(eps,alp),"adfde")
        testch1("dfde", hdfde(eps,alp), adfde(eps,alp), ndfde(eps,alp))
        testch1("dfda", hdfda(eps,alp), adfda(eps,alp), ndfda(eps,alp))
        testch1("d2fdede", hd2fdede(eps,alp), ad2fdede(eps,alp), nd2fdede(eps,alp))
        testch1("d2fdeda", hd2fdeda(eps,alp), ad2fdeda(eps,alp), nd2fdeda(eps,alp))
        testch1("d2fdade", hd2fdade(eps,alp), ad2fdade(eps,alp), nd2fdade(eps,alp))
        testch1("d2fdada", hd2fdada(eps,alp), ad2fdada(eps,alp), nd2fdada(eps,alp))
    if hasattr(hm, "g"):
        print("Checking derivatives of g...")
        testch1("dgds", hdgds(sig,alp), adgds(sig,alp), ndgds(sig,alp))
        testch1("dgda", hdgda(sig,alp), adgda(sig,alp), ndgda(sig,alp))
        testch1("d2gdsds", hd2gdsds(sig,alp), ad2gdsds(sig,alp), nd2gdsds(sig,alp))
        testch1("d2gdsda", hd2gdsda(sig,alp), ad2gdsda(sig,alp), nd2gdsda(sig,alp))
        testch1("d2gdads", hd2gdads(sig,alp), ad2gdads(sig,alp), nd2gdads(sig,alp))
        testch1("d2gdada", hd2gdada(sig,alp), ad2gdada(sig,alp), nd2gdada(sig,alp))
        
    if hasattr(hm, "y"):
        print("Checking derivatives of y...")
        testch1("dyde", hdyde(eps,sig,alp,chi), adyde(eps,sig,alp,chi), ndyde(eps,sig,alp,chi))
        testch1("dyds", hdyds(eps,sig,alp,chi), adyds(eps,sig,alp,chi), ndyds(eps,sig,alp,chi))
        testch1("dyda", hdyda(eps,sig,alp,chi), adyda(eps,sig,alp,chi), ndyda(eps,sig,alp,chi))
        testch1("dydc", hdydc(eps,sig,alp,chi), adydc(eps,sig,alp,chi), ndydc(eps,sig,alp,chi))
                
    if hasattr(hm, "w"):
        print("Checking derivative of w...")
        testch1("dwdc", hdwdc(eps,sig,alp,chi), adwdc(eps,sig,alp,chi), ndwdc(eps,sig,alp,chi))
        
    print("Checks complete for:",model+",",npass,"passed,",nfail,"failed,",nmiss,"missed checks")
    if not hasattr(hm, "f"): print("hm.f not present")
    if not hasattr(hm, "g"): print("hm.g not present")
    if not hasattr(hm, "y"): print("hm.y not present")
    if not hasattr(hm, "w"): print("hm.w not present")
    if nfail > 0:
        print("Summary of fails:")
        for text in fails: print("  ",text)
    if nzero > 0:
        print("Warning - zero values (may not have been rigorously checked):")
        for text in zeros: print("  ",text)
            
def princon(hm):
    if quiet:
        return()
    print("Constants for model:",hm.file)
    print(hm.const)
    print("Derived values:")
    for i in range(len(hm.const)):
        print(hm.name_const[i] + " =", hm.const[i])

def calcsum(vec):
    global eps, sig, chi, alp, nstep, nsub, deps, sigrec, const, hm, n_int, var
    sumerr = 0.0
    print(vec)
    hm.const = copy.deepcopy(const)
    for i in range(n_int):
       hm.const[2+2*i] = vec[i]
    if hasattr(hm, "deriv"): hm.deriv()
    eps = 0.0
    sig = 0.0
    alp = np.zeros(n_int)
    chi = np.zeros(n_int)
    if "_h" in var:
        alp = np.zeros(n_int+1)
        chi = np.zeros(n_int+1)
    if "_cbh" in var:
        alp = np.zeros(n_int+1)
        chi = np.zeros(n_int+1)
    for step in range(nstep):
        for i in range(nsub):
            strain_inc_f_spec(deps)
        error = sig - sigrec[step]
        sumerr += error**2
    return sumerr

def solve_L(yo, Lmatp, Lrhsp):
    global Lmat, Lrhs, L
    ytol = -0.00001                          # tolerance on triggering yield
    Lmat = np.eye(hm.n_y)                    # initialise matrix and RHS for elastic solution
    Lrhs = np.zeros(hm.n_y)
    for N in range(hm.n_y):                  # loop over yield surfaces
        if yo[N] > ytol:                     # if this surface yielding ...
            Lmat[N] = Lmatp[N]               # over-write line in matrix and RHS with plastic solution
            Lrhs[N] = Lrhsp[N]
 #           print("Lmat =",Lmat)
 #           print("Lrhs =",Lrhs)
    L = np.linalg.solve(Lmat, Lrhs)          # solve for plastic multipliers
    L = np.array([max(Lv, 0.0) for Lv in L]) # make plastic multipliers non-negative
#    L = np.array([abs(Lv) for Lv in L])       #make plastic multipliers positive (for test only)
    return L

def optim(base, variant):
    global sigrec, hm, nstep, nsub, n_int, const, deps, var
    global eps, sig, alp, chi   
    var = variant
    sigrec = np.zeros(nstep)
    print("calculate base curve from", base)
    hm = importlib.import_module(base)    
    hm.setvals()
    hm.const = copy.deepcopy(const[:2+2*n_int])
    if hasattr(hm, "deriv"): hm.deriv()
    eps = 0.0
    sig = 0.0
    alp = np.zeros(n_int)
    chi = np.zeros(n_int)
    for step in range(nstep):
       for i in range(nsub): strain_inc_f_spec(deps)
       sigrec[step] = sig       
    print("optimise", variant)
    hm = importlib.import_module(variant)
    hm.setvals()
    vec = np.zeros(n_int)
    for i in range(n_int): vec[i] = const[2*i+2]
    print(vec,calcsum(vec))
    #optimize.Bounds(0.0,np.inf)
    bnds = optimize.Bounds(0.0001,np.inf)
    resultop = optimize.minimize(calcsum, vec, method='L-BFGS-B', bounds=bnds)
    vec = resultop.x
    print(vec,calcsum(vec))
    for i in range(n_int): const[2+2*i] = resultop.x[i]
    return const

def derive_from_points(modeltype, epsin, sigin, Einf=0.0, epsmax=0.0, HARM_R=0.0):
    global nstep, nsub, n_int, const, deps
    eps = np.array(epsin)
    sig = np.array(sigin)
    le = len(eps)
    ls = len(sig)
    if ls != le:
        print("Unequal numbers of values")
        le = min(le, ls)
    n_int = le
    E = np.zeros(n_int+1)
    E[0] = sig[0] / eps[0]
    for i in range(1,n_int):
        E[i] = (sig[i]-sig[i-1]) / (eps[i]-eps[i-1])
    E[n_int] = Einf
    print("eps =",eps)
    print("sig =",sig)
    print("E   =",E)
    k = np.zeros(n_int)
    H = np.zeros(n_int)
    if "ser" in modeltype:
        print("Series parameters")
        E0 = E[0]
        for i in range(n_int):
            k[i] = sig[i]
            H[i] = E[i+1]*E[i]/(E[i] - E[i+1])
        const = [E0, n_int]
        for i in range(n_int):
            const.append(round(k[i],6))
            const.append(round(H[i],6))
        base = "h1epmk_ser"
    elif "par" in modeltype:
        print("Parallel parameters")
        for i in range(n_int):
            H[i] = E[i] - E[i+1]                
            k[i] = eps[i]*H[i]
        const = [Einf, n_int]
        for i in range(n_int):
            const.append(round(k[i],6))
            const.append(round(H[i],6))
        base = "h1epmk_par"
    elif "nest" in modeltype:
        print("Nested parameters")
        E0 = E[0]
        for i in range(n_int):
            k[i] = sig[i] - sig[i-1]
            H[i] = E[i+1]*E[i]/(E[i] - E[i+1])                
        k[0] = sig[0]
        const = [E0, n_int]
        for i in range(n_int):
            const.append(round(k[i],6))
            const.append(round(H[i],6))
        base = "h1epmk_nest"
    if "_b" in modeltype: #now optimise for bounding surface model
        print("Optimise parameters for _b option")
        nstep = 100
        nsub = 10
        if epsmax == 0.0:
            epsmax = 1.5*eps[n_int-1]
            print("setting epsmax =",epsmax)
        deps = epsmax / float(nstep*nsub)
        const = optim(base, modeltype)
        for i in range(2,2+2*n_int):
            const[i] = round(const[i],6)
    if "_h" in modeltype: #now optimise for HARM model
        print("Optimise parameters for _h option")
        nstep = 100
        nsub = 10
        if epsmax == 0.0:
            epsmax = 1.5*eps[n_int-1]
            print("setting epsmax =",epsmax)
        deps = epsmax / float(nstep*nsub)
        const.append(HARM_R)
        const = optim(base, modeltype)
        for i in range(2,2+2*n_int):
            const[i] = round(const[i],6)
    if "_cbh" in modeltype: #now optimise for bounding HARM model
        print("Optimise parameters for _cbh option")
        nstep = 100
        nsub = 10
        if epsmax == 0.0:
            epsmax = 1.5*eps[n_int-1]
            print("setting epsmax =",epsmax)
        deps = epsmax / float(nstep*nsub)
        const.append(HARM_R)
        const = optim(base, modeltype)
        for i in range(2,2+2*n_int):
            const[i] = round(const[i],6)
    return const

def numdiff_1(mode, n_dim, fun, var, alp, vari):
    if mode == 0: 
        var1 = var - vari  
        var2 = var + vari
        f1 = fun(var1,alp)
        f2 = fun(var2,alp)
        num = (f2 - f1) / (2.0*vari)
    elif mode == 1: 
        num = np.zeros([n_dim])
        for i in range(n_dim):
            var1 = copy.deepcopy(var)
            var2 = copy.deepcopy(var)
            var1[i] = var[i] - vari  
            var2[i] = var[i] + vari    
            f1 = fun(var1,alp)
            f2 = fun(var2,alp)
            num[i] = (f2 - f1) / (2.0*vari)
    else:
        num = np.zeros([n_dim,n_dim])
        for i in range(n_dim):
            for j in range(n_dim):
                var1 = copy.deepcopy(var)
                var2 = copy.deepcopy(var)
                var1[i,j] = var[i,j] - vari  
                var2[i,j] = var[i,j] + vari    
                f1 = fun(var1,alp)
                f2 = fun(var2,alp)
                num[i,j] = (f2 - f1) / (2.0*vari)
    return num
def numdiff_2(mode, n_dim, n_int, fun, var, alp, alpi):
    if mode == 0:
        num = np.zeros([n_int])
        for k in range(n_int):
            alp1 = copy.deepcopy(alp)
            alp2 = copy.deepcopy(alp)
            alp1[k] = alp[k] - alpi  
            alp2[k] = alp[k] + alpi
            f1 = fun(var,alp1)
            f2 = fun(var,alp2)
            num[k] = (f2 - f1) / (2.0*alpi)
    elif mode == 1: 
        num = np.zeros([n_int,n_dim])
        for k in range(n_int):
            for i in range(n_dim):
                alp1 = copy.deepcopy(alp)
                alp2 = copy.deepcopy(alp)
                alp1[k,i] = alp[k,i] - alpi  
                alp2[k,i] = alp[k,i] + alpi
                f1 = fun(var,alp1)
                f2 = fun(var,alp2)
                num[k,i] = (f2 - f1) / (2.0*alpi)
    else: 
        num = np.zeros([n_int,n_dim,n_dim])
        for k in range(n_int):
            for i in range(n_dim):
                for j in range(n_dim):
                    alp1 = copy.deepcopy(alp)
                    alp2 = copy.deepcopy(alp)
                    alp1[k,i,j] = alp[k,i,j] - alpi  
                    alp2[k,i,j] = alp[k,i,j] + alpi
                    f1 = fun(var,alp1)
                    f2 = fun(var,alp2)
                    num[k,i,j] = (f2 - f1) / (2.0*alpi)
    return num
def numdiff_3(mode, n_dim, n_int, n_y, fun, eps,sig,alp,chi, chii):
    if mode == 0:
        num = np.zeros([n_y,n_int])
        for k in range(n_int):
            chi1 = copy.deepcopy(chi)
            chi2 = copy.deepcopy(chi)
            chi1[k] = chi[k] - chii  
            chi2[k] = chi[k] + chii
            f1 = fun(eps,sig,alp,chi1)
            f2 = fun(eps,sig,alp,chi2)
            for l in range(n_y):
                num[l,k] = (f2[l] - f1[l]) / (2.0*chii)
    elif mode == 1:
        num = np.zeros([n_y,n_int,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                chi1 = copy.deepcopy(chi)
                chi2 = copy.deepcopy(chi)
                chi1[k,i] = chi[k,i] - chii  
                chi2[k,i] = chi[k,i] + chii
                f1 = fun(eps,sig,alp,chi1)
                f2 = fun(eps,sig,alp,chi2)
                for l in range(n_y):
                    num[l,k,i] = (f2[l] - f1[l]) / (2.0*chii)
    else:
        num = np.zeros([n_y,n_int,n_dim,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                for j in range(n_dim):   
                    chi1 = copy.deepcopy(chi)
                    chi2 = copy.deepcopy(chi)
                    chi1[k,i,j] = chi[k,i,j] - chii  
                    chi2[k,i,j] = chi[k,i,j] + chii
                    f1 = fun(eps,sig,alp,chi1)
                    f2 = fun(eps,sig,alp,chi2)
                    for l in range(n_y):
                        num[l,k,i,j] = (f2[l] - f1[l]) / (2.0*chii)
    return num
def numdiff_4e(mode, n_dim, n_int, n_y, fun, eps,sig,alp,chi, vari):
    if mode == 0: 
        num = np.zeros([n_y])
        var1 = eps - vari  
        var2 = eps + vari
        f1 = fun(var1,sig,alp,chi)
        f2 = fun(var2,sig,alp,chi)
        for l in range(n_y):
            num[l] = (f2[l] - f1[l]) / (2.0*vari)
    elif mode == 1: 
        num = np.zeros([n_y,n_dim])
        for i in range(n_dim):   
            var1 = copy.deepcopy(eps)
            var2 = copy.deepcopy(eps)
            var1[i] = eps[i] - vari  
            var2[i] = eps[i] + vari
            f1 = fun(var1,sig,alp,chi)
            f2 = fun(var2,sig,alp,chi)
            for l in range(n_y):
                num[l,i] = (f2[l] - f1[l]) / (2.0*vari)
    else: 
        num = np.zeros([n_y,n_dim,n_dim])
        for i in range(n_dim):   
            for j in range(n_dim):   
                var1 = copy.deepcopy(eps)
                var2 = copy.deepcopy(eps)
                var1[i,j] = eps[i,j] - vari  
                var2[i,j] = eps[i,j] + vari
                f1 = fun(var1,sig,alp,chi)
                f2 = fun(var2,sig,alp,chi)
                for l in range(n_y):
                    num[l,i,j] = (f2[l] - f1[l]) / (2.0*vari)
    return num         
def numdiff_4s(mode, n_dim, n_int, n_y, fun, eps,sig,alp,chi, vari):
    if mode == 0: 
        num = np.zeros([n_y])
        var1 = sig - vari  
        var2 = sig + vari
        f1 = fun(eps,var1,alp,chi)
        f2 = fun(eps,var2,alp,chi)
        for l in range(n_y):
            num[l] = (f2[l] - f1[l]) / (2.0*vari)
    elif mode == 1: 
        num = np.zeros([n_y,n_dim])
        for i in range(n_dim):   
            var1 = copy.deepcopy(sig)
            var2 = copy.deepcopy(sig)
            var1[i] = sig[i] - vari  
            var2[i] = sig[i] + vari
            f1 = fun(eps,var1,alp,chi)
            f2 = fun(eps,var2,alp,chi)
            for l in range(n_y):
                num[l,i] = (f2[l] - f1[l]) / (2.0*vari)
    else: 
        num = np.zeros([n_y,n_dim,n_dim])
        for i in range(n_dim):   
            for j in range(n_dim):   
                var1 = copy.deepcopy(sig)
                var2 = copy.deepcopy(sig)
                var1[i,j] = sig[i,j] - vari  
                var2[i,j] = sig[i,j] + vari
                f1 = fun(eps,var1,alp,chi)
                f2 = fun(eps,var2,alp,chi)
                for l in range(n_y):
                    num[l,i,j] = (f2[l] - f1[l]) / (2.0*vari)
    return num         
def numdiff_5(mode, n_dim, n_int, n_y, fun, eps,sig,alp,chi, alpi):
    if mode == 0: 
        num = np.zeros([n_y,n_int])
        for k in range(n_int):
            alp1 = copy.deepcopy(alp)
            alp2 = copy.deepcopy(alp)
            alp1[k] = alp[k] - alpi  
            alp2[k] = alp[k] + alpi
            f1 = fun(eps,sig,alp1,chi)
            f2 = fun(eps,sig,alp2,chi)
            for l in range(n_y):
                num[l,k] = (f2[l] - f1[l]) / (2.0*alpi)
    elif mode == 1: 
        num = np.zeros([n_y,n_int,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                alp1 = copy.deepcopy(alp)
                alp2 = copy.deepcopy(alp)
                alp1[k,i] = alp[k,i] - alpi  
                alp2[k,i] = alp[k,i] + alpi
                f1 = fun(eps,sig,alp1,chi)
                f2 = fun(eps,sig,alp2,chi)
                for l in range(n_y):
                    num[l,k,i] = (f2[l] - f1[l]) / (2.0*alpi)
    else: 
        num = np.zeros([n_y,n_int,n_dim,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                for j in range(n_dim):   
                    alp1 = copy.deepcopy(alp)
                    alp2 = copy.deepcopy(alp)
                    alp1[k,i,j] = alp[k,i,j] - alpi  
                    alp2[k,i,j] = alp[k,i,j] + alpi
                    f1 = fun(eps,sig,alp1,chi)
                    f2 = fun(eps,sig,alp2,chi)
                    for l in range(n_y):
                        num[l,k,i,j] = (f2[l] - f1[l]) / (2.0*alpi)
    return num                
def numdiff_6(mode, n_dim, n_int, fun, eps,sig,alp,chi, chii):
    if mode == 0: 
        num = np.zeros([n_int])
        for k in range(n_int):
            chi1 = copy.deepcopy(chi)
            chi2 = copy.deepcopy(chi)
            chi1[k] = chi[k] - chii  
            chi2[k] = chi[k] + chii
            f1 = fun(eps,sig,alp,chi1)
            f2 = fun(eps,sig,alp,chi2)
            num[k] = (f2 - f1) / (2.0*chii)
    elif mode == 1: 
        num = np.zeros([n_int,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                chi1 = copy.deepcopy(chi)
                chi2 = copy.deepcopy(chi)
                chi1[k,i] = chi[k,i] - chii  
                chi2[k,i] = chi[k,i] + chii
                f1 = fun(eps,sig,alp,chi1)
                f2 = fun(eps,sig,alp,chi2)
                num[k,i] = (f2 - f1) / (2.0*chii)
    else: 
        num = np.zeros([n_int,n_dim,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                for j in range(n_dim):   
                    chi1 = copy.deepcopy(chi)
                    chi2 = copy.deepcopy(chi)
                    chi1[k,i,j] = chi[k,i,j] - chii  
                    chi2[k,i,j] = chi[k,i,j] + chii
                    f1 = fun(eps,sig,alp,chi1)
                    f2 = fun(eps,sig,alp,chi2)
                num[k,i,j] = (f2 - f1) / (2.0*chii)
    return num
def numdiff_6a(mode, n_dim, n_int, fun, eps,sig,alp,chi, chii):
    f0 = fun(eps,sig,alp,chi)
    if mode == 0: 
        num = np.zeros([n_int])
        for k in range(n_int):
            chi1 = copy.deepcopy(chi)
            chi2 = copy.deepcopy(chi)
            chi1[k] = chi[k] - chii  
            chi2[k] = chi[k] + chii
            f1 = fun(eps,sig,alp,chi1)
            f2 = fun(eps,sig,alp,chi2)
            if abs(f2-f0) > abs(f1-f0):
                num[k] = (f2 - f0) / chii
            else:
                num[k] = (f0 - f1) / chii
    elif mode == 1: 
        num = np.zeros([n_int,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                chi1 = copy.deepcopy(chi)
                chi2 = copy.deepcopy(chi)
                chi1[k,i] = chi[k,i] - chii  
                chi2[k,i] = chi[k,i] + chii
                f1 = fun(eps,sig,alp,chi1)
                f2 = fun(eps,sig,alp,chi2)
                if abs(f2-f0) > abs(f1-f0):
                    num[k,i] = (f2 - f0) / chii
                else:
                    num[k,i] = (f0 - f1) / chii
    else: 
        num = np.zeros([n_int,n_dim,n_dim])
        for k in range(n_int):
            for i in range(n_dim):   
                for j in range(n_dim):   
                    chi1 = copy.deepcopy(chi)
                    chi2 = copy.deepcopy(chi)
                    chi1[k,i,j] = chi[k,i,j] - chii  
                    chi2[k,i,j] = chi[k,i,j] + chii
                    f1 = fun(eps,sig,alp,chi1)
                    f2 = fun(eps,sig,alp,chi2)
                    if abs(f2-f0) > abs(f1-f0):
                        num[k,i,j] = (f2 - f0) / chii
                    else:
                        num[k,i,j] = (f0 - f1) / chii
    return num
def numdiff2_1(mode, n_dim, fun, var, alp, vari):
    if mode == 0: 
        num = 0.0
        var1 = var - vari  
        var3 = var + vari
        f1 = fun(var1,alp)
        f2 = fun(var,alp)
        f3 = fun(var3,alp)
        num = (f1 - 2.0*f2 + f3) / (vari**2)
    elif mode == 1: 
        num = np.zeros([n_dim,n_dim])
        for i in range(n_dim):
            for j in range(n_dim):
                if i==j:
                    var1 = copy.deepcopy(var)
                    var3 = copy.deepcopy(var)
                    var1[i] = var[i] - vari  
                    var3[i] = var[i] + vari
                    f1 = fun(var1,alp)
                    f2 = fun(var,alp)
                    f3 = fun(var3,alp)
                    num[i,i] = (f1 - 2.0*f2 + f3) / (vari**2)
                else:
                    var1 = copy.deepcopy(var)
                    var2 = copy.deepcopy(var)
                    var3 = copy.deepcopy(var)
                    var4 = copy.deepcopy(var)
                    var1[i] = var[i] - vari  
                    var1[j] = var[j] - vari  
                    var2[i] = var[i] - vari  
                    var2[j] = var[j] + vari  
                    var3[i] = var[i] + vari  
                    var3[j] = var[j] - vari  
                    var4[i] = var[i] + vari  
                    var4[j] = var[j] + vari  
                    f1 = fun(var1,alp)
                    f2 = fun(var2,alp)
                    f3 = fun(var3,alp)
                    f4 = fun(var4,alp)
                    num[i,j] = (f1 - f2 - f3 + f4) / (4.0*(vari**2))
    else: 
        num = np.zeros([n_dim,n_dim,n_dim,n_dim])
        for i in range(n_dim):
            for j in range(n_dim):
                for k in range(n_dim):
                    for l in range(n_dim):
                        if i==k and j==l:
                            var1 = copy.deepcopy(var)
                            var3 = copy.deepcopy(var)
                            var1[i,j] = var[i,j] - vari  
                            var3[i,j] = var[i,j] + vari
                            f1 = fun(var1,alp)
                            f2 = fun(var,alp)
                            f3 = fun(var3,alp)
                            num[i,j,i,j] = (f1 - 2.0*f2 + f3) / (vari**2)
                        else:
                            var1 = copy.deepcopy(var)
                            var2 = copy.deepcopy(var)
                            var3 = copy.deepcopy(var)
                            var4 = copy.deepcopy(var)
                            var1[i,j] = var[i,j] - vari  
                            var1[k,l] = var[k,l] - vari  
                            var2[i,j] = var[i,j] - vari  
                            var2[k,l] = var[k,l] + vari  
                            var3[i,j] = var[i,j] + vari  
                            var3[k,l] = var[k,l] - vari  
                            var4[i,j] = var[i,j] + vari  
                            var4[k,l] = var[k,l] + vari  
                            f1 = fun(var1,alp)
                            f2 = fun(var2,alp)
                            f3 = fun(var3,alp)
                            f4 = fun(var4,alp)
                            num[i,j,k,l] = (f1 - f2 - f3 + f4) / (4.0*(vari**2))
    return num
def numdiff2_2(mode, n_dim, n_int, fun, var, alp, vari, alpi):
    if mode == 0: 
        num = np.zeros(n_int)
        for N in range(n_int):
            alp1 = copy.deepcopy(alp)
            alp2 = copy.deepcopy(alp)
            var1 = var - vari  
            var2 = var + vari
            alp1[N] = alp[N] - alpi  
            alp2[N] = alp[N] + alpi
            f1 = fun(var1,alp1)
            f2 = fun(var2,alp1)
            f3 = fun(var1,alp2)
            f4 = fun(var2,alp2)
            num[N] = (f1 - f2 - f3 + f4) / (4.0*vari*alpi)
    elif mode == 1: 
        num = np.zeros([n_dim,n_int,n_dim])
        for N in range(n_int):
            for i in range(n_dim):
                for j in range(n_dim):
                    var1 = copy.deepcopy(var)
                    var2 = copy.deepcopy(var)
                    alp1 = copy.deepcopy(alp)
                    alp2 = copy.deepcopy(alp)
                    var1[i] = var[i] - vari  
                    var2[i] = var[i] + vari
                    alp1[N,j] = alp[N,j] - alpi  
                    alp2[N,j] = alp[N,j] + alpi
                    f1 = fun(var1,alp1)
                    f2 = fun(var2,alp1)
                    f3 = fun(var1,alp2)
                    f4 = fun(var2,alp2)
                    num[i,N,j] = (f1 - f2 - f3 + f4) / (4.0*vari*alpi)
    else: 
        num = np.zeros([n_dim,n_dim,n_int,n_dim,n_dim])
        for N in range(n_int):
            for i in range(n_dim):
                for j in range(n_dim):
                    for k in range(n_dim):
                        for l in range(n_dim):
                            var1 = copy.deepcopy(var)
                            var2 = copy.deepcopy(var)
                            alp1 = copy.deepcopy(alp)
                            alp2 = copy.deepcopy(alp)
                            var1[i,j] = var[i,j] - vari  
                            var2[i,j] = var[i,j] + vari
                            alp1[N,k,l] = alp[N,k,l] - alpi  
                            alp2[N,k,l] = alp[N,k,l] + alpi
                            f1 = fun(var1,alp1)
                            f2 = fun(var2,alp1)
                            f3 = fun(var1,alp2)
                            f4 = fun(var2,alp2)
                            num[i,j,N,k,l] = (f1 - f2 - f3 + f4) / (4.0*vari*alpi)
    return num
def numdiff2_3(mode, n_dim, n_int, fun, var, alp, vari, alpi):
    if mode == 0: 
        num = np.zeros(n_int)
        for N in range(n_int):
            alp1 = copy.deepcopy(alp)
            alp2 = copy.deepcopy(alp)
            var1 = var - vari  
            var2 = var + vari
            alp1[N] = alp[N] - alpi  
            alp2[N] = alp[N] + alpi
            f1 = fun(var1,alp1)
            f2 = fun(var2,alp1)
            f3 = fun(var1,alp2)
            f4 = fun(var2,alp2)
            num[N] = (f1 - f2 - f3 + f4) / (4.0*vari*alpi)
    elif mode == 1: 
        num = np.zeros([n_int,n_dim,n_dim])
        for N in range(n_int):
            for i in range(n_dim):
                for j in range(n_dim):
                    var1 = copy.deepcopy(var)
                    var2 = copy.deepcopy(var)
                    alp1 = copy.deepcopy(alp)
                    alp2 = copy.deepcopy(alp)
                    var1[i] = var[i] - vari  
                    var2[i] = var[i] + vari
                    alp1[N,j] = alp[N,j] - alpi  
                    alp2[N,j] = alp[N,j] + alpi
                    f1 = fun(var1,alp1)
                    f2 = fun(var2,alp1)
                    f3 = fun(var1,alp2)
                    f4 = fun(var2,alp2)
                    num[N,j,i] = (f1 - f2 - f3 + f4) / (4.0*vari*alpi)
    else: 
        num = np.zeros([n_int,n_dim,n_dim,n_dim,n_dim])
        for N in range(n_int):
            for i in range(n_dim):
                for j in range(n_dim):
                    for k in range(n_dim):
                        for l in range(n_dim):
                            var1 = copy.deepcopy(var)
                            var2 = copy.deepcopy(var)
                            alp1 = copy.deepcopy(alp)
                            alp2 = copy.deepcopy(alp)
                            var1[i,j] = var[i,j] - vari  
                            var2[i,j] = var[i,j] + vari
                            alp1[N,k,l] = alp[N,k,l] - alpi  
                            alp2[N,k,l] = alp[N,k,l] + alpi
                            f1 = fun(var1,alp1)
                            f2 = fun(var2,alp1)
                            f3 = fun(var1,alp2)
                            f4 = fun(var2,alp2)
                            num[N,k,l,i,j] = (f1 - f2 - f3 + f4) / (4.0*vari*alpi)
    return num
def numdiff2_4(mode, n_dim, n_int, fun, var, alp, alpi):
    if mode == 0: 
        num = np.zeros([n_int,n_int])
        for N in range(n_int):
            for M in range(n_int):
                if N==M:
                    alp1 = copy.deepcopy(alp)
                    alp3 = copy.deepcopy(alp)
                    alp1[N] = alp[N] - alpi  
                    alp3[N] = alp[N] + alpi
                    f1 = fun(var,alp1)
                    f2 = fun(var,alp)
                    f3 = fun(var,alp3)
                    num[N,N] = (f1 - 2.0*f2 + f3) / (alpi**2)
                else:
                    alp1 = copy.deepcopy(alp)
                    alp2 = copy.deepcopy(alp)
                    alp3 = copy.deepcopy(alp)
                    alp4 = copy.deepcopy(alp)
                    alp1[N] = alp[N] - alpi  
                    alp1[M] = alp[M] - alpi  
                    alp2[N] = alp[N] - alpi  
                    alp2[M] = alp[M] + alpi  
                    alp3[N] = alp[N] + alpi  
                    alp3[M] = alp[M] - alpi  
                    alp4[N] = alp[N] + alpi  
                    alp4[M] = alp[M] + alpi
                    f1 = fun(var,alp1)
                    f2 = fun(var,alp2)
                    f3 = fun(var,alp3)
                    f4 = fun(var,alp4)
                    num[N,M] = (f1 - f2 - f3 + f4) / (4.0*(alpi**2))
    elif mode == 1: 
        num = np.zeros([n_int,n_dim,n_int,n_dim])
        for N in range(n_int):
            for M in range(n_int):
                for i in range(n_dim):
                    for j in range(n_dim):
                        if N==M and i==j:
                            alp1 = copy.deepcopy(alp)
                            alp3 = copy.deepcopy(alp)
                            alp1[N,i] = alp[N,i] - alpi  
                            alp3[N,i] = alp[N,i] + alpi
                            f1 = fun(var,alp1)
                            f2 = fun(var,alp)
                            f3 = fun(var,alp3)
                            num[N,i,N,i] = (f1 - 2.0*f2 + f3) / (alpi**2)
                        else:
                            alp1 = copy.deepcopy(alp)
                            alp2 = copy.deepcopy(alp)
                            alp3 = copy.deepcopy(alp)
                            alp4 = copy.deepcopy(alp)
                            alp1[N,i] = alp[N,i] - alpi  
                            alp1[M,j] = alp[M,j] - alpi  
                            alp2[N,i] = alp[N,i] - alpi  
                            alp2[M,j] = alp[M,j] + alpi  
                            alp3[N,i] = alp[N,i] + alpi  
                            alp3[M,j] = alp[M,j] - alpi  
                            alp4[N,i] = alp[N,i] + alpi  
                            alp4[M,j] = alp[M,j] + alpi  
                            f1 = fun(var,alp1)
                            f2 = fun(var,alp2)
                            f3 = fun(var,alp3)
                            f4 = fun(var,alp4)
                            num[N,i,M,j] = (f1 - f2 - f3 + f4) / (4.0*(alpi**2))
    else: 
        num = np.zeros([n_int,n_dim,n_dim,n_int,n_dim,n_dim])
        for N in range(n_int):
            for M in range(n_int):
                for i in range(n_dim):
                    for j in range(n_dim):
                        for k in range(n_dim):
                            for l in range(n_dim):
                                if N==M and i==k and j==l:
                                    alp1 = copy.deepcopy(alp)
                                    alp3 = copy.deepcopy(alp)
                                    alp1[N,i,j] = alp[N,i,j] - alpi  
                                    alp3[N,i,j] = alp[N,i,j] + alpi
                                    f1 = fun(var,alp1)
                                    f2 = fun(var,alp)
                                    f3 = fun(var,alp3)
                                    num[N,i,j,N,i,j] = (f1 - 2.0*f2 + f3) / (alpi**2)
                                else:
                                    alp1 = copy.deepcopy(alp)
                                    alp2 = copy.deepcopy(alp)
                                    alp3 = copy.deepcopy(alp)
                                    alp4 = copy.deepcopy(alp)
                                    alp1[N,i,j] = alp[N,i,j] - alpi  
                                    alp1[M,k,l] = alp[M,k,l] - alpi  
                                    alp2[N,i,j] = alp[N,i,j] - alpi  
                                    alp2[M,k,l] = alp[M,k,l] + alpi  
                                    alp3[N,i,j] = alp[N,i,j] + alpi  
                                    alp3[M,k,l] = alp[M,k,l] - alpi  
                                    alp4[N,i,j] = alp[N,i,j] + alpi  
                                    alp4[M,k,l] = alp[M,k,l] + alpi  
                                    f1 = fun(var,alp1)
                                    f2 = fun(var,alp2)
                                    f3 = fun(var,alp3)
                                    f4 = fun(var,alp4)
                                    num[N,i,j,M,k,l] = (f1 - f2 - f3 + f4) / (4.0*(alpi**2))
    return num

def startup():
    global title, model, quiet, auto_t_step
    global mode, n_dim
    global fform, gform, rate
    global epsi, sigi, alpi, chii
    global recording, test, history, history_rec
    global acc, colour
    global t, high, test_high
    global pref
    # set up default values
    title = ""
    quiet = False
    auto_t_step = False
    model = "undefined"
    mode = 0
    n_dim = 1
    fform = False
    gform = False
    rate = False
    epsi = 0.0001
    sigi = 0.001
    alpi = 0.0001
    chii = 0.001
    recording = True
    test = False
    history = False
    history_rec = []
    acc = 0.5 # rate acceleration factor
    colour = "b"
    # initialise some key variables
    t = 0.0
    high = [[]]
    test_high = []
    pref = ["analytical", "automatic", "numerical"]
def histrec():
    global history, history_rec, line
    if history: history_rec.append(line)
           
def readv(args):
    global mode, n_dim
    if mode == 0: 
        return float(args[1])
    elif mode == 1: 
        return np.array([float(i) for i in args[1:n_dim+1]])
    elif mode == 2: 
        temp = np.array([float(i) for i in args[1:n_dim**2+1]]).reshape(n_dim,n_dim)
        return temp
def readvs(args):
    global mode, n_dim
    if mode == 0: 
        return float(args[0])
    elif mode == 1: 
        return np.array([float(i) for i in args[:n_dim]])
    elif mode == 2: 
        temp = np.array([float(i) for i in args[:n_dim**2]]).reshape(n_dim,n_dim)
        return temp
    
def printderivs():
    qprint("eps    ", eps)
    qprint("sig    ", sig)
    qprint("alp    ", alp)
    qprint("chi    ", chi)
    qprint("f      ", hm.f(eps,alp))
    qprint("dfde   ", hm.dfde(eps,alp))
    qprint("dfda   ", hm.dfda(eps,alp))
    qprint("d2fdede", hm.d2fdede(eps,alp))
    qprint("d2fdeda", hm.d2fdeda(eps,alp))
    qprint("d2fdade", hm.d2fdade(eps,alp))
    qprint("d2fdada", hm.d2fdada(eps,alp))
    qprint("g      ", hm.g(sig,alp))
    qprint("dgds   ", hm.dgds(sig,alp))
    qprint("dgda   ", hm.dgda(sig,alp))
    qprint("d2gdsds", hm.d2gdsds(sig,alp))
    qprint("d2gdsda", hm.d2gdsda(sig,alp))
    qprint("d2gdads", hm.d2gdads(sig,alp))
    qprint("d2gdada", hm.d2gdada(sig,alp))
    qprint("y      ", hm.y(eps,sig,alp,chi))
    qprint("dyde   ", hm.dyde(eps,sig,alp,chi))
    qprint("dyds   ", hm.dyds(eps,sig,alp,chi))
    qprint("dyda   ", hm.dyda(eps,sig,alp,chi))
    qprint("dydc   ", hm.dydc(eps,sig,alp,chi))
    qprint("w      ", hm.w(eps,sig,alp,chi))
    qprint("dwdc   ", hm.dwdc(eps,sig,alp,chi))

def run_general_inc(Smat, Emat, Tdt, dt, nprint, nsub):
    global eps, sig
    global start_inc
    if voigt:
        Smat[:,3:6] = Utils.rooth*Smat[:,3:6]
        Emat[:,3:6] = Utils.root2*Emat[:,3:6]
    dTdt = Tdt / float(nprint*nsub)
    ddt  = dt  / float(nprint*nsub)
    start_inc = True
    for iprint in range(nprint):
        for isub in range(nsub): 
            general_inc(Smat, Emat, dTdt, ddt)
        record(eps, sig)
    qprint("Increment complete")

def strain_inc(deps, dt):
    if rate:
        if gform: 
            strain_inc_g_r(deps, dt)
        else: 
            strain_inc_f_r(deps, dt) #default is f-form for this case, even if not set explicitly
    else:
        if gform: 
            strain_inc_g(deps, dt)
        else: 
            strain_inc_f(deps, dt) #default is f-form for this case, even if not set explicitly
def stress_inc(dsig, dt):
    if rate:
        if fform: 
            stress_inc_f_r(dsig, dt)
        else: 
            stress_inc_g_r(dsig, dt) #default is g-form for this case, even if not set explicitly
    else:
        if fform: 
            stress_inc_f(dsig, dt)
        else: 
            stress_inc_g(dsig, dt) #default is g-form for this case, even if not set explicitly
def general_inc(Smat, Emat, dTdt, dt):
    if fform:
        if rate:
            general_inc_f_r(Smat, Emat, dTdt, dt)
        else:
            general_inc_f(Smat, Emat, dTdt, dt)
    elif gform: 
        if rate:
            general_inc_g_r(Smat, Emat, dTdt, dt)
        else:
            if use_RKDP:
                general_inc_g_RKDP(Smat, Emat, dTdt, dt)
            else:
                general_inc_g(Smat, Emat, dTdt, dt)
    else: 
        error("Error in general_inc: f-form or g-form needs to be specified")

def update_f(dt, deps, dalp):
    global t, eps, sig, alp, chi
    t   = t   + dt
    eps = eps + deps
    alp = alp + dalp
    sigold = copy.deepcopy(sig)
    chiold = copy.deepcopy(chi)
    sig  = sig_f(eps,alp)
    chi  = chi_f(eps,alp)
    dsig = sig - sigold
    dchi = chi - chiold
    if hasattr(hm,"update"): hm.update(t,eps,sig,alp,chi,dt,deps,dsig,dalp,dchi)
def update_g(dt, dsig, dalp):
    global t, eps, sig, alp, chi
    t   = t   + dt
    sig = sig + dsig
    alp = alp + dalp
    epsold = copy.deepcopy(eps)
    chiold = copy.deepcopy(chi)
    eps  = eps_g(sig,alp)
    chi  = chi_g(sig,alp)
    deps = eps - epsold
    dchi = chi - chiold
    if hasattr(hm,"update"): hm.update(t,eps,sig,alp,chi,dt,deps,dsig,dalp,dchi)

def inc_mess(routine):
    global start_inc
    if start_inc: 
        qprint("Using "+routine+" for increment")
        start_inc = False
    
def general_inc_f_r(Smat, Emat, dTdt, dt): # for mode 1 only at present # updated for weights
    global eps, sig, alp, chi
    inc_mess("general_inc_f_r")
    dwdc_ = dwdc(eps,sig,alp,chi)
    d2fdede_ = d2fdede(eps,alp)
    d2fdeda_ = d2fdeda(eps,alp)
    P = np.linalg.inv(Emat + np.einsum("ij,jk->ik", Smat, d2fdede_))
    dalp  = np.einsum(ein_chi, rwt, dwdc_)*dt
    deps = np.einsum(ein_b, P, (dTdt - np.einsum("ij,mjk,mk->i", Smat, d2fdeda_, dalp)))
    update_f(dt, deps, dalp)
def general_inc_g_r(Smat, Emat, dTdt, dt): # for mode 1 only at present # updated for weights
    global eps, sig, alp, chi
    inc_mess("general_inc_g_r")
    dwdc_ = dwdc(eps,sig,alp,chi)
    d2gdsds_ = d2gdsds(sig,alp)
    d2gdsda_ = d2gdsda(sig,alp)
    Q = np.linalg.inv(Smat - np.einsum("ij,jk->ik", Emat, d2gdsds_))
    dalp  = np.einsum(ein_chi, rwt, dwdc_)*dt
    dsig = np.einsum(ein_b, Q, (dTdt + np.einsum("ij,mjk,mk->i", Emat, d2gdsda_, dalp)))
    update_g(dt, dsig, dalp)
    
def strain_inc_f_r(deps, dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("strain_inc_f_r")
    dwdc_ = dwdc(eps,sig,alp,chi)
    #print(ein_chi)
    #print(rwt)
    #print(dwdc_)
    dalp  = np.einsum(ein_chi, rwt, dwdc_)*dt
    update_f(dt, deps, dalp)
def stress_inc_g_r(dsig,dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("stress_inc_g_r")
    dwdc_ = dwdc(eps,sig,alp,chi)
    dalp  = np.einsum(ein_chi, rwt, dwdc_)*dt
    update_g(dt, dsig, dalp)

def strain_inc_g_r(deps, dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("strain_inc_g_r")
    dwdc_ = dwdc(eps,sig,alp,chi)
    d2gdsds_ = d2gdsds(sig,alp)
    d2gdsda_ = d2gdsda(sig,alp)
    if mode == 0: 
        D = -1.0 / d2gdsds_
    else: 
        D = -np.linalg.inv(d2gdsds_)
    dalp = np.einsum(ein_chi, rwt, dwdc_)*dt
    dsig = np.einsum(ein_b, D, (deps + np.einsum(ein_c, d2gdsda_, dalp)))
    update_g(dt, dsig, dalp)
def stress_inc_f_r(dsig, dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("stress_inc_f_r")
    dwdc_ = dwdc(eps,sig,alp,chi)
    d2fdede_ = d2fdede(eps,alp)
    d2fdeda_ = d2fdeda(eps,alp)
    if mode == 0: 
        C = 1.0 / d2fdede_
    else: 
        C = np.linalg.inv(d2fdede_)
    dalp = np.einsum(ein_chi, rwt, dwdc_)*dt
    deps = np.einsum(ein_b, C, (dsig - np.einsum(ein_c, d2fdeda_, dalp)))
    update_f(dt, deps, dalp)

def general_inc_f(Smat, Emat, dTdt, dt): # for mode 1 only at present # updated for weights
    global eps, sig, alp, chi
    #global P, temp, Lmatp, Lrhsp, Lmat, Lrhs, L, yo
    inc_mess("general_inc_f")
    yo = hm.y(eps,sig,alp,chi)
    dyde_ = dyde(eps,sig,alp,chi)
    dyds_ = dyds(eps,sig,alp,chi)
    dyda_ = dyda(eps,sig,alp,chi)
    rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    d2fdede_ = d2fdede(eps,alp)
    d2fdeda_ = d2fdeda(eps,alp)
    d2fdade_ = d2fdade(eps,alp)
    d2fdada_ = d2fdada(eps,alp)
    dyde_minus = dyde_ + np.einsum(ein_f,dyds_,d2fdede_) - np.einsum(ein_h,rwtdydc_,d2fdade_)
    dyda_minus = dyda_ + np.einsum(ein_g,dyds_,d2fdeda_) - np.einsum(ein_i,rwtdydc_,d2fdada_)
    P = np.linalg.inv(Emat + np.einsum("ij,jk->ik", Smat, d2fdede_))
    temp1 = np.einsum(ein_f, dyde_minus, P)
    temp2 = np.einsum("pr,rl,lnk->pnk", temp1, Smat, d2fdeda_)
    Lmatp = np.einsum(ein_j, (temp2 - dyda_minus), rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, temp1, dTdt)    
    L = solve_L(yo, Lmatp, Lrhsp)
    dalp  = np.einsum(ein_d, L, rwtdydc_)
    temp3 = np.einsum("ij,jmk,mk->i", Smat, d2fdeda_, dalp)
    deps  = np.einsum(ein_b, P, (dTdt - temp3))
    update_f(dt, deps, dalp)
def general_inc_g(Smat, Emat, dTdt, dt): #for mode 1 only at present # updated for weights
    global eps, sig, alp, chi
    #global Q, temp, Lmatp, Lrhsp, Lmat, Lrhs, L, yo
    inc_mess("general_inc_g")
    yo = hm.y(eps,sig,alp,chi)
    dyde_ = dyde(eps,sig,alp,chi)
    dyds_ = dyds(eps,sig,alp,chi)
    dyda_ = dyda(eps,sig,alp,chi)
    rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    d2gdsds_ = d2gdsds(sig,alp)
    d2gdsda_ = d2gdsda(sig,alp)
    d2gdads_ = d2gdads(sig,alp)
    d2gdada_ = d2gdada(sig,alp)
    dyds_minus = dyds_ - np.einsum(ein_f,dyde_,d2gdsds_) - np.einsum(ein_h, rwtdydc_, d2gdads_)
    dyda_minus = dyda_ - np.einsum(ein_g,dyde_,d2gdsda_) - np.einsum(ein_i, rwtdydc_, d2gdada_)
    Q = np.linalg.inv(Smat - np.einsum("ij,jk->ik", Emat, d2gdsds_))
    temp1 = np.einsum(ein_f, dyds_minus, Q)
    temp2 = np.einsum("pr,rl,lnk->pnk",temp1, Emat, d2gdsda_)
    Lmatp = np.einsum(ein_j, (-temp2 - dyda_minus), rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, temp1, dTdt)    
    L = solve_L(yo, Lmatp, Lrhsp)
    dalp  = np.einsum(ein_d, L, rwtdydc_)
    temp3 = np.einsum("ij,jmk,mk->i", Emat, d2gdsda_, dalp)
    dsig  = np.einsum(ein_b, Q, (dTdt + temp3))
    update_g(dt, dsig, dalp)
def general_inc_g_RKDP(Smat, Emat, dTdt, dt): #experimental
    global eps, sig, alp, chi, t
    #global Q, temp, Lmatp, Lrhsp, Lmat, Lrhs, L, yo
    errtol = 0.01
    inc_mess("general_inc_g_RKDP")
    def dstep(siga, alpa, Smat, Emat, dTdt):
        epsa = eps_g(siga, alpa)
        chia = chi_g(siga, alpa)
        yo    = hm.y(epsa,siga,alpa,chia)
        dyde_ = dyde(epsa,siga,alpa,chia)
        dyds_ = dyds(epsa,siga,alpa,chia)
        dyda_ = dyda(epsa,siga,alpa,chia)
        rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(epsa,siga,alpa,chia))
        d2gdsds_ = d2gdsds(siga,alpa)
        d2gdsda_ = d2gdsda(siga,alpa)
        d2gdads_ = d2gdads(siga,alpa)
        d2gdada_ = d2gdada(siga,alpa)
        dyds_minus = dyds_ - np.einsum(ein_f,dyde_,d2gdsds_) - np.einsum(ein_h, rwtdydc_, d2gdads_)
        dyda_minus = dyda_ - np.einsum(ein_g,dyde_,d2gdsda_) - np.einsum(ein_i, rwtdydc_, d2gdada_)
        Q = np.linalg.inv(Smat - np.einsum("ij,jk->ik", Emat, d2gdsds_))
        temp1 = np.einsum(ein_f, dyds_minus, Q)
        temp2 = np.einsum("pr,rl,lnk->pnk",temp1, Emat, d2gdsda_)
        Lmatp = np.einsum(ein_j, (-temp2 - dyda_minus), rwtdydc_)
        Lrhsp = acc*yo + np.einsum(ein_e, temp1, dTdt)    
        L = solve_L(yo, Lmatp, Lrhsp)
        dalp  = np.einsum(ein_d, L, rwtdydc_)
        temp3 = np.einsum("ij,jmk,mk->i", Emat, d2gdsda_, dalp)
        dsig  = np.einsum(ein_b, Q, (dTdt + temp3))
        if hasattr(hm,"update"): 
            print("RKDM routine cannot at present handle hm.update")
            error()
        return dsig, dalp
    err = 1.0
    substeps = 1
    sig0 = copy.deepcopy(sig) # make copies in case this step is abandoned
    alp0 = copy.deepcopy(alp)
    for i in range(10):
        sig = copy.deepcopy(sig0)
        alp = copy.deepcopy(alp0)
        dTdtsub = dTdt / float(substeps)
        sig_err_sum = 0.0
        sig_inc_sum = 0.0
        alp_err_sum = 0.0
        alp_inc_sum = 0.0
        for sub in range(substeps):
            sig1 = copy.deepcopy(sig)
            alp1 = copy.deepcopy(alp)
            dsig1, dalp1 = dstep(sig1, alp1, Smat, Emat, dTdtsub)
            
            sig2 = sig1 + dsig1*tfac[0,0]    
            alp2 = alp1 + dalp1*tfac[0,0]    
            dsig2, dalp2 = dstep(sig2, alp2, Smat, Emat, dTdtsub)
            
            sig3 = sig1 + dsig1*tfac[1,0] + dsig2*tfac[1,1]
            alp3 = alp1 + dalp1*tfac[1,0] + dalp2*tfac[1,1]
            dsig3, dalp3 = dstep(sig3, alp3, Smat, Emat, dTdtsub)
            
            sig4 = sig1 + dsig1*tfac[2,0] + dsig2*tfac[2,1] + dsig3*tfac[2,2]
            alp4 = alp1 + dalp1*tfac[2,0] + dalp2*tfac[2,1] + dalp3*tfac[2,2]
            dsig4, dalp4 = dstep(sig4, alp4, Smat, Emat, dTdtsub)
            
            sig5 = sig1 + dsig1*tfac[3,0] + dsig2*tfac[3,1] + dsig3*tfac[3,2] + dsig4*tfac[3,3]
            alp5 = alp1 + dalp1*tfac[3,0] + dalp2*tfac[3,1] + dalp3*tfac[3,2] + dalp4*tfac[3,3]
            dsig5, dalp5 = dstep(sig5, alp5, Smat, Emat, dTdtsub)
            
            sig6 = sig1 + dsig1*tfac[4,0] + dsig2*tfac[4,1] + dsig3*tfac[4,2] + dsig4*tfac[4,3] + dsig5*tfac[4,4]
            alp6 = alp1 + dalp1*tfac[4,0] + dalp2*tfac[4,1] + dalp3*tfac[4,2] + dalp4*tfac[4,3] + dalp5*tfac[4,4]
            dsig6, dalp6 = dstep(sig6, alp6, Smat, Emat, dTdtsub)
            
            sig_err = dsig1*terr[0] + dsig2*terr[1] + dsig3*terr[2] + dsig4*terr[3] + dsig5*terr[4] + dsig6*terr[5]    
            alp_err = dalp1*terr[0] + dalp2*terr[1] + dalp3*terr[2] + dalp4*terr[3] + dalp5*terr[4] + dalp6*terr[5]    
            sig_inc = dsig1*t5th[0] + dsig2*t5th[1] + dsig3*t5th[2] + dsig4*t5th[3] + dsig5*t5th[4] + dsig6*t5th[5]    
            alp_inc = dalp1*t5th[0] + dalp2*t5th[1] + dalp3*t5th[2] + dalp4*t5th[3] + dalp5*t5th[4] + dalp6*t5th[5]
            
            sig_err_norm = np.sqrt(np.einsum("i,i->",   sig_err,sig_err)/float(n_dim)) 
            alp_err_norm = np.sqrt(np.einsum("ij,ij->", alp_err,alp_err)/float(n_dim*n_int))    
            sig_inc_norm = np.sqrt(np.einsum("i,i->",   sig_inc,sig_inc)/float(n_dim))    
            alp_inc_norm = np.sqrt(np.einsum("ij,ij->", alp_inc,alp_inc)/float(n_dim*n_int))
            sig_err_sum += sig_err_norm
            sig_inc_sum += sig_inc_norm
            alp_err_sum += alp_err_norm
            alp_inc_sum += alp_inc_norm
#            if sig_inc_norm == 0.0: sig_inc_norm = Utils.small
#            if alp_inc_norm == 0.0: alp_inc_norm = Utils.small
#            sig_err_rel = sig_err_norm / sig_inc_norm    
#            alp_err_rel = alp_err_norm / alp_inc_norm    
            print("{:10.6f}".format(sig_inc_norm), "{:10.6f}".format(sig_err_norm), 
                  "{:21.6f}".format(alp_inc_norm), "{:10.6f}".format(alp_err_norm))    
#            errstep = np.maximum(sig_err_rel, alp_err_rel)
#            err     = np.maximum(err, errstep)
            sig = sig + sig_inc
            alp = alp + alp_inc

        if sig_inc_sum == 0.0: sig_inc_sum = Utils.small
        if alp_inc_sum == 0.0: alp_inc_sum = Utils.small
        sig_err_rel = sig_err_sum / sig_inc_sum    
        alp_err_rel = alp_err_sum / alp_inc_sum    
        print("{:10.6f}".format(sig_inc_sum), "{:10.6f}".format(sig_err_sum), "{:10.6f}".format(sig_err_rel), 
              "{:10.6f}".format(alp_inc_sum), "{:10.6f}".format(alp_err_sum), "{:10.6f}".format(alp_err_rel))    
#        err = np.maximum(sig_err_rel, alp_err_rel)
        err = sig_err_rel
        if err < errtol:
            eps = eps_g(sig, alp)
            chi = chi_g(sig, alp)
            t   = t + dt
            return
        else:
#            substeps = int(1.25*(err/errtol)**0.2) + 1
            substeps = substeps*2
            print("Substeps:",substeps)
    print("Too many trials in general_inc_RKDP")
    error()
def strain_inc_f_spec(deps): # updated for weights
    global eps, sig, alp, chi, acc
    acc_rec = acc + 0.0 # record global acc factor
    acc = 0.5
    yo = hm.y(eps,sig,alp,chi)
    rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    dyde_ = hm.dyde(eps,sig,alp,chi)
    dyds_ = hm.dyds(eps,sig,alp,chi)
    dyda_ = hm.dyda(eps,sig,alp,chi)
    d2fdede_ = hm.d2fdede(eps,alp)
    d2fdade_ = hm.d2fdade(eps,alp)
    d2fdeda_ = hm.d2fdeda(eps,alp)
    d2fdada_ = hm.d2fdada(eps,alp)
    dyde_minus = dyde_ + np.einsum(ein_f,dyds_,d2fdede_) - np.einsum(ein_h, rwtdydc_, d2fdade_)
    dyda_minus = dyda_ + np.einsum(ein_g,dyds_,d2fdeda_) - np.einsum(ein_i, rwtdydc_, d2fdada_)
    Lmatp = -np.einsum(ein_j, dyda_minus, rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, dyde_minus, deps)
    L     = solve_L(yo, Lmatp, Lrhsp)
    dalp  = np.einsum(ein_d, L, rwtdydc_)
    eps = eps + deps
    alp = alp + dalp
    sig = hm.dfde(eps,alp)
    chi = chi_f(eps,alp)
    acc = acc_rec + 0.0 # restore global acc factor
 
def strain_inc_f(deps, dt): # updated for weights
    global eps, sig, alp, chi
#    print("eps =",eps)
#    print("sig =",sig)
#    print("alp =",alp)
#    print("chi =",chi)
    inc_mess("stress_inc_f")
    yo = hm.y(eps,sig,alp,chi)
#    print(ein_dydc)
#    print(rwt)
#    print(dydc_f(chi,eps,alp))
    rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    dyde_ = dyde(eps,sig,alp,chi)
    dyds_ = dyds(eps,sig,alp,chi)
    dyda_ = dyda(eps,sig,alp,chi)
    d2fdede_ = d2fdede(eps,alp)
    d2fdade_ = d2fdade(eps,alp)
    d2fdeda_ = d2fdeda(eps,alp)
    d2fdada_ = d2fdada(eps,alp)
#    print("yo =",yo)
#    print("dydc =",dydc_)
#    print("dyde =",dyde_)
#    print("dyda =",dyda_)
#    print("d2fdade =",d2fdade_)
#    print("d2fdada =",d2fdade_)
    dyde_minus = dyde_ + np.einsum(ein_f,dyds_,d2fdede_) - np.einsum(ein_h, rwtdydc_, d2fdade_)
    dyda_minus = dyda_ + np.einsum(ein_g,dyds_,d2fdeda_) - np.einsum(ein_i, rwtdydc_, d2fdada_)
#    print("dyda_minus =",dyda_minus)
#    print("dyde_minus =",dyde_minus)
    Lmatp =         -np.einsum(ein_j, dyda_minus, rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, dyde_minus, deps)
#    print("Lmatp =",Lmatp)
#    print("Lrhsp =",Lrhsp)
    L    = solve_L(yo, Lmatp, Lrhsp)
    dalp = np.einsum(ein_d, L, rwtdydc_)
    update_f(dt, deps, dalp)
def stress_inc_g(dsig, dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("stress_inc_g")
    yo = hm.y(eps,sig,alp,chi)
    #print("yo",yo)
    rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    dyde_ = dyde(eps,sig,alp,chi)
    dyds_ = dyds(eps,sig,alp,chi)
    dyda_ = dyda(eps,sig,alp,chi)
    d2gdsds_ = d2gdsds(sig,alp)
    d2gdads_ = d2gdads(sig,alp)
    d2gdsda_ = d2gdsda(sig,alp)
    d2gdada_ = d2gdada(sig,alp)
    dyds_minus = dyds_ - np.einsum(ein_f,dyde_,d2gdsds_) - np.einsum(ein_h, rwtdydc_, d2gdads_)
    dyda_minus = dyda_ - np.einsum(ein_g,dyde_,d2gdsda_) - np.einsum(ein_i, rwtdydc_, d2gdada_)
    #print("dyda_minus",dyda_minus)
    #print("dyds_minus",dyds_minus)
    #print("dydc",dydc_)
    Lmatp = -np.einsum(ein_j, dyda_minus, rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, dyds_minus, dsig)
    #print("Lmatp",Lmatp)
    #print("Lrhsp",Lrhsp)
    L = solve_L(yo, Lmatp, Lrhsp)
    #print("L",L)
    dalp = np.einsum(ein_d, L, rwtdydc_)
    update_g(dt, dsig, dalp)
def strain_inc_g(deps, dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("strain_inc_g")
    yo = hm.y(eps,sig,alp,chi)
    dyde_ = dyde(eps,sig,alp,chi)
    dyds_ = dyds(eps,sig,alp,chi)
    dyda_ = dyda(eps,sig,alp,chi)
    # rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    rwtdydc_ = dydc(eps,sig,alp,chi)
    d2gdsds_ = d2gdsds(sig,alp)
    d2gdsda_ = d2gdsda(sig,alp)
    d2gdads_ = d2gdads(sig,alp)
    d2gdada_ = d2gdada(sig,alp)
    if mode == 0: 
        D = -1.0 / d2gdsds_
    elif mode == 1:
        D = -np.linalg.inv(d2gdsds_)
    else:
        print(d2gdsds(sig,alp))
        print(d2gdsds(sig,alp).reshape(9,9))
        D = -np.linalg.inv(d2gdsds(sig,alp).reshape(9,9))
        D = D.reshape(3,3,3,3)
    dyds_minus = dyds_ - np.einsum(ein_f,dyde_,d2gdsds_) - np.einsum(ein_h, rwtdydc_, d2gdads_)
    dyda_minus = dyda_ - np.einsum(ein_g,dyde_,d2gdsda_) - np.einsum(ein_i, rwtdydc_, d2gdada_)
    temp = np.einsum(ein_f, dyds_minus, D)
    Lmatp = np.einsum(ein_j, (-dyda_minus - np.einsum(ein_g, temp, d2gdsda_)), rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, temp, deps)
    L = solve_L(yo, Lmatp, Lrhsp)
    dalp = np.einsum(ein_d, L, rwtdydc_)
    dsig = np.einsum(ein_b, D, (deps + np.einsum(ein_c, d2gdsda_, dalp)))
    update_g(dt, dsig, dalp)
def stress_inc_f(dsig, dt): # updated for weights
    global eps, sig, alp, chi
    inc_mess("stress_inc_f")
    yo = hm.y(eps,sig,alp,chi)
    dyde_ = dyde(eps,sig,alp,chi)
    dyds_ = dyds(eps,sig,alp,chi)
    dyda_ = dyda(eps,sig,alp,chi)
    rwtdydc_ = np.einsum(ein_dydc, rwt, dydc(eps,sig,alp,chi))
    d2fdede_ = d2fdede(eps,alp)
    d2fdeda_ = d2fdeda(eps,alp)
    d2fdade_ = d2fdade(eps,alp)
    d2fdada_ = d2fdada(eps,alp)
    if mode == 0: 
        C = 1.0 / d2fdede_
    else: 
        C = np.linalg.inv(d2fdede_)
    dyde_minus = dyde_ + np.einsum(ein_f,dyds_,d2fdede_) - np.einsum(ein_h, rwtdydc_, d2fdade_)
    dyda_minus = dyda_ + np.einsum(ein_g,dyds_,d2fdeda_) - np.einsum(ein_i, rwtdydc_, d2fdada_)
    temp = np.einsum(ein_f, dyde_minus, C)
    Lmatp = np.einsum(ein_j, (np.einsum(ein_g, temp, d2fdeda_) - dyda_minus), rwtdydc_)
    Lrhsp = acc*yo + np.einsum(ein_e, temp, dsig)
    L = solve_L(yo, Lmatp, Lrhsp)
    dalp = np.einsum(ein_d, L, rwtdydc_)
    deps = np.einsum(ein_b, C, (dsig - np.einsum(ein_c, d2fdeda_, dalp)))
    update_f(dt, deps, dalp)

def record(eps, sig):
    result = True
    epso = eps_to_voigt(eps) # convert if using Voigt option
    sigo = sig_to_voigt(sig)
    if hasattr(hm,"step_print"): hm.step_print(t,epso,sigo,alp,chi)
    if recording:
        if mode == 0:
            if np.isnan(epso) or np.isnan(sigo): 
                result = False
            else: 
                rec[curr_test].append(np.concatenate(([t],[epso],[sigo],alp,chi)))
        else:
            if np.isnan(epso).any() or np.isnan(sigo).any(): 
                result = False
            else: 
                rec[curr_test].append(np.concatenate(([t],epso.flatten(),
                                                          sigo.flatten(),
                                                          alp.flatten(),
                                                          chi.flatten())))
    return result
def delrec():
    del rec[curr_test][-1]
    
def recordt(eps, sig): #record test data
    if mode == 0: 
        test_rec.append(np.concatenate(([t],[eps],[sig])))
    else: 
        test_rec.append(np.concatenate(([t],eps,sig)))

def results_print(oname):
    print("")
    if oname[-4:] != ".csv": oname = oname + ".csv"
    out_file = open(oname, 'w')
    names = names_()
    units = units_()
    for recline in rec:
        if mode == 0:
            print("{:>10} {:>14} {:>14}".format(*names))
            print("{:>10} {:>14} {:>14}".format(*units))
            out_file.write(",".join(names)+"\n")
            out_file.write(",".join(units)+"\n")
            for item in recline:
                print("{:10.4f} {:14.8f} {:14.8f}".format(*item[:3]))
                out_file.write(",".join([str(num) for num in item])+"\n")
        elif mode == 1:
            print(("{:>10} "+"{:>14} "*12).format(*names))
            print(("{:>10} "+"{:>14} "*12).format(*units))
            out_file.write(",".join(names)+"\n")
            out_file.write(",".join(units)+"\n")
            for item in recline:
                print(("{:10.4f} "+"{:14.8f} "*12).format(*item[:1+2*n_dim]))
                out_file.write(",".join([str(num) for num in item])+"\n")
    out_file.close()

def results_csv(oname):
    if oname[-4:] != ".csv": oname = oname + ".csv"
    out_file = open(oname, 'w')
    #names = names_()
    #units = units_()
    for recline in rec:
        if mode == 0:
            #out_file.write(",".join(names)+"\n")
            #out_file.write(",".join(units)+"\n")
            for item in recline:
                out_file.write(",".join([str(num) for num in item[:3]])+"\n")
        elif mode == 1:
            #out_file.write(",".join(names)+"\n")
            #out_file.write(",".join(units)+"\n")
            for item in recline:
                out_file.write(",".join([str(num) for num in item])+"\n")
    out_file.close()

def plothigh(plt,x,y,col,highl,ax1,ax2,lw=1):
    plt.plot(x, y, col, linewidth=lw)
    for item in highl: 
        plt.plot(x[item[0]:item[1]], y[item[0]:item[1]], 'r',lw)
    plt.plot(0.0,0.0)            
    # plt.set_xlabel(greek(ax1))
    # plt.set_ylabel(greek(ax2))

def greek(name):
    gnam = name.replace("eps",   r"$\epsilon$")
    gnam = gnam.replace("theta", r"$\theta$")
    gnam = gnam.replace("sig",   r"$\sigma$")
    gnam = gnam.replace("1", r"$_1$")
    gnam = gnam.replace("2", r"$_2$")
    gnam = gnam.replace("3", r"$_3$")
    gnam = gnam.replace("4", r"$_4$")
    return gnam
    
def results_graph(pname, axes, xsize, ysize):
    global test_col, high
    if pname[-4:] != ".png": pname = pname + ".png"
    names = names_()
    plt.rcParams["figure.figsize"]=(xsize,ysize)
    fig, ax = plt.subplots()
    plt.title(title)
    for i in range(len(rec)):
        recl = rec[i]
        for j in range(len(names)):
            if axes[0] == names[j]: 
                ix = j
                x = [item[j] for item in recl]
                if test:
                    xt = [item[j] for item in test_rec]
            if axes[1] == names[j]: 
                iy = j
                y = [item[j] for item in recl]
                if test:
                    yt = [item[j] for item in test_rec]
        plothigh(ax, x, y, test_col[i], high[i], nunit(ix), nunit(iy))
        if test: plothigh(ax, xt, yt, "g", test_high, "", "")
    print("Graph of",axes[1],"v.",axes[0])
    plt.title(title)
    if pname != "null.png": plt.savefig(pname)
    plt.show()

def names_():
    if hasattr(hm,"names"): 
        return hm.names
    elif mode == 0: 
        return ["t","eps","sig"]
    elif mode == 1:
        if n_dim == 1: 
            return ["t","eps","sig"]   
        elif n_dim == 2: 
            return ["t","eps1","eps2","sig1","sig2"]
        elif n_dim == 3: 
            return ["t","eps1","eps2","eps3","sig1","sig2","sig3"]
        elif n_dim == 6: 
            return ["t","eps11","eps22","eps33","gam23","gam31","gam12",
                        "sig11","sig22","sig33","tau23","tau31","tau12"]
    else:
        return ["t","eps11","eps12","eps13",
                    "eps21","eps22","eps23",
                    "eps31","eps32","eps33",
                    "sig11","sig12","sig13",
                    "sig21","sig22","sig23",
                    "sig31","sig32","sig33"]
def units_():
    if hasattr(hm,"units"): 
        return hm.units
    elif mode == 0: 
        return ["s","-","Pa"]
    elif mode == 1: 
        if n_dim == 1: 
            return ["s","-","Pa"]   
        elif n_dim == 2: 
            return ["s","-","-","Pa","Pa"]   
        elif n_dim == 3: 
            return ["s","-","-","-","Pa","Pa","Pa"]
        elif n_dim == 6: 
            return ["s","-","-","-","-","-","-","Pa","Pa","Pa","Pa","Pa","Pa"]
    else:
        return ["t","-","-","-",
                    "-","-","-",
                    "-","-","-",
                    "Pa","Pa","Pa",
                    "Pa","Pa","Pa",
                    "Pa","Pa","Pa"]
def nunit(i):
    return names_()[i] + " (" + units_()[i] + ")"

def results_plot(pname):
    global test_col, high
    if pname[-4:] != ".png": pname = pname + ".png"
    if mode == 0:
        plt.rcParams["figure.figsize"]=(6.0,6.0)
        fig, ((ax1)) = plt.subplots(1,1)
        ax1 = plt.subplot(1, 1, 1)
        plt.title(title)
        for i in range(len(rec)):
            recl = rec[i]
            e = [item[1] for item in recl]
            s = [item[2] for item in recl]
            if test:
                et = [item[1] for item in test_rec]
                st = [item[2] for item in test_rec]
            plothigh(ax1, e, s, test_col[i], high[i], nunit(1), nunit(2))
            if test: 
                plothigh(ax1, et, st, 'g', test_high, "", "")
        plt.title(title)
        if pname != "null.png": plt.savefig(pname)
        plt.show()
    elif mode == 1:
        plt.rcParams["figure.figsize"]=(8.2,8.0)
        fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2,2)
        ax1 = plt.subplot(2, 2, 1)
        plt.title(title)
        ax2 = plt.subplot(2, 2, 2)
        ax3 = plt.subplot(2, 2, 3)
        ax4 = plt.subplot(2, 2, 4)
        plt.subplots_adjust(wspace=0.5,hspace=0.3)
        for i in range(len(rec)):
            recl = rec[i]
            e1 = [item[1] for item in recl]
            e2 = [item[2] for item in recl]
            s1 = [item[3] for item in recl]
            s2 = [item[4] for item in recl]
            if test:
                et1 = [item[1] for item in test_rec]
                et2 = [item[2] for item in test_rec]
                st1 = [item[3] for item in test_rec]
                st2 = [item[4] for item in test_rec]
            plothigh(ax1, s1, s2, test_col[i], high[i], nunit(3), nunit(4))
            if test: 
                plothigh(ax1, st1, st2, 'g', test_high, nunit(3), nunit(4))
            plothigh(ax2, e2, s2, test_col[i], high[i], nunit(2), nunit(4))
            if test: 
                plothigh(ax2, et2, st2, 'g', test_high, nunit(2), nunit(4))
            plothigh(ax3, s1, e1, test_col[i], high[i], nunit(3), nunit(1))
            if test: 
                plothigh(ax3, st1, et1, 'g', test_high, nunit(3), nunit(1))
            plothigh(ax4, e2, e1, test_col[i], high[i], nunit(2), nunit(1))
            if test: 
                plothigh(ax4, et2, et1, 'g', test_high, nunit(2), nunit(1))
        if pname != "null.png": plt.savefig(pname)
        plt.show()

def results_plotCS(pname):
    global test_col, high
    if pname[-4:] != ".png": pname = pname + ".png"
    plt.rcParams["figure.figsize"]=(13.0,8.0)
    #fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(2,3)
    ax2 = plt.subplot(2, 3, 2)
    plt.title(title)
    ax3 = plt.subplot(2, 3, 3)
    ax4 = plt.subplot(2, 3, 4)
    ax5 = plt.subplot(2, 3, 5)
    ax6 = plt.subplot(2, 3, 6)
    plt.subplots_adjust(wspace=0.45,hspace=0.3)
    ax4.set_xscale("log")
    ax4.invert_yaxis()
    ax5.invert_yaxis()
    ax6.invert_yaxis()
    for i in range(len(rec)):
        recl = rec[i]
        e1 = [item[1] for item in recl]
        e2 = [item[2] for item in recl]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        if test:
            et1 = [item[1] for item in test_rec]
            et2 = [item[2] for item in test_rec]
            st1 = [item[3] for item in test_rec]
            st2 = [item[4] for item in test_rec]
        if hasattr(hm,"M"):
            maxp = np.max(s1)
            maxq = np.max(s2)
            maxp = np.min([maxp,maxq*1.3/hm.M])
            minp = np.max(s1)
            minq = -np.min(s2)
            minp = np.min([minp,minq*1.3/hm.M])
            cslp = np.array([minp, 0.0, maxp])
            cslq = np.array([-minp*hm.M, 0.0, maxp*hm.M])
            ax2.plot(cslp,cslq,"red")
        plothigh(ax2, s1, s2, test_col[i], high[i], nunit(3), nunit(4))
        if test: 
            plothigh(ax2, st1, st2, 'g', test_high, nunit(3), nunit(4))
        plothigh(ax3, e2, s2, test_col[i], high[i], nunit(2), nunit(4))
        if test: 
            plothigh(ax3, et2, st2, 'g', test_high, nunit(2), nunit(4))
        plothigh(ax4, s1, e1, test_col[i], high[i], nunit(3)+" (log scale)", nunit(1))
        if test: 
            plothigh(ax4, st1, et1, 'g', test_high, nunit(3), nunit(1))
        plothigh(ax5, s1, e1, test_col[i], high[i], nunit(3), nunit(1))
        if test: 
            plothigh(ax5, st1, et1, 'g', test_high, nunit(3), nunit(1))
        plothigh(ax6, e2, e1, test_col[i], high[i], nunit(2), nunit(1))
        if test: 
            plothigh(ax6, et2, et1, 'g', test_high, nunit(2), nunit(1))
    if pname != "null.png": plt.savefig(pname)
    plt.show()

def modestart(): # initialise the calculation for the specified mode
    global ein_a, ein_b, ein_c, ein_d, ein_e, ein_f, ein_g, ein_h, ein_i, ein_j
    global sig, sig_inc, sig_targ, sig_cyc, dsig
    global eps, eps_inc, eps_targ, eps_cyc, deps
    global chi, alp, dalp
    global yo
    global n_int, n_y
    global rwt, ein_chi, ein_dydc

    #print(hm.n_int)
    n_int = max(1,hm.n_int)
    n_y = max(1,hm.n_y)
    if mode == 0:           # typical usage below
        ein_a = ",->"       # sig eps -> (scalar)
        ein_b = ",->"       # D deps -> dsig
        ein_c = "m,m->"     # d2fdeda dalp ->
        ein_d = "N,Nm->m"   # L dydc ->
        ein_e = "N,->N"     # dyde deps ->
        ein_f = "N,->N"     # dyde C ->
        ein_g = "N,n->Nn"   # temp d2fdeda ->
        ein_h = "Nm,m->N"   # dydc d2fdade -> dyde
        ein_i = "Nm,mn->Nn" # dydc d2fdada -> dyda
        ein_j = "Nn,Mn->NM" # dyda dydc ->
        sig = 0.0
        eps = 0.0
        alp = np.zeros(n_int)
        chi = np.zeros(n_int)
        if hasattr(hm, "rwt"):
            rwt = hm.rwt
        else:
            rwt  = np.ones(n_int)
        ein_chi = "N,N->N"
        ein_dydc = "N,MN->MN"
        dsig = 0.0
        deps = 0.0
        dalp = np.zeros(n_int)
        sig_inc  = 0.0
        sig_targ = 0.0
        sig_cyc  = 0.0
        eps_inc  = 0.0
        eps_targ = 0.0
        eps_cyc  = 0.0
        yo = np.zeros(n_y)
    elif mode == 1:
        ein_a = "i,i->"
        ein_b = "ki,i->k"
        ein_c = "imj,mj->i"
        ein_d = "N,Nmi->mi"
        ein_e = "Ni,i->N"
        ein_f = "Nj,jk->Nk"
        ein_g = "Nk,knl->Nnl"
        ein_h = "Nmi,mij->Nj"
        ein_i = "Nmi,minj->Nnj"
        ein_j = "Nni,Mni->NM"
        sig = np.zeros(n_dim)
        eps = np.zeros(n_dim)
        alp = np.zeros([n_int,n_dim])
        chi = np.zeros([n_int,n_dim])
        if hasattr(hm, "rwt"):
            rwt = hm.rwt
        else:
            rwt  = np.ones(n_int)
        ein_chi = "N,Ni->Ni"
        ein_dydc = "N,MNi->MNi"
        dsig = np.zeros(n_dim)
        deps = np.zeros(n_dim)
        dalp = np.zeros([n_int,n_dim])
        sig_inc  = np.zeros(n_dim)
        sig_targ = np.zeros(n_dim)
        sig_cyc  = np.zeros(n_dim)
        eps_inc  = np.zeros(n_dim)
        eps_targ = np.zeros(n_dim)
        eps_cyc  = np.zeros(n_dim)
        yo = np.zeros(n_y)
    elif mode == 2:
        ein_a = "ij,ij->"
        ein_b = "klij,ij->kl"
        ein_c = "ijmkl,mkl->ij"
        ein_d = "N,Nmij->mij"
        ein_e = "Nij,ij->N"
        ein_f = "Nij,ijkl->Nkl"
        ein_g = "Nkl,klnij->Nnij"
        ein_h = "Nmij,mijkl->Nkl"
        ein_i = "Nmij,mijnkl->Nnkl"
        ein_j = "Nnij,Mnij->NM"
        sig = np.zeros([n_dim,n_dim])
        eps = np.zeros([n_dim,n_dim])
        alp = np.zeros([n_int,n_dim,n_dim])
        chi = np.zeros([n_int,n_dim,n_dim])
        if hasattr(hm, "rwt"):
            rwt = hm.rwt
        else:
            rwt  = np.ones(n_int)
        ein_chi = "N,Nij->Nij"
        ein_dydc = "N,MNij->MNij"
        dsig = np.zeros([n_dim,n_dim])
        deps = np.zeros([n_dim,n_dim])
        dalp = np.zeros([n_int,n_dim,n_dim])
        sig_inc  = np.zeros([n_dim,n_dim])
        sig_targ = np.zeros([n_dim,n_dim])
        sig_cyc  = np.zeros([n_dim,n_dim])
        eps_inc  = np.zeros([n_dim,n_dim])
        eps_targ = np.zeros([n_dim,n_dim])
        eps_cyc  = np.zeros([n_dim,n_dim])
        yo = np.zeros(n_y)
    else:
        error("Mode not recognised:" + mode)        

def modestart_short():
    global rwt, ein_chi, ein_dydc
    print("Modestart: n_int =",n_int,", n_dim =",n_dim)
    an_int = max(1,n_int)
    if mode == 0:
        sig = 0.0
        eps = 0.0
        alp = np.zeros(an_int)
        chi = np.zeros(an_int)
        rwt = np.ones(an_int)
        ein_chi = "N,N->N"
        ein_dydc = "N,MN->MN"
        ein_inner = ",->"
        ein_contract = ",->"
    elif mode == 1:
        sig = np.zeros(n_dim)
        eps = np.zeros(n_dim)
        alp = np.zeros([an_int,n_dim])
        chi = np.zeros([an_int,n_dim])
        rwt = np.ones(an_int)
        ein_chi = "N,Ni->Ni"
        ein_dydc = "N,MNi->MNi"
        ein_inner = "ij,jk->ik"
        ein_contract = "i,i->"
    elif mode == 2:
        sig = np.zeros([n_dim,n_dim])
        eps = np.zeros([n_dim,n_dim])
        alp = np.zeros([an_int,n_dim,n_dim])
        chi = np.zeros([an_int,n_dim,n_dim])
        rwt = np.ones(an_int)
        ein_chi = "N,Nij->Nij"
        ein_dydc = "N,MNij->MNij"
        ein_inner = "ijkl,klmn->ijmn"
        ein_contract = "ij,ij->"
    else:
        print("Mode not recognised:", mode)
        error()
    return ein_inner, ein_contract, sig, eps, chi, alp

# Drive - driving routine for hyperplasticity models
def drive(arg="hyper.dat"): 
    print("")
    print("+-----------------------------------------------------------------------------+")
    print("| HyperDrive: driving routine for hyperplasticity models                      |")
    print("| (c) G.T. Houlsby, 2018-2021                                                 |")
    print("|                                                                             |")
    print("| \x1b[1;31mThis program is provided in good faith, but with no warranty of correctness\x1b[0m |")
    print("+-----------------------------------------------------------------------------+")
    
    print("Current directory: " + os.getcwd())
    if os.path.isfile("HyperDrive.cfg"):
        input_file = open("HyperDrive.cfg", 'r')
        last = input_file.readline().split(",")
        drive_file_name = last[0]
        check_file_name = last[1]
        input_file.close()
    else:
        drive_file_name = "hyper.dat"
        check_file_name = "unknown"
    if arg != "hyper.dat":
        drive_file_name = arg
    input_found = False
    count = 0
    while not input_found:
        count += 1
        if count > 5: error("Too many tries")
        drive_file_name = read_def("Enter input filename", drive_file_name)
        if drive_file_name[-4:] != ".dat": 
            drive_file_name = drive_file_name + ".dat"
        if os.path.isfile(drive_file_name):
            print("Reading from file: " + drive_file_name)
            print("")
            input_file = open(drive_file_name, 'r')
            input_found = True
        else:
            print("File not found:", drive_file_name)
    last_file = open("HyperDrive.cfg", 'w')
    last_file.writelines([drive_file_name+","+check_file_name])
    last_file.close()
    
    startup()
    process(input_file)
    
    print("End of (processing from) file: " + drive_file_name)
    input_file.close()

#############
###################
############################
###################################
###### ADDITION #############################
# XXX
def returnrecext(args=[]):
    global rec
    ultest = rec
    print('Your file now contains last record')
    return(ultest)

def c(): check()
def d(): drive()
def cd(direc): 
    os.chdir(direc)
    print("Current directory: " + os.getcwd())
def pd(): 
    print("Current directory: " + os.getcwd())

if __name__ == "__main__":
    print("\nHyperDrive routines loaded")
    print("(c) G.T. Houlsby 2018-2021\n")
    print("Usage:")
    print("  drive('datafile') - run HyperDrive taking input from 'datafile.dat'")
    print("  drive()           - run HyperDrive taking input from last datafile used")
    print("  check('model')    - run HyperCheck on code in 'model.py'")
    print("  check()           - run HyperCheck on last code tested")
    print("  cd('dir')         - change to directory 'dir'")
    print("  pd()              - print current directory\n")
    print("Current directory: " + os.getcwd())