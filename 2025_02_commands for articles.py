# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:01:41 2025

@author: lsimonin

newp.plotD([],tmd1+tmd2d+tmd3d+tmd4d+tmd5d+tmd6d,['k','b','g','r','c','m'],[500,-0,500,0,0.25,-0.00,0.04,-0.06],['1','2','3','4','5','6'])
newp.plotD([],tmd1+tmd2l+tmd3l+tmd4l+tmd5l+tmd6l,['k','b','g','r','c','m'],[500,-0,500,0,0.25,-0.00,0.04,-0.06],['1','2','3','4','5','6'])
newp.plotU([],tmu1+tmu2d+tmu3d+tmu4d+tmu5d+tmu6d,['k','b','g','r','c','m'],[400,-0,250,0,0.05,-0.00,0.04,-0.06],['1','2','3','4','5','6'],yt1=4,yt2=4)

newp.plotD([],tmd1+tmd1bis,['k','k--'],[500,-0,500,0,0.25,-0.00,0.025,-0.08],['     m=0.7\n'+'non-linear elastic','     m=0\n'+'linear elastic'])
newp.plotU([],tmu1+tmu1bis,['k','k--'],[500,-0,500,0,0.05,-0.00,0.025,-0.08],['     m=0.7\n'+'non-linear elastic','     m=0\n'+'linear elastic'])

"""

from HyperDrive import Commands as H
import newplot as newp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
TC = __import__('2025_02_test_commands')

H.init()
H.title(["HySand test run for theory article"])
H.g_form()
H.acc([0.5])
H.mode([1,2])
H.quiet()

prec = 20
bug = 0     # bug detector

#Parameters
bet = 1.0   # Dilation ratio at densest state
A0 = 40     # Anisotropy rate development at densest state
r = 0.3     # Fusion exponent between friction and consolidation
h = 1000    # Hardening modulus for most inner yield surface at densest state
b = 3       # Hardening exponent
a0 = 0      # Initial anisotropy
NN = 1      # number of yield surfaces
K = 50000   # Bulk modulus at pr
G = 40000   # Shear modulus at pr
ind = 0.7   # Pressure index for non-linear elasticity and consolidation
phic = 33   # Friction angle at critical state
Be = 2.05   # Specific volume of loosest state at pr
Ga = 1.95   # Specific volume of critical state at pr
De = 1.7    # Specific volume of densest state at pr
l0 = 0.007  # Normal consolidation slope in loosest state
l1 = 0.003  # Normal consolidation slope in densest state

# For monotonic
eqd = 0.3   # 30% deviatoric strain command for monotonic drained tests
equ = 0.05  # 5% deviatoric strain command for monotonic undrained tests
# For cyclic
qc = 30             # Deviatoric stress amplitude
ncyc = 10              # Number of cycles
pcut = 5            # p at the end of a quarter cycle at which the test is stopped if reached
eqcut = 0.05        # deviatoric strain at which test is stopped if reached at end of half cycle
eqprec = 0.000001   # deviatoric strain step size
firstd = 0          # indicator of first drained cycle

anidil = 0
if anidil == 1:
    H.model(["2025_02_anidil"])
    NN = 1
    pt,qt=100,0 # Initial stress values
    vt = 2.0
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd4,tmd4l = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    # ##### Monotonic undrained #####
    wu4,tmu4l = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    vt = 1.9
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd4,tmd4m = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu4,tmu4m = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    vt = 1.8
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd4,tmd4d = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu4,tmu4d = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    ##### Cyclic stress-controlled undrained #####
    wc4,tcui4,TCUINh4,TCUINm14,TCUINm24,TCUINl4 = TC.TCUS_test(pt,qt,vt,qc,ncyc,pcut,eqcut,prec,eqprec,params,bug,firstd)
# plotD([],[tmd4l[0][100:]]+[tmd4m[0][100:]]+[tmd4d[0][100:]],['k-.','r','b--'],[500,-0,500,0,0.25,-0.00,0.04,-0.06],['v'+'$_t$'+' = 2.0','v'+'$_t$'+' = 1.9','v'+'$_t$'+' = 1.8'],lab=['',''],ncl=1,lw=1.75,ttl='a) Anisotropic dilation')
# plotU([],[tmu4l[0][100:]]+[tmu4m[0][100:]]+[tmu4d[0][100:]],['k-.','r','b--'],[400,-0,250,0,0.05,-0.00,0.025,-0.1],['v'+'$_t$'+' = 2.0','v'+'$_t$'+' = 1.9','v'+'$_t$'+' = 1.8'],lab=['',''],yt1=4,yt2=4,ncl=1,lw=1.75,ttl='a) Anisotropic dilation')
# plotU([],[tcui4[0][100:]],['k'],[45,-45,120,0,0.01,-0.01,0.025,-0.1],['v'+'$_t$'+' = 1.8'],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=4,ncl=1,lw=1.2,ttl='a) Anisotropic dilation')

HySand1 = 0
if HySand1 == 1:            # One Surface HySand
    H.model(["2025_02_HySand"])
    NN = 1
    pt,qt=100,0 # Initial stress values
    vt = 2.0
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd5,tmd5l = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu5,tmu5l = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    vt = 1.9
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd5,tmd5m = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu5,tmu5m = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    vt = 1.8
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd5,tmd5d = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu5,tmu5d = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    ##### Cyclic stress-controlled undrained #####
    wc5,tcui5,TCUINh5,TCUINm15,TCUINm25,TCUINl5 = TC.TCUS_test(pt,qt,vt,qc,ncyc,pcut,eqcut,prec,eqprec,params,bug,firstd)    
# plotD([],[tmd5l[0][100:]]+[tmd5m[0][100:]]+[tmd5d[0][100:]],['k-.','r','b--'],[500,-0,500,0,0.25,-0.00,0.04,-0.06],['v'+'$_t$'+' = 2.0','v'+'$_t$'+' = 1.9','v'+'$_t$'+' = 1.8'],lab=['',''],ncl=1,lw=1.75,ttl='b) plus Consolidation')
# plotU([],[tmu5l[0][100:]]+[tmu5m[0][100:]]+[tmu5d[0][100:]],['k-.','r','b--'],[400,-0,250,0,0.05,-0.00,0.025,-0.1],['v'+'$_t$'+' = 2.0','v'+'$_t$'+' = 1.9','v'+'$_t$'+' = 1.8'],lab=['',''],yt1=4,yt2=4,ncl=1,lw=1.75,ttl='b) plus Consolidation')
# plotU([],[tcui5[0][100:]],['k'],[45,-45,120,0,0.01,-0.01,0.025,-0.1],['v'+'$_t$'+' = 1.8'],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=4,ncl=1,lw=1.2,ttl='b) plus Consolidation')

HySand10 = 0
if HySand10 == 1:            # Ten Surfaces HySand
    H.model(["2025_02_HySand"])
    NN = 10
    pt,qt = 100,0 # Initial stress values
    vt = 2.0
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd6,tmd6l = TC.drained_test(pt,qt,vt,eqd,prec*10,params,bug)
    ##### Monotonic undrained #####
    wu6,tmu6l = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    vt = 1.9
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd6,tmd6m = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu6,tmu6m = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    vt = 1.8
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    ##### Monotonic drained #####
    wd6,tmd6d = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    ##### Monotonic undrained #####
    wu6,tmu6d = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
    ##### Cyclic stress-controlled undrained #####
    wc6,tcui6,TCUINh6,TCUINm16,TCUINm26,TCUINl56= TC.TCUS_test(pt,qt,vt,qc,ncyc,pcut,eqcut,prec,eqprec,params,bug,firstd)    
# plotD([],[tmd6l[0][100:]]+[tmd6m[0][100:]]+[tmd6d[0][100:]],['k-.','r','b--'],[500,-0,500,0,0.25,-0.00,0.04,-0.06],['v'+'$_t$'+' = 2.0','v'+'$_t$'+' = 1.9','v'+'$_t$'+' = 1.8'],lab=['',''],ncl=1,lw=1.75,ttl='c) plus Multi-surface')
# plotU([],[tmu6l[0][100:-12]]+[tmu6m[0][100:]]+[tmu6d[0][100:]],['k-.','r','b--'],[400,-0,250,0,0.05,-0.00,0.025,-0.1],['v'+'$_t$'+' = 2.0','v'+'$_t$'+' = 1.9','v'+'$_t$'+' = 1.8'],lab=['',''],yt1=4,yt2=4,ncl=1,lw=1.75,ttl='c) plus Multi-surface')
# plotU([],[tcui6[0][100:]],['k'],[45,-45,120,0,0.01,-0.01,0.025,-0.1],['v'+'$_t$'+' = 1.8'],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=4,ncl=1,lw=1.2,ttl='c) plus Multi-surface')

# Simulation of strain cycling with increasing strain amplitude at a constant mean effective stress
st_cy = 0
if st_cy == 1:
    vt = 1.8
    pt = 100
    qt = 0
    H.model(["2025_02_anidil"])
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([1,0,0,1,0,0,0,0,0.0,0.01,1,10,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,0.01,1,100,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,-0.02,1,100,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,0.03,1,100,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,-0.04,1,25,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,0.05,1,25,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,-0.06,1,25,prec])
    H.end()
    tcd4=H.returnrec()
    H.model(["2025_02_HySand"])
    NN = 1
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([1,0,0,1,0,0,0,0,0.0,0.01,1,10,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,0.01,1,100,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,-0.02,1,100,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,0.03,1,100,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,-0.04,1,25,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,0.05,1,25,prec])
    H.general_inc([1,0,0,0,0,0,0,1,0.0,-0.06,1,25,prec])
    H.end()
    tcd5=H.returnrec()    
# plotD([],[tcd4[0][100:]]+[tcd5[0][100:]],['k','b-.'],[180,-120,500,0,0.03,-0.03,0.02,-0.01],[' Anisotropic\n    dilation','      plus\nConsolidation'],lab=['',''],lw=1.5,xt1=4,yt1=3,ori='h')

prec = 10

# Simulation of preloaded test
pre = 1
if pre == 1:
    H.model(["2025_02_HySand"])
    NN=1
    pt,qt,vt = 100,0,1.8
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.01,1,100,5*prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.04,1,100,prec])
    H.end()
    tmupl0=H.returnrec()
    # OCR = 3
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    lstate = H.returnstate()
    eps_t = lstate[1]
    H.general_inc([0,1,1,0,0,0,0,0,0.0,200,1,100,prec])
    H.general_inc([0,1,1,0,0,0,0,0,0.0,-200,1,100,prec])
    lstate = H.returnstate()
    eps_p = lstate[1]
    print("\np=",lstate[3])
    print("Such that upon shear loading v=",vt*np.exp(-(eps_p-eps_t)))
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.01,1,50,2*prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.04,1,100,prec])
    H.end()
    tmupl1=H.returnrec()
    # Drained preshearing
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([1,0,0,1,0,0,0,0,0.0,50,1,100,prec])
    H.general_inc([1,0,0,1,0,0,0,0,0.0,-50,1,100,prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.01,1,50,2*prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.04,1,100,prec])
    H.end()
    tmupl2=H.returnrec()
# plotU([],[tmupl0[0][100:]]+[tmupl1[0][300:]]+[tmupl2[0][300:]],['k','b-.','r--',],[250,-0,200,0,0.03,-0.00,0.02,-0.01],['    Normal\nconsolidation','   OCR = 3','Loading after\n    drained\n preshearing'],lab=['',''],lw=1.5,xt1=4,xt2=3,ori='h')

###############################################################################
###############################################################################
###############################################################################
###############################   PERFORMANCE   ###############################
###############################################################################
###############################################################################
###############################################################################

H.model(["2025_02_HySand"])
prec = 20
[NN,vt ,  K  ,  G  ,ind ,phic, bet ,  Be , Ga ,  De ,  lb  ,  ld  , A0 ,  r  ,  h ,b,a0] = \
[10,1.8,49300,41700,0.73,33.9,0.995,2.054,1.98,1.677,0.0072,0.0027,55.2,0.263,1334,3,0 ]
params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]

# vt = 1.677
# params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
# H.const(params)
# H.start()
# H.init_stress([10**-6,0])
# H.general_inc([0,1,1,0,0,0,0,0,0.0,100,1,100,prec])
# H.end()
# isotestd = H.returnrec()

##### ISOTROPIC #####
isot = 0
if isot == 1:
    prec = 20
    def iso_load(vt):
        pt,qt = 50,0
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        H.const(params)
        H.start()
        TC.start_t(pt,qt,vt,params,prec)
        H.general_inc([0,1,1,0,0,0,0,0,0.0,750,1,100,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,-750,1,100,prec])
        H.end()
    iso = [0]*7
    vtl = [0,1.974,1.823,1.69,1.963,1.824,1.68]
    for i in [1,2,3]:
        print('\nStarting simulation of ISO%i'%i)
        iso_load(vtl[i])
        iso[i] = H.returnrec()
    def iso_loadc(vt):
        pt,qt = 50,0
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        H.const(params)
        H.start()
        TC.start_t(pt,qt,vt,params,prec)
        H.general_inc([0,1,1,0,0,0,0,0,0.0,350,1,35,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,-200,1,20,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,300,1,30,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,-200,1,20,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,300,1,30,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,-200,1,20,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,300,1,30,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,-200,1,20,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,300,1,30,prec])
        H.general_inc([0,1,1,0,0,0,0,0,0.0,-750,1,75,prec])
        H.end()
    for i in [4,5,6]:
        print('\nStarting simulation of ISO%i'%i)
        iso_loadc(vtl[i])
        iso[i] = H.returnrec()
    
    # newp.plotISO(ISO[1],[iso[1][0][100:]],['k:','b'],[100,-0,900,0,0.05,-0.0,0.03,-0.0],['data','simulation'],lab=['',''],yt1=6,xt1=9,ncl=1,lw=1.2,ttl='a) ISO1')
    # newp.plotISO(ISO[2],[iso[2][0][100:]],['k:','b'],[100,-0,900,0,0.05,-0.0,0.02,-0.0],['data','simulation'],lab=['',''],yt1=4,xt1=9,ncl=1,lw=1.2,ttl='b) ISO2')
    # newp.plotISO(ISO[3],[iso[3][0][100:]],['k:','b'],[100,-0,900,0,0.05,-0.0,0.015,-0.0],['data','simulation'],lab=['',''],yt1=3,xt1=9,ncl=1,lw=1.2,ttl='c) ISO3')
    # newp.plotISO(ISO[4],[iso[4][0][100:]],['k:','b'],[100,-0,900,0,0.05,-0.0,0.03,-0.0],['data','simulation'],lab=['',''],yt1=6,xt1=9,ncl=1,lw=1.2,ttl='d) ISO4')
    # newp.plotISO(ISO[5],[iso[5][0][100:]],['k:','b'],[100,-0,900,0,0.05,-0.0,0.02,-0.0],['data','simulation'],lab=['',''],yt1=4,xt1=9,ncl=1,lw=1.2,ttl='e) ISO5')
    # newp.plotISO(ISO[6],[iso[6][0][100:]],['k:','b'],[100,-0,900,0,0.05,-0.0,0.015,-0.0],['data','simulation'],lab=['',''],yt1=3,xt1=9,ncl=1,lw=1.2,ttl='f) ISO6')

##### OEDOMETRIC #####
oedot = 0
if oedot == 1:
    prec = 1
    nit = 5000
    def oedo_load(vt):
        pt,qt = 10**-6,0
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        H.const(params)
        H.start()
        H.init_stress([pt,qt])
        H.general_inc([1,2/3,0,0,0,0,2,-3,400,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,-399.9,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,399.9,0.00,1,nit,prec])
        H.end()
    oe = [0]*13
    vtl = [0,2.0386,2.0294,1.9899,1.9711,1.9462,1.9077,1.8462,1.83285,1.8076,1.7773,1.7404,1.7215]
    for i in range(1,13):
        print('\nStarting simulation of OE%i'%i)
        oedo_load(vtl[i])
        oe[i]=H.returnrec()
    # newp.plotOE([OE[i][0] for i in range(1,13)],[],[2.039,2.029,1.99,1.971,1.946,1.908,1.846,1.833,1.808,1.777,1.740,1.721],['k:','c--','b','g-.','r:','m--','k','c-.','b:','g--','r','m-.'],[1000,0.1,1.05,0.65],['OE1','OE2','OE3','OE4','OE5','OE6','OE7','OE8','OE9','OE10','OE11','OE12'],lab=['',''],yt1=9,yt2=6,xt1=10,xt2=4,ncl=1,lw=1.5,ttl='Oedometric tests: data')
    # newp.plotOE([],[oe[i][0] for i in range(1,13)],[2.039,2.029,1.99,1.971,1.946,1.908,1.846,1.833,1.808,1.777,1.740,1.721],['k:','c--','b','g-.','r:','m--','k','c-.','b:','g--','r','m-.'],[1000,0.1,1.05,0.65],['OE1','OE2','OE3','OE4','OE5','OE6','OE7','OE8','OE9','OE10','OE11','OE12'],lab=['',''],yt1=9,yt2=6,xt1=10,xt2=4,ncl=1,lw=1.5,ttl='Oedometric tests: HySand_base simulation')

oedoct = 0
if oedoct == 1:
    # Cyclic Oedometric tests  
    prec = 1
    nit = 2500
    def oedo_load_c(vt):
        pt,qt = 10**-6,0
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        H.const(params)
        H.start()
        H.init_stress([pt,qt])
        H.general_inc([1,2/3,0,0,0,0,2,-3,20,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,-19,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,55,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,-55,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,140,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,-140,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,400,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,-400,0.00,1,nit,prec])
        H.general_inc([1,2/3,0,0,0,0,2,-3,400,0.00,1,nit,prec])
        H.end()
    oec = [0]*5
    vtl = [0,1.9976,1.8898,1.7919,1.7563]
    for i in range(1,5):
        print('\nStarting simulation of OEC%i'%i)
        oedo_load_c(vtl[i])
        oec[i]=(H.returnrec())
    # plotOE([OEC[i][0] for i in [1,2,3,4]],[],[1.998,1.89,1.792,1.756],['k:','b','g-.','r--'],[1000,0.1,1.0,0.7],['OE1','OE2','OE3','OE4'],lab=['',''],yt1=7,yt2=6,xt1=10,xt2=4,ncl=1,lw=1.5,ttl='Oedometric cyclic tests: data')
    # plotOE([],[oec[i][0] for i in [1,2,3,4]],[1.998,1.89,1.792,1.756],['k:','b','g-.','r--'],[1000,0.1,1.0,0.7],['OE1','OE2','OE3','OE4'],lab=['',''],yt1=7,yt2=6,xt1=10,xt2=4,ncl=1,lw=1.5,ttl='Oedometric cyclic tests: HySand_base ')
    # plotU([],oec[3],['b','g-.','r--'],[250,-50,250,0,0.015,-0.00,0.02,-0.08],['OEC3'],lab=['',''],yt1=7,yt2=6,xt1=5,xt2=3,ncl=1,lw=1.5,ttl='Oedometric cyclic test: HySand_base',ori='h')
##### MONOTONIC DRAINED #####
tmdt = 0
if tmdt == 1:
    eqd = 0.25
    prec = 20
    bug = 0
    TMDlist = [[],
            [1.996,50],[1.975,100],[1.975,200],[1.970,300],[1.960,400],
            [1.880,50],[1.862,100],[1.859,200],[1.848,300],[1.847,400],
            [1.840,50],[1.819,100],[1.824,200],[1.822,300],[1.814,400],
            [1.743,50],[1.758,100],[1.748,200],[1.734,300],[1.753,400],
            [1.734,50],[1.735,100],[1.706,200],[1.697,300],[1.718,400],]
    tmd = [0]*26 
    for i in range(1,26):
        print('\nStarting simulation of TMD%i'%i)
        qt = 0
        vt,pt = TMDlist[i] # Initial stress values
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wd,tmd[i] = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
    # plotD([TMD[i][0] for i in [11,12,13,14,15]],[],['k:','b--','g','r-.','m:'],[1250,-0,500,0,0.25,-0.00,0.025,-0.1],['p'+'$_t$'+'=%i kPa'%i for i in [50,100,200,300,400]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='a) Drained for varied p'+'$_t$'+': data')
    # plotD([TMD[i][0] for i in [3,8,13,18,23]],[],['k:','b--','g','r-.','m:'],[1250,-0,500,0,0.25,-0.00,0.025,-0.1],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.975,1.859,1.824,1.748,1.706]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='b) Drained for varied v'+'$_t$'+' (p'+'$_t$'+'=200 kPa): data')
    # plotD([TMD[i][0] for i in [3,8,13,18,23]],[],['k:','b--','g','r-.','m:'],[250,-0,500,0,0.0025,-0.00,0.0025,-0.0],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.975,1.859,1.824,1.748,1.706]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='c) Zoom on drained for varied v'+'$_t$'+': data')
    # plotD([],[tmd[i][0][100:] for i in [11,12,13,14,15]],['k:','b--','g','r-.','m:'],[1250,-0,500,0,0.25,-0.00,0.025,-0.1],['p'+'$_t$'+'=%i kPa'%i for i in [50,100,200,300,400]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='d) Drained for varied p'+'$_t$'+': HySand_base')
    # plotD([],[tmd[i][0][100:] for i in [3,8,13,18,23]],['k:','b--','g','r-.','m:'],[1250,-0,500,0,0.25,-0.00,0.025,-0.1],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.975,1.859,1.824,1.748,1.706]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='e) Drained for varied v'+'$_t$'+': HySand_base')
    # plotD([],[tmd[i][0][100:] for i in [3,8,13,18,23]],['k:','b--','g','r-.','m:'],[250,-0,500,0,0.0025,-0.00,0.0025,-0.0],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.975,1.859,1.824,1.748,1.706]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='f) Zoom for varied v'+'$_t$'+': HySand_base')
    
##### MONOTONIC UNDRAINED #####
tmut = 0
if tmut == 1:
    eqd = 0.05
    prec = 20
    bug = 0
    TMUlist = [[],
            [1.828,100],[1.814,200],[1.822,300],[1.819,400],[1.946,200],[1.728,200],
            [1.828,100],[1.853,200],[1.828,300],[1.827,400],[1.964,200],[1.698,200]]
    tmu = [0]*13
    for i in range(1,7):
        print('\nStarting simulation of TMU%i'%i)
        pt = TMUlist[i][1]
        qt = 0
        vt = TMUlist[i][0]
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        if i==5:
            wd,tmu[i] = TC.undrained_test(pt,qt,vt,3*eqd,prec,params,bug)
        else:
            wd,tmu[i] = TC.undrained_test(pt,qt,vt,eqd,prec,params,bug)
    for i in range(7,13):
        print('\nStarting simulation of TMU%i'%i)
        eqdm = -eqd
        if i == 11:
            eqdm *=7
        pt = TMUlist[i][1]
        qt = 0
        vt = TMUlist[i][0]
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wd,tmu[i] = TC.undrained_test(pt,qt,vt,eqdm,prec,params,bug)
    
    # plotU([TMU[i][0] for i in [1,2,3,4]],[],['k:','b--','g','r-.','m:'],[500,-0,500,0,0.05,-0.00,0.025,-0.1],['p'+'$_t$'+'=%i kPa'%i for i in [100,200,300,400]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='a) Undrained for varied p'+'$_t$'+': data')
    # plotU([TMU[i][0] for i in [5,2,6]],[],['k:','b--','g','r-.','m:'],[500,-0,500,0,0.05,-0.00,0.025,-0.1],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.946,1.814,1.728]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='b) Undrained for varied v'+'$_t$'+' (p'+'$_t$'+'=200 kPa): data')
    # plotU([TMU[i][0] for i in [11,8,12]],[],['k:','b--','g','r-.','m:'],[0,-500,500,0,0.0,-0.05,0.025,-0.1],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.946,1.814,1.728]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='c) Undrained for varied v'+'$_t$'+' (extension): data')
    # plotU([],[tmu[i][0][100:] for i in [1,2,3,4]],['k:','b--','g','r-.','m:'],[500,-0,500,0,0.05,-0.00,0.025,-0.1],['p'+'$_t$'+'=%i kPa'%i for i in [100,200,300,400]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='d) Undrained for varied p'+'$_t$'+': HySand_base')
    # plotU([],[tmu[i][0][100:] for i in [5,2,6]],['k:','b--','g','r-.','m:'],[500,-0,500,0,0.05,-0.00,0.025,-0.1],['p'+'$_t$'+'=%i kPa'%i for i in [100,200,300,400]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='e) Undrained for varied v'+'$_t$'+': HySand_base')
    # plotU([],[tmu[i][0][100:] for i in [11,8,12]],['k:','b--','g','r-.','m:'],[0,-500,500,0,0.0,-0.05,0.025,-0.1],['v'+'$_t$'+'={:.3f}'.format(i) for i in [1.946,1.814,1.728]],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=5,ncl=2,lw=1.5,ttl='f) Undrained for varied v'+'$_t$'+': HySand_base')

##### CYCLIC STRESS CONTROLLED ISOTROPIC #####
tcuit = 0
if tcuit == 1:
    eqcut = 0.025
    eqprec = 0.000004
    bug = 0
    TCUIlist = [[],
                [1.952,200,30,125,0], #1
                [1.962,100,15,200,1], #2
                [1.980,100,20,50 ,1], #3
                [1.963,100,25,25 ,1], #4
                [1.956,200,30,150,1], #5
                [1.964,300,45,150,1], #6
                [1.800,200,60,15 ,0], #7
                [1.822,100,20,250,1], #8
                [1.799,100,25,100,1], #9
                [1.826,100,30,20 ,1], #10
                [1.843,200,40,150,1], #11
                [1.814,200,50,100,1], #12
                [1.834,200,60,20 ,1], #13
                [1.847,300,60,250,1], #14
                [1.809,300,75,100,1], #15
                [1.817,300,90,50 ,1], #16
                # [1.726,200,60,200,0], #17
                [1.726,200,60,50,0], #17
                [1.760,100,30,50 ,1], #18
                [1.762,100,40,15 ,1], #19
                [1.757,100,50,8  ,1], #20
                [1.755,200,60,70 ,1], #21
                [1.745,300,90,250,1]] #22
    tcui = [0]*23
    TCUINh = [0]*23
    TCUINm1 = [0]*23
    TCUINm2 = [0]*23
    TCUINl = [0]*23
    for i in range(1,23):
        print('\nStarting simulation of TCUI%i'%i)
        pt = TCUIlist[i][1]
        qt = 0
        vt = TCUIlist[i][0]
        qc = TCUIlist[i][2] # amplitude of deviatoric stress
        ncyc = TCUIlist[i][3] # Number of cycles
        firstd = TCUIlist[i][4] # first drained cycle (1 yes, 0 no)
        pcut = pt/7.5
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wc,tcui[i],TCUINh[i],TCUINm1[i],TCUINm2[i],TCUINl[i] = TC.TCUS_test(pt,qt,vt,qc,ncyc,pcut,eqcut,prec,eqprec,params,bug,firstd)
# plotU(TCUI[1] ,[],['k'],[75,-75,225,0,0.01,-0.01,0.002,-0.1],['v'+'$_t$'+'=1.952'],lab=['',''],yt1=7,yt2=7,xt1=9,xt2=4,ncl=2,lw=1.2,ttl='a) Loose sand: data')
# plotU(TCUI[7] ,[],['k'],[75,-75,225,0,0.01,-0.01,0.002,-0.1],['v'+'$_t$'+'=1.8  '],lab=['',''],yt1=7,yt2=7,xt1=9,xt2=4,ncl=2,lw=1.2,ttl='b) Medium dense sand: data')
# plotU(TCUI[17],[],['k'],[75,-75,225,0,0.01,-0.01,0.002,-0.1],['v'+'$_t$'+'=1.726'],lab=['',''],yt1=7,yt2=7,xt1=9,xt2=4,ncl=2,lw=1.2,ttl='c) Dense sand: data')
# plotU([TCUI[7][0][200:320]] ,[],['k'],[75,-75,225,0,0.001,-0.0005,0.002,-0.1],['v'+'$_t$'+'=1.8  '],lab=['',''],yt1=7,yt2=7,xt1=9,xt2=3,ncl=2,lw=1.2,ttl='a) Medium dense sand: data')

# plotU([],[tcui[7][0][100:]],['b'],[75,-75,225,0,0.05,-0.05,0.002,-0.1],['v'+'$_t$'+'=1.8  '],lab=['',''],yt1=7,yt2=7,xt1=9,xt2=4,ncl=2,lw=1.2,ttl='', ori = 'h')

##### CYCLIC STRESS CONTROLLED ANISOTROPIC #####
tcuat = 0
if tcuat == 1:
    eqcut = 0.75
    bug = 0
    prec = 5 #25 # Changed from 10
    eqprec = 0.000005 # instead of 0.000001
    TCUAlist = [[],
                [1.945,300,150 ,120,10,0,1 ], #1
                [1.812,300,150 ,120,999,0,2 ], #2
                [1.727,300,150 ,120,999,0,3 ], #3
                [1.817,300,150 ,60 ,999,1,4 ], #4
                [1.820,300,150 ,90 ,999,1,5 ], #5
                [1.825,300,150 ,120,683,1,6 ], #6
                [1.916,300,150 ,120,100,1,7 ], #7
                [1.732,300,150 ,120,2237,1,8 ], #8
                [1.841,100,75  ,25 ,999,1,9 ], #9
                [1.806,200,150 ,50 ,999,1,10], #10
                [1.819,300,225 ,75 ,687,1,11], #11
                [1.827,200,200 ,60 ,999,1,12], #12 
                [1.814,200,150 ,60 ,1981,1,13], #13
                [1.840,200,100 ,60 ,999,1,14], #14
                [1.808,200,60  ,60 ,130,1,15], #15
                [1.844,200,50  ,60 ,100,1,16], #16
                [1.810,200,-50 ,60 ,100,1,17], #17
                [1.824,200,-50 ,40 ,100,1,18], #18
                [1.808,200,-100,40 ,1781,1,19]] #19

    # tcua = [0]*20
    # TCUANh = [0]*20
    # TCUANm1 = [0]*20
    # TCUANm2 = [0]*20
    # TCUANl = [0]*20
    TCUANavg = [0]*20
    tcuansimul = [[]]
    for i in range(1,20):
        print('\nStarting simulation of TCUA%i'%i)
        vt = TCUAlist[i][0]
        pt = TCUAlist[i][1]
        qt = TCUAlist[i][2]
        qc = TCUAlist[i][3]
        n  = TCUAlist[i][4]
        firstd = TCUAlist[i][5]
        pcut = pt/20
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wc1,tcua[i],TCUANh[i],TCUANm1[i],TCUANm2[i],TCUANl[i] = TC.TCUS_test(pt,qt,vt,qc,n,pcut,eqcut,prec,eqprec,params,bug,firstd,opt=2)
        TCUANavg[i] = [0]*len(TCUANh[i])
        eqini = tcua[i][0][100+int(firstd)*60][2]
        tcor = []
        for j in range(len(TCUANh[i])):
            TCUANavg[i][j] = [0]*3
            for k in range(2):
                TCUANavg[i][j][k] = (TCUANh[i][j][k]+TCUANl[i][j][k]+TCUANm1[i][j][k]+TCUANm2[i][j][k])/4
            TCUANavg[i][j][2] = (TCUANh[i][j][2]+TCUANl[i][j][2]+TCUANm1[i][j][2]+TCUANm2[i][j][2])/4-tcua[i][0][100][2]
            if i<17:
                tcor.append([j+1,TCUANm2[i][j][1],TCUANm2[i][j][2]-eqini])
            else:
                tcor.append([j+1,TCUANm2[i][j][1],-TCUANm2[i][j][2]-eqini])
        tcor.insert(0,[0,float(pt),0])
        tcuansimul.append([tcor])
    # plotU(TCUA[1],[],['k'],[300,-0,350,0,0.24,-0.,0.002,-0.1],['v'+'$_t$'+'=1.945'],lab=['',''],yt1=6,yt2=6,xt1=7,xt2=6,ncl=2,lw=1.2,ttl='a) Loose sand: data')
    # plotU(TCUA[2],[],['k'],[300,-0,350,0,0.12,-0.,0.002,-0.1],['v'+'$_t$'+'=1.812'],lab=['',''],yt1=6,yt2=6,xt1=7,xt2=3,ncl=2,lw=1.2,ttl='b) Medium dense sand: data')
    # plotU(TCUA[3],[],['k'],[300,-0,350,0,0.12,-0.,0.002,-0.1],['v'+'$_t$'+'=1.727'],lab=['',''],yt1=6,yt2=6,xt1=7,xt2=3,ncl=2,lw=1.2,ttl='c) Dense sand: data')
    # plotU([],[tcua[1][0][100:]],['b'],[300,-0,350,0,0.24,-0.,0.002,-0.1],[''],lab=['',''],yt1=6,yt2=6,xt1=7,xt2=6,ncl=2,lw=1.2,ttl='d) Loose sand: HySand_base')
    # plotU([],[tcua[2][0][100:]],['b'],[300,-0,350,0,0.12,-0.,0.002,-0.1],[''],lab=['',''],yt1=6,yt2=6,xt1=7,xt2=3,ncl=2,lw=1.2,ttl='e) Medium dense sand: HySand_base')
    # plotU([],[tcua[3][0][100:]],['b'],[300,-0,350,0,0.08,-0.,0.002,-0.1],[''],lab=['',''],yt1=6,yt2=6,xt1=7,xt2=2,ncl=2,lw=1.2,ttl='f) Dense sand: HySand_base')
# tcuant = 1
# if tcuant == 1:
#     tcuansimul = [[]]
#     for i in range(1,20):
#         pt = TCUAlist[i][1]
#         firstd = TCUAlist[i][5]
#         eqini = tcua[i][0][100+int(firstd)*60][2]
#         tcor = []
#         for j in range(len(TCUANm2[i])):
#             if i<17:
#                 tcor.append([j+1,TCUANm2[i][j][1],TCUANm2[i][j][2]-eqini])
#             else:
#                 tcor.append([j+1,TCUANm2[i][j][1],-TCUANm2[i][j][2]-eqini])
#         tcor.insert(0,[0,float(pt),0])
#         tcuansimul.append([tcor])
# plotTCUAN(TCUAN[2],tcuansimul[2],['k','b--'],limits=[1000,0,0.12,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=6,xt2=5,yt2=5,ttl="a) TCUA2")      
# plotTCUAN(TCUAN[3],tcuansimul[3],['k','b--'],limits=[1000,0,0.06,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=3,xt2=5,yt2=5,ttl="b) TCUA3")  
# plotTCUAN(TCUAN[4],tcuansimul[4],['k','b--'],limits=[1000,0,0.04,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=2,xt2=5,yt2=5,ttl="c) TCUA4")
# plotTCUAN(TCUAN[5],tcuansimul[5],['k','b--'],limits=[1000,0,0.08,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=4,xt2=5,yt2=5,ttl="d) TCUA5")
# plotTCUAN(TCUAN[6],tcuansimul[6],['k','b--'],limits=[750,0,0.12,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=3,yt1=6,xt2=5,yt2=5,ttl="e) TCUA6")
# plotTCUAN(TCUAN[7],tcuansimul[7],['k','b--'],limits=[250,0,0.16,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=1,yt1=8,xt2=5,yt2=5,ttl="f) TCUA7")
# plotTCUAN(TCUAN[8],tcuansimul[8],['k','b--'],limits=[2500,0,0.08,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=10,yt1=4,xt2=5,yt2=5,ttl="g) TCUA8")
# plotTCUAN(TCUAN[9],tcuansimul[9],['k','b--'],limits=[1000,0,0.04,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=2,xt2=5,yt2=5,ttl="h) TCUA9")
# plotTCUAN(TCUAN[10],tcuansimul[10],['k','b--'],limits=[1000,0,0.02,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=1,xt2=5,yt2=5,ttl="i) TCUA10")
# plotTCUAN(TCUAN[11],tcuansimul[11],['k','b--'],limits=[750,0,0.04,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=3,yt1=2,xt2=5,yt2=5,ttl="j) TCUA11")
# plotTCUAN(TCUAN[12],tcuansimul[12],['k','b--'],limits=[1000,0,0.04,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=2,xt2=5,yt2=5,ttl="k) TCUA12")
# plotTCUAN(TCUAN[13],tcuansimul[13],['k','b--'],limits=[2000,0,0.04,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=8,yt1=2,xt2=5,yt2=5,ttl="l) TCUA13")
# plotTCUAN(TCUAN[14],tcuansimul[14],['k','b--'],limits=[1000,0,0.08,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=4,xt2=5,yt2=5,ttl="m) TCUA14")
# plotTCUAN(TCUAN[15],tcuansimul[15],['k','b--'],limits=[250,0,0.06,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=1,yt1=3,xt2=5,yt2=5,ttl="n) TCUA15")
# plotTCUAN(TCUAN[16],tcuansimul[16],['k','b--'],limits=[250,0,0.08,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=1,yt1=3,xt2=5,yt2=5,ttl="o) TCUA16")
# plotTCUAN(TCUAN[17],tcuansimul[17],['k','b--'],limits=[250,-0.04,0,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=1,yt1=2,xt2=5,yt2=5,ttl="p) TCUA17")
# plotTCUAN(TCUAN[18],tcuansimul[18],['k','b--'],limits=[250,-0.02,0,1,0],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=1,yt1=1,xt2=5,yt2=5,ttl="q) TCUA18")

    
##### CYCLIC Strain CONTROLLED #####
tcuet = 0
if tcuet == 1:
    eqcut = 0.5
    bug = 0
    #TCUE list    vt   pt  qt    eqc  n firstd   CORRECTED STARTING V0 CORRECTED number of cycles
    TCUElist = [[],
                [1.825,200,150 ,0.0006,250,0], #1
                [1.802,200,150 ,0.0006,350,1], #2
                [1.863,200,-100,0.0006,150,0], #3
                [1.844,200,150 ,0.0004,1000,1], #4
                [1.807,200,150 ,0.0008,250,1], #5
                [1.939,200,150 ,0.0006,200,1], #6
                [1.757,200,150 ,0.0006,700,1], #7
                [1.805,200,230 ,0.0006,450,1], #8
                [1.804,200,0   ,0.0006,200,1], #9
                [1.854,200,-100,0.0006,300,1], #10
                [1.851,200,-150,0.0006,250,1], #11
                [1.847,50,37.5 ,0.0006,100,1], #12
                [1.800,100,75  ,0.0006,150,1], #13
                [1.820,300,225 ,0.0006,700,1], #14
                
                [1.944,200,0.0 ,0.01  ,10,0], #15    #Large strains
                [1.804,200,0.0 ,0.01  ,10,0], #16
                [1.698,200,0.0 ,0.01  ,20,0], #17
                [1.812,100,0.0 ,0.01  ,10,0], #18
                [1.814,700,0.0 ,0.01  ,10,0], #19
                [1.816,200,0.0 ,0.005 ,10,0], #20
                [1.841,200,150 ,0.01  ,10,0], #21
                [1.686,100,0.0 ,0.01  ,20,0], #22
                [1.674,700,0.0 ,0.01  ,20,0]] #23
    tcue = [0]*24
    tcuenp = [0]*24
    tcuensimul = [[]]
    for i in range(1,24):
        prec = 10
        if i>=15: prec*=10
        print('\nStarting simulation of TCUE%i'%i)
        vt = TCUElist[i][0]
        pt = TCUElist[i][1]
        qt = TCUElist[i][2]
        eqc = TCUElist[i][3]
        n = TCUElist[i][4]
        firstd = TCUElist[i][5]
        pcut = 1
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wc1,tcue[i],tcuenp[i] = TC.TCUE_test(pt,qt,vt,eqc,n,pcut,prec,params,bug,firstd)
        tcueni = []
        for j in range(len(TCUENp[i])):
            tcueni.append([j,TCUENp[i][j]])
        tcuensimul.append([tcueni])
# plotU(TCUE[9],[],['k'],[100,-100,220,0,0.0008,-0.0008,0.002,-0.1],['v'+'$_t$'+'=1.804'],lab=['',''],yt1=9,yt2=9,xt1=5,xt2=4,ncl=2,lw=1.2,ttl='a) TCUE9: data')
# plotU(TCUE[1],[],['k'],[175,-25,220,0,0.0008,-0.0008,0.002,-0.1],['v'+'$_t$'+'=1.825'],lab=['',''],yt1=9,yt2=9,xt1=5,xt2=4,ncl=2,lw=1.2,ttl='b) TCUE1: data')

# plotU([],[tcue[9][0][200:]],['b'],[100,-100,220,0,0.0008,-0.0008,0.002,-0.1],[''],lab=['',''],yt1=9,yt2=9,xt1=5,xt2=4,ncl=2,lw=1.2,ttl='d) TCUE9: simulation')
# plotU([],[tcue[1][0][100:]],['b'],[175,-25,220,0,0.0008,-0.0008,0.002,-0.1],[''],lab=['',''],yt1=9,yt2=9,xt1=5,xt2=4,ncl=2,lw=1.2,ttl='e) TCUE1: simulation')
# plotU([],[tcue[10][0][200:]],['b'],[50,-150,220,0,0.0008,-0.0008,0.002,-0.1],[''],lab=['',''],yt1=9,yt2=9,xt1=5,xt2=4,ncl=2,lw=1.2,ttl='f) TCUE10: simulation')

# plotU(TCUE[15],[],['k'],[200,-100,300,0,0.012,-0.012,0.002,-0.1],['v'+'$_t$'+'=1.944'],lab=['',''],yt1=6,yt2=6,xt1=3,xt2=5,ncl=2,lw=1.2,ttl='a) TCUE15: data')
# plotU(TCUE[16],[],['k'],[200,-100,300,0,0.012,-0.012,0.002,-0.1],['v'+'$_t$'+'=1.804'],lab=['',''],yt1=6,yt2=6,xt1=3,xt2=5,ncl=2,lw=1.2,ttl='b) TCUE16: data')
# plotU(TCUE[17],[],['k'],[400,-200,300,0,0.012,-0.012,0.002,-0.1],['v'+'$_t$'+'=1.698'],lab=['',''],yt1=12,yt2=12,xt1=3,xt2=5,ncl=2,lw=1.2,ttl='c) TCUE17: data')
# plotU([],[tcue[15][0][100:]],['b'],[200,-100,300,0,0.012,-0.012,0.002,-0.1],[''],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=5,ncl=2,lw=1.2,ttl='d) TCUE15: HySand_base')
# plotU([],[tcue[16][0][100:]],['b'],[200,-100,300,0,0.012,-0.012,0.002,-0.1],[''],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=5,ncl=2,lw=1.2,ttl='e) TCUE16: HySand_base')
# plotU([],[tcue[17][0][100:]],['b'],[400,-200,300,0,0.012,-0.012,0.002,-0.1],[''],lab=['',''],yt1=12,yt2=6,xt1=6,xt2=5,ncl=2,lw=1.2,ttl='f) TCUE17: HySand_base')

# plotU(TCUE[18],[tcue[18][0][100:]],['k','b--'],[100,-50,150,0,0.012,-0.012,0.002,-0.1],['data','HySand_base'],lab=['',''],yt1=3,yt2=3,xt1=3,xt2=5,ncl=2,lw=1.2,ttl='a) TCUE18')
# plotU(TCUE[19],[tcue[19][0][100:]],['k','b--'],[450,-200,800,0,0.012,-0.012,0.002,-0.1],['data','HySand_base'],lab=['',''],yt1=13,yt2=10,xt1=16,xt2=5,ncl=2,lw=1.2,ttl='b) TCUE19')
# plotU(TCUE[20],[tcue[20][0][100:]],['k','b--'],[200,-100,300,0,0.012,-0.012,0.002,-0.1],['data','HySand_base'],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=5,ncl=2,lw=1.2,ttl='c) TCUE20')
# plotU(TCUE[21],[tcue[21][0][100:]],['k','b--'],[250,-50,300,0,0.012,-0.012,0.002,-0.1],['data','HySand_base'],lab=['',''],yt1=6,yt2=6,xt1=6,xt2=5,ncl=2,lw=1.2,ttl='d) TCUE21')
# plotU(TCUE[22],[tcue[22][0][100:]],['k','b--'],[300,-200,200,0,0.012,-0.012,0.002,-0.1],['data','HySand_base'],lab=['',''],yt1=10,yt2=10,xt1=3,xt2=5,ncl=2,lw=1.2,ttl='e) TCUE22')
# plotU(TCUE[23],[tcue[23][0][100:]],['k','b--'],[1000,-600,800,0,0.012,-0.012,0.002,-0.1],['data','HySand_base'],lab=['',''],yt1=32,yt2=32,xt1=16,xt2=5,ncl=2,lw=1.2,ttl='f) TCUE23')
        
# plotTCUEN(TCUEN[1],tcuensimul[1],['k','b--'],limits=[250,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=4,ttl="a) TCUE1")
# plotTCUEN(TCUEN[2],tcuensimul[2],['k','b--'],limits=[300,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=6,yt1=4,ttl="b) TCUE2")
# plotTCUEN(TCUEN[3],tcuensimul[3],['k','b--'],limits=[150,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=3,yt1=4,ttl="c) TCUE3")
# plotTCUEN(TCUEN[4],tcuensimul[4],['k','b--'],limits=[1000,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=20,yt1=4,ttl="d) TCUE4")
# plotTCUEN(TCUEN[5],tcuensimul[5],['k','b--'],limits=[200,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=4,ttl="e) TCUE5")
# plotTCUEN(TCUEN[6],tcuensimul[6],['k','b--'],limits=[200,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=4,ttl="f) TCUE6")
# plotTCUEN(TCUEN[7],tcuensimul[7],['k','b--'],limits=[700,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=14,yt1=4,ttl="g) TCUE7")
# plotTCUEN(TCUEN[8],tcuensimul[8],['k','b--'],limits=[450,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=9,yt1=4,ttl="h) TCUE8")
# plotTCUEN(TCUEN[9],tcuensimul[9],['k','b--'],limits=[200,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=4,ttl="i) TCUE9")
# plotTCUEN(TCUEN[10],tcuensimul[10],['k','b--'],limits=[250,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=4,ttl="j) TCUE10")
# plotTCUEN(TCUEN[11],tcuensimul[11],['k','b--'],limits=[200,0,200],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=4,yt1=4,ttl="k) TCUE11")
# plotTCUEN(TCUEN[12],tcuensimul[12],['k','b--'],limits=[100,0,50],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=2,yt1=1,ttl="l) TCUE12")
# plotTCUEN(TCUEN[13],tcuensimul[13],['k','b--'],limits=[150,0,100],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=3,yt1=2,ttl="m) TCUE13")
# plotTCUEN(TCUEN[14],tcuensimul[14],['k','b--'],limits=[700,0,300],legend0=['data','simulation'],lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=14,yt1=6,ttl="n) TCUE14")

##### DRAINED WITH UNLOAD-RELOADS #####
tmcdt = 0
if tmcdt == 1:
    prec = 10
    bug = 0    
    TMCDlist = [[],
                [1.962,100,0.02,12,[0.0177,0.0373,0.05778,0.0788,0.100,0.1226,0.1457,0.1691,0.19323,0.2185,0.2436,0.2698]], #1
                [1.829,100,0.02,12,[0.01877,0.0403,0.062558,0.0855,0.1092,0.133,0.1575,0.1823,0.20795,0.2336,0.2601,0.2869]], #2
                [1.701,100,0.02,12,[0.0218,0.04696,0.07284,0.09873,0.1246,0.1505,0.1766,0.2025,0.2286,0.25496,0.28106,0.30745]], #3
                [1.820,100,0.06,4 ,[0.0641,0.13675,0.21284,0.292194]], #4
                [1.821,100,0.12,2 ,[0.1374,0.2934]], #5
                [1.810,200,0.02,12,[0.01877,0.04034,0.0627335,0.08596,0.1096,0.1337,0.1582,0.1831,0.2085,0.2346,0.2606,0.2867]], #6
                [1.814,50 ,0.02,12,[0.01958,0.0418,0.0646,0.0881,0.112,0.136,0.1610,0.1856,0.2112,0.237,0.2633,0.2891]]] #7
    tmcd = [0]*8
    for i in range(1,8):
        print('\nStarting simulation of TMCD%i'%i)
        vt  = TMCDlist[i][0]
        pt  = TMCDlist[i][1]
        qt = 0
        e1c = TMCDlist[i][2]
        n   = TMCDlist[i][3]
        eql = TMCDlist[i][4]
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wc1,tmcd[i] = TC.TMCD_test(pt,qt,vt,e1c,n,prec,params,eql,bug)
# plotD(TMCD[1],[tmcd[1][0][100:]],['k','b--'],[500,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=10,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='a) TMCD1, v'+'$_t$'+'=1.962')
# plotD(TMCD[2],[tmcd[2][0][100:]],['k','b--'],[500,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=10,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='b) TMCD2, v'+'$_t$'+'=1.829')
# plotD(TMCD[3],[tmcd[3][0][100:]],['k','b--'],[500,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=10,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='c) TMCD3, v'+'$_t$'+'=1.701')
# plotD(TMCD[4],[tmcd[4][0][100:]],['k','b--'],[500,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=10,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='d) TMCD4, v'+'$_t$'+'=1.820')
# plotD(TMCD[5],[tmcd[5][0][100:]],['k','b--'],[500,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=10,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='e) TMCD5, v'+'$_t$'+'=1.821')
# plotD(TMCD[6],[tmcd[6][0][100:]],['k','b--'],[750,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=15,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='f) TMCD6, v'+'$_t$'+'=1.810')
# plotD(TMCD[7],[tmcd[7][0][100:]],['k','b--'],[250,-0,220,0,0.3,-0,0.04,-0.12],['Data','Simulation'],lab=['',''],yt1=8,yt2=5,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='g) TMCD7, v'+'$_t$'+'=1.814',ori='h')

##### UNDRAINED CYCLIC STRAIN CONTROLLED + MONOTONIC#####
tmcut = 0
if tmcut == 1:
    tmcu = [0]*7
    prec = 5
    # TMCU1->4
    def tmuctstrain(pt,qt,vt,eqc,eqa,n): # n: number of cycles
        pcut = 1
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        H.const(params)
        H.start()
        TC.start_t(pt,qt,vt,params,prec)
        for i in range(n):
            print("cycle:",i+1)
            H.general_inc([0,0,0,0,1,0,0,1,0.0,eqc+eqa/4,1,25,prec])
            H.general_inc([0,0,0,0,1,0,0,1,0.0,-2*eqc+eqa/2,1,50,prec])
            H.general_inc([0,0,0,0,1,0,0,1,0.0,eqc+eqa/4,1,25,prec])
            lstate = H.returnstate()
            psig = lstate[3]
            if psig < pcut:
                break
        H.end()
    eqc = 0.0006
    pt = 200
    qt = 150
    n = 500
    tmuctstrainl = [ ['vt','eqa'],[1.827,0.00002],[1.825,0.00005],[1.834,0.0001],[1.816,0.0003]]
    for i in range(1,5):
        print('Calculations for TMCU%i'%i)
        vt = tmuctstrainl[i][0]
        eqa = tmuctstrainl[i][1]
        tmuctstrain(pt,qt,vt,eqc,eqa,n)
        tmcu[i] = H.returnrec()
    # TMCU5
    print('Calculations for TMC5')
    vt = 1.803
    pt = 200
    qt = 0
    qc = 50
    pcut = 60
    eqcut = 0.025
    eqprec = 0.000004
    n = 100
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    for i in range(n):
        print("cycle:",i+1)
        H.und_stress_from_strain([1,qc,eqprec,5])
        lstate = H.returnstate()
        epsq = abs(lstate[2])
        if epsq > eqcut:
            break
        H.und_stress_from_strain([1,0,-eqprec,5])
        lstate = H.returnstate()
        psig = lstate[3]
        if psig < pcut:
            break
        H.und_stress_from_strain([1,-qc,-eqprec,5])
        lstate = H.returnstate()
        epsq = abs(lstate[2])
        if epsq > eqcut:
            break
        H.und_stress_from_strain([1,0,eqprec,5])
        lstate = H.returnstate()
        psig = lstate[3]
        if psig < pcut:
            break
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.05,1,200,prec])
    H.end()
    tmcu[5]=H.returnrec()    
    # TMCU6
    print('Calculations for TMC6')
    vt = 1.684
    pt = 200
    qt = 0
    qc = 100
    pcut = 0.05*pt
    eqcut = 0.025
    eqprec = 0.000004
    n = 10
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([0,0,0,1,1,0,0,0,0.0,400,1,200,prec])
    H.general_inc([0,0,0,1,1,0,0,0,0.0,-350,1,200,prec])
    H.general_inc([0,0,0,1,1,0,0,0,0.0,200,1,200,prec])
    H.general_inc([0,0,0,1,1,0,0,0,0.0,-250,1,200,prec])
    for i in range(n):
        print("cycle:",i+1)
        H.und_stress_from_strain([1,qc,eqprec,5])
        lstate = H.returnstate()
        epsq = abs(lstate[2])
        if epsq > eqcut:
            break
        H.und_stress_from_strain([1,0,-eqprec,5])
        lstate = H.returnstate()
        psig = lstate[3]
        if psig < pcut:
            break
        H.und_stress_from_strain([1,-qc,-eqprec,5])
        lstate = H.returnstate()
        epsq = abs(lstate[2])
        if epsq > eqcut:
            break
        H.und_stress_from_strain([1,0,eqprec,5])
        lstate = H.returnstate()
        psig = lstate[3]
        if psig < pcut:
            break
    H.end()
    tmcu[6]=H.returnrec()

# plotU(TMCU[1],[],['k'],[200,-50,250,0,0.12,-0.005,0.04,-0.12],['v'+'$_t$'+'=1.827'],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='a) TMCU1: data')
# plotU(TMCU[2],[],['k'],[200,-50,250,0,0.12,-0.005,0.04,-0.12],['v'+'$_t$'+'=1.825'],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='b) TMCU2: data')
# plotU(TMCU[3],[],['k'],[200,-50,250,0,0.12,-0.005,0.04,-0.12],['v'+'$_t$'+'=1.834'],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='c) TMCU3: data')
# plotU([],[tmcu[1][0][100:]],['b'],[200,-50,250,0,0.12,-0.005,0.04,-0.12],[''],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='d) TMCU1: HySand_base')
# plotU([],[tmcu[2][0][100:]],['b'],[200,-50,250,0,0.12,-0.005,0.04,-0.12],[''],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='e) TMCU2: HySand_base')
# plotU([],[tmcu[3][0][100:40000]],['b'],[200,-50,250,0,0.12,-0.005,0.04,-0.12],[''],lab=['',''],yt1=5,yt2=5,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='f) TMCU3: HySand_base')

# plotU(TMCU[4],[],['k'],[1750,-250,1250,0,0.12,-0.005,0.04,-0.12],['v'+'$_t$'+'=1.816'],lab=['',''],yt1=8,yt2=8,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='a) TMCU4: data')
# plotU(TMCU[5],[],['k'],[500,-100,300,0,0.04,-0.01,0.05,-0.12],['v'+'$_t$'+'=1.803'],lab=['',''],yt1=12,yt2=12,xt1=6,xt2=6,ncl=2,lw=1.2,ttl='b) TMCU5: data')
# plotU(TMCU[6],[],['k'],[450,-150,300,0,0.02,-0.04,0.05,-0.12],['v'+'$_t$'+'=1.684'],lab=['',''],yt1=12,yt2=12,xt1=6,xt2=7,ncl=2,lw=1.2,ttl='c) TMCU6: data')
# plotU([],[tmcu[4][0][100:22000]],['b'],[1750,-250,1250,0,0.12,-0.005,0.04,-0.12],[''],lab=['',''],yt1=8,yt2=8,xt1=5,xt2=7,ncl=2,lw=1.2,ttl='d) TMCU4: HySand_base')

# plotU([],[tmcu[6][0][100:]],['b'],[450,-150,300,0,0.02,-0.04,0.05,-0.12],[''],lab=['',''],yt1=12,yt2=12,xt1=6,xt2=7,ncl=2,lw=1.2,ttl='f) TMCU6: HySand_base')


##### PRE-LOADED UNDRAINED #####

tpt = 0
if tpt == 1 :
    qt = 0
    prec = 10
    tp = [0]*12
    
    vt = 1.940
    pt = 300
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.3,1,600,prec])
    H.end()
    tp[10] = H.returnrec()
    
    v0 =1.939
    pc0=300
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([0,1,1,0,0,0,0,0,0.0,500,1,50,prec])
    H.general_inc([0,1,1,0,0,0,0,0,0.0,-500,1,50,prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.3,1,600,prec])
    H.end()
    tp[9] = H.returnrec()
    
    v0 =1.938
    pc0=300
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    H.general_inc([3,-1,0,1,0,0,0,0,0.0,0.01,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.002,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.003,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.005,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.03,1,150,prec])
    lstate = H.returnstate()
    qsig = lstate[4]
    H.general_inc([3,-1,0,1,0,0,0,0,0.0,-qsig,1,100,prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.05,1,50,10*prec])
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.25,1,550,prec])
    H.end()
    tp[8] = H.returnrec()
    
# plotU(TP[10]+TP[9]+TP[8],[],['k','b--','r:'],[1600,-0,1000,0,0.3,-0.0,0.05,-0.12],['Normal\nconsolidation','Over-\nconsolidation\nto 800 kPa','Drained shear\npreloading'],lab=['',''],yt1=8,yt2=8,xt1=5,xt2=6,ncl=2,lw=1.2,ttl='a) Pre-loading influence: data',ori='h')
# plotU([],[tp[10][0][100:]]+[tp[9][0][100:]]+[tp[8][0][100:]],['k','b--','r:'],[1600,-0,1000,0,0.3,-0.0,0.05,-0.12],['Normal\nconsolidation','Over-\nconsolidation\nto 800 kPa','Drained shear\npreloading'],lab=['',''],yt1=8,yt2=8,xt1=5,xt2=6,ncl=2,lw=1.2,ttl='b) Pre-loading influence: HySand_base',ori='h')



##### DUQUE 7 COMMON LIMITATIONS #####


##### Over/Under shooting #####

ovrshtt = 0
if ovrshtt == 1 :
    prec = 10
    pt = 100
    qt = 0
    vt = Be - (Be-De) * 0.16
    oversht = [0]*2
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    lstate = H.returnstate()
    e10 = (lstate[1]+3*lstate[2]) / 3
    ep0 = lstate[1]
    H.general_inc([3,-1,0,1,0,0,0,0,0.0,80,1,100,prec])
    lstate = H.returnstate()
    e1 = (lstate[1]+3*lstate[2]) / 3 - e10
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.2-e1,1,500,prec])
    H.end()
    oversht[0] = H.returnrec()
    
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    lstate = H.returnstate()
    e10 = (lstate[1]+3*lstate[2]) / 3
    H.general_inc([3,-1,0,1,0,0,0,0,0.0,80,1,100,prec])
    for i in range(5):
        H.general_inc([0,0,0,0,1,0,0,1,0.0,0.003,1,50,prec])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,-0.0003,1,10,prec])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,0.0003,1,10,prec])
    e1 = (lstate[1]+3*lstate[2]) / 3 - e10
    H.general_inc([0,0,0,0,1,0,0,1,0.0,0.2-e1,1,500,prec])
    H.end()
    oversht[1] = H.returnrec()
# plotU([],[oversht[0][0][100:]]+[oversht[1][0][100:]],['m','b--'],[200,-0,200,0,0.2,-0.0,0.05,-0.12],['Monotonic','Cyclic'],lab=['',''],yt1=4,yt2=4,xt1=4,xt2=4,ncl=2,lw=1.2,ttl='',ori='h')

##### Cyclic strength curves #####
csct = 0
if csct == 1:
    CSR      = [0,   0.075,0.075,0.1,0.125,0.075,0.075,   0.15,0.1,0.125,0.15,0.1,0.125,0.15,0.1,0.125,0.15,      0.15,0.15,0.2,0.25,0.15,0.15]
    Nliqdata = [0,     72 , 87  , 15,  5  ,  86 , 144 ,    11 ,249, 100 , 15 ,146,  77 , 15 ,257, 110 , 24 ,       185, 54 , 15,  6 , 61 , 269]
    Nliqsimu = [0,    106 , 119 , 29, 13  , 141 , 137 ,    10 ,106,  49 , 17 ,106,  50 , 18 ,113,  57 , 22 ,        25, 28 ,  8,  4 , 35 ,  42]
    fig, ax = plt.subplots(figsize=(8, 4))
    lw0,mw0 = 1.2,5
    ax.plot(  Nliqdata[2:5],  CSR[2:5], lw=lw0,linestyle='-' ,color = 'c'   , marker='s',markersize=mw0, label='Loose, p'+'$_t$'+'=100kPa: data') 
    ax.plot(  Nliqsimu[2:5],  CSR[2:5], lw=lw0,linestyle='-' ,color = 'b'   , marker='s',markersize=mw0, label='Loose, p'+'$_t$'+'=100kPa: HySand_base') 
    ax.plot( Nliqdata[8:11], CSR[8:11], lw=lw0,linestyle='-' ,color = 'pink', marker='s',markersize=mw0, label='Medium, p'+'$_t$'+'=100kPa: data',fillstyle='none') 
    ax.plot( Nliqsimu[8:11], CSR[8:11], lw=lw0,linestyle='-' ,color = 'r'   , marker='s',markersize=mw0, label='Medium, p'+'$_t$'+'=100kPa: HySand_base',fillstyle='none') 
    ax.plot(Nliqdata[11:14],CSR[11:14], lw=lw0,linestyle='--',color = 'pink', marker='o',markersize=mw0, label='Medium, p'+'$_t$'+'=200kPa: data',fillstyle='none') 
    ax.plot(Nliqsimu[11:14],CSR[11:14], lw=lw0,linestyle='--',color = 'r'   , marker='o',markersize=mw0, label='Medium, p'+'$_t$'+'=200kPa: HySand_base',fillstyle='none') 
    ax.plot(Nliqdata[14:17],CSR[14:17], lw=lw0,linestyle=':' ,color = 'pink', marker='D',markersize=mw0, label='Medium, p'+'$_t$'+'=300kPa: data',fillstyle='none') 
    ax.plot(Nliqsimu[14:17],CSR[14:17], lw=lw0,linestyle=':' ,color = 'r'   , marker='D',markersize=mw0, label='Medium, p'+'$_t$'+'=300kPa: HySand_base',fillstyle='none')
    ax.plot(Nliqdata[18:21],CSR[18:21], lw=lw0,linestyle='-' ,color = 'lime', marker='s',markersize=mw0, label='Dense, p'+'$_t$'+'=100kPa: data') 
    ax.plot(Nliqsimu[18:21],CSR[18:21], lw=lw0,linestyle='-' ,color = 'g'   , marker='s',markersize=mw0, label='Dense, p'+'$_t$'+'=100kPa: HySand_base') 
    # ax.plot(Nliqdata[7:17], Nliqsimu[7:17], lw=0, color = '0.5', marker='p', label='Medium') 
    # ax.plot(Nliqdata[17:23], Nliqsimu[17:23], lw=0, color = '0.25', marker='D', label='Dense') 
    # ax.plot([1,500],[1,500],lw=1,color = '0.5', linestyle = ':')

    plt.xlabel('N'+'$_{cyc}$'+' to liquefaction', labelpad=10, fontsize=12)
    plt.ylabel('CSR', labelpad=10, fontsize=12)
    fig.legend(ncol=1, loc='center', prop={'size': 10}, bbox_to_anchor=(0.825, 0.5),facecolor='whitesmoke')
    # plt.legend(ncl=ncl0, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.5, 0.975-(nrow-1)*0.025),facecolor='whitesmoke')

    ax.set_xscale('log')
    plt.xlim([3,300])
    plt.ylim([0,0.25])
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_ticks([0,0.05,0.1,0.15,0.2,0.25])
    ax.xaxis.set_ticks([3,10,100,300])
    plt.grid(linestyle=':',linewidth=1)
    plt.subplots_adjust(bottom=0.14,top=0.95,left=0.1,right=0.65,wspace=0.25, hspace=0.3)
    ax.tick_params(direction="in")
    # plt.grid(linestyle=':',linewidth=1)
    plt.minorticks_on()
    plt.grid(which='minor',linestyle=':',linewidth=0.2)
    plt.show()  




##### Stress loops by Poblete #####
pobt = 0
if pobt == 1 :
    prec = 10
    vt = 1.857
    pt = 200
    qt = 100+20/(2/3)**0.5
    params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
    H.const(params)
    H.start()
    TC.start_t(pt,qt,vt,params,prec)
    nn=30
    for i in range(200*nn):
        print('Progress: %i/30'%int(i/200))
        dp = 20/3**0.5 * (np.pi/100) * np.cos(np.pi/100*i)
        dq = -20*(3/2)**0.5 * (np.pi/100) * np.sin(np.pi/100*i)
        H.general_inc([1,0,0,1,0,0,0,0,dp,dq,1,1,prec])
    H.end()
    pob = H.returnrec()
# plotpob([],[pob[0][100:]],test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncol=5)


###############################################################################
###############################################################################
###############################################################################
###############################   Calibration   ###############################
###############################################################################
###############################################################################
###############################################################################

# plotISOloglog([ISO[1][0][:185]],[],test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncol=5)
# plotISOpow([ISO[1][0][:185]],[],test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncol=5)