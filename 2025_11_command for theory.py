# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 11:09:31 2025

@author: lsimonin
"""

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
pre = 0
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


# For theory article, comparison of number of surfaces
compN=0
if compN == 1:
    H.model(["2025_02_HySand"])
    prec = 10
    [NN,vt ,  K  ,  G  ,ind ,phic, bet ,  Be , Ga ,  De ,  lb  ,  ld  , A0 ,  r  ,  h ,b,a0] = \
    [10,1.8,50000,400000,0.7,33,1,2.05,1.95,1.7,0.007,0.003,40,0.3,1000,3,0 ]
    eqd = 0.3
    i=0
    tmdb = [0]*5
    for NN in [5,10,20,50,100]:
        print('\nStarting simulation of TMD%i'%i)
        qt = 0
        vt,pt = 1.8,100 # Initial stress values
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wd,tmdb[i] = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
        i+=1
    i = 0
    tmub = [0]*5
    equ = 0.05
    for NN in [5,10,20,50,100]:
        print('\nStarting simulation of TMU%i'%i)
        qt = 0
        vt,pt = 1.8,100 # Initial stress values
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wd,tmub[i] = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
        i+=1

# plotD([],[tmdb[i][0][100:] for i in [0,1,2,3,4]],['k:','b--','g','r-.','m'],[400,-0,225,0,0.3,-0.0,0.01,-0.06],['N=%i'%i for i in [5,10,20,50,100]],lab=['',''],yt1=7,yt2=4,xt1=5,xt2=5,ncl=3,lw=1.2,ttl='a) Drained test',save=1)

# For theory article, comparison of number of surfaces, with different  combinations of r and m
compN2 = False
if compN2:
    H.model(["2025_02_HySand"])
    prec = 10
    [NN,vt ,  K  ,  G  ,ind ,phic, bet ,  Be , Ga ,  De ,  lb  ,  ld  , A0 ,  r  ,  h ,b,a0] = \
    [10,1.8,50000,400000,0.99,33   , 1  ,2.05 ,1.95, 1.7 ,0.007 ,0.003 , 40 , 1 ,1000,3,0 ]
    eqd = 0.3
    i=0
    tmdb = [0]*5
    for NN in [5,10,20,50]:
        print('\nStarting simulation of TMD%i'%i)
        qt = 0
        vt,pt = 1.8,100 # Initial stress values
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wd,tmdb[i] = TC.drained_test(pt,qt,vt,eqd,prec,params,bug)
        i+=1
    i = 0
    tmub = [0]*5
    equ = 0.05
    for NN in [5,10,20,50]:
        print('\nStarting simulation of TMU%i'%i)
        qt = 0
        vt,pt = 1.8,100 # Initial stress values
        params = [NN,vt,K,G,ind,phic,bet,Be,Ga,De,l0,l1,A0,r,h,b,a0]
        wd,tmub[i] = TC.undrained_test(pt,qt,vt,equ,prec,params,bug)
        i+=1

# plotD([],[tmdb[i][0][100:] for i in range(4)],['k:','b--','g','r-.','m'],[400,-0,225,0,0.3,-0.0,0.01,-0.06],['N=%i'%i for i in [5,10,20,50]],lab=['',''],yt1=7,yt2=4,xt1=5,xt2=5,ncl=3,lw=1.2,ttl='m=0.5 and r=0.5',save=1)



