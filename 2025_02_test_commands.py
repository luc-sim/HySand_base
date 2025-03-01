# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:41:42 2025

@author: lsimonin

Modified to have the initiation routine

"""

import numpy as np
from HyperDrive import Commands as H

def vi_calc(pt,qt,vt,params,prec,pr=0):
    vi = vt
    if pr==1:
        print("Initial specific volume:",vi)
    params[1] = vi
    H.const(params)
    H.start()
    p0=10**-6
    H.init_stress([p0,0])
    lstate = H.returnstate()
    eps0 = lstate[1]
    if pr==1:
        print('eps0:',eps0)
    H.general_inc([0,1,1,0,0,0,0,0,qt,pt-p0,1,100,prec])
    lstate = H.returnstate()
    eps_p = lstate[1]
    vt_c = params[1]*np.exp(-(eps_p-eps0))
    diff_v = vt_c - vt
    if pr==1:
        print("\np=",lstate[3])
        print("q="  ,lstate[4])
        print("After start (substracting initial value), eps_p=",eps_p-eps0)
        print("Such that at (pt,qt), vt_c=",vt_c)
        print("\n")
    vi -= diff_v
    i=0
    while abs(diff_v)>10**-6:
        if pr==1:
            print("Initial specific volume adjusted:",vi)
            print("\n Step number",i+1)
        params[1] = vi
        H.const(params)
        H.start()
        p0=10**-6
        H.init_stress([p0,0])
        lstate = H.returnstate()
        eps0 = lstate[1]
        if pr==1:
            print('eps0:',eps0)
        H.general_inc([0,1,1,0,0,0,0,0,qt,pt-p0,1,100,prec])
        lstate = H.returnstate()
        eps_p = lstate[1]
        vt_c = params[1]*np.exp(-(eps_p-eps0))
        diff_v = vt_c - vt
        if pr==1:
            print("\np=",lstate[3])
            print("q="  ,lstate[4])
            print("After start (substracting initial value), eps_p=",eps_p-eps0)
            print("Such that at (pt,qt), vt_c=",vt_c)
            print("\n")
        vi -= diff_v
        i+=1
    H.end()
    return vi

def start_t(pt,qt,vt,params,prec):
    H.const(params)
    H.start()
    # lstate = H.returnstate()
    # print("BEFORE START lstate=",lstate)
    print("\nITERATIVE vi CALCULATION STARTING:\n\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/")
    vi = vi_calc(pt,qt,vt,params,prec,pr=0)
    # lstate = H.returnstate()
    # print("p=",lstate[3])
    # print("q=",lstate[4])
    # print("AFTER vi lstate=",lstate)
    print("/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\nvi=",vi)
    print("\n")
    p0=10**-6
    H.start()
    H.init_stress([p0,0])
    lstate = H.returnstate()
    # print("p=",lstate[3])
    # print("q=",lstate[4])
    # print("After init stress lstate=",lstate)
    # print("\nvi=",vi)
    eps0 = lstate[1]
    H.general_inc([0,1,1,0,0,0,0,0,qt,pt-p0,1,100,prec])
    lstate = H.returnstate()
    print("p=",lstate[3])
    print("q=",lstate[4])
    print("\nCHECK OF INITIAL SPECIFIC VOLUME: vt_c-vt=",vi*np.exp(-(lstate[1]-eps0))-vt)

def drained_test(pt,qt,vt,eq,prec,params,bug):
    H.const(params)
    H.start()
    start_t(pt,qt,vt,params,prec)
    H.general_inc([3,-1,0,1,0,0,0,0,0.0,0.01,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.002,1,10,prec])
    itern = H.returnitern()
    if bug == 1:
        if itern<21: #if bug
            return (0,[])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.003,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.005,1,10,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.02,1,20,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.04,1,40,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.06,1,60,prec])
    H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.08,1,80,prec])
    itern = H.returnitern()
    if bug == 1:
        if itern<241: #if bug
            return (0,[])
    if eq>0.2:
        H.general_inc([3,-1,0,0,0,0,0,1,0.0,eq-0.2,1,100,prec])
    H.end()
    recmd=H.returnrec()
    w = 1
    return(w,recmd)

def undrained_test(pt,qt,vt,eq,prec,params,bug):
    H.const(params)
    H.start()
    start_t(pt,qt,vt,params,prec)
    if eq<0.0:
        H.general_inc([0,0,0,0,1,0,0,1,0.0,eq,1,120,prec])
    else:
        H.general_inc([0,0,0,0,1,0,0,1,0.0,0.002,1,10,prec])
        itern = H.returnitern()
        if bug == 1:
            if itern<11: #if bug
                return (0,[])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,0.03,1,100,prec])
        itern = H.returnitern()
        if bug == 1:
            if itern<111: #if bug
                return (0,[])
        if eq>0.03:
            H.general_inc([0,0,0,0,1,0,0,1,0.0,eq-0.03,1,100,prec])
    H.end()
    recmu=H.returnrec()
    w = 1
    return(w,recmu)

def undrained_test_ext(pt,qt,vt,eq,prec,params,bug): #in extension
    H.const(params)
    H.start()
    start_t(pt,qt,vt,params,prec)
    if eq<0.0:
        H.general_inc([0,0,0,0,1,0,0,1,0.0,-0.002,1,10,prec])
        itern = H.returnitern()
        if bug == 1:
            if itern<11: #if bug
                return (0,[])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,-0.03,1,100,prec])
        itern = H.returnitern()
        if bug == 1:
            if itern<111: #if bug
                return (0,[])
        if eq<-0.03:
            H.general_inc([0,0,0,0,1,0,0,1,0.0,eq+0.03,1,100,prec])
    else:
        H.general_inc([0,0,0,0,1,0,0,1,0.0,0.002,1,10,prec])
        itern = H.returnitern()
        if bug == 1:
            if itern<11: #if bug
                return (0,[])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,0.03,1,100,prec])
        itern = H.returnitern()
        if bug == 1:
            if itern<111: #if bug
                return (0,[])
        if eq>0.03:
            H.general_inc([0,0,0,0,1,0,0,1,0.0,eq-0.03,1,100,prec])
    H.end()
    recmu=H.returnrec()
    w = 1
    return(w,recmu)

def TCUS_test(pt,qt,vt,qc,n,pcut,eqcut,prec,eqprec,params,bug,firstd,opt=0):
    # opt=1 is designed for TCUA to have just the top kpas being strain controlled
    nsub = 5
    H.const(params)
    H.start()
    start_t(pt,qt,vt,params,prec)
    if firstd == 1:
        H.general_inc([3,-1,0,1,0,0,0,0,0.0,1.0,1,2,prec])
        H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc-1,1,18,prec])
        H.general_inc([3,-1,0,1,0,0,0,0,0.0,-2*qc,1,20,prec])
        H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc,1,20,prec])
    testinfoh = []
    testinfom1 = []
    testinfom2 = []
    testinfol = []
    if opt == 0:
        for i in range(n):
            print("cycle:",i+1)
            if i==0:
                H.und_stress_from_strain([1,qt+qc,eqprec,nsub*4])
            else:
                H.und_stress_from_strain([1,qt+qc,eqprec,nsub])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfoh.append([i,psig,epsq])
            if epsq > eqcut:
                break
            H.und_stress_from_strain([1,qt,-eqprec,nsub])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfom1.append([i,psig,epsq])
            if psig < pcut:
                break
            H.und_stress_from_strain([1,qt-qc,-eqprec,nsub])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfol.append([i,psig,epsq])
            if epsq > eqcut:
                break
            H.und_stress_from_strain([1,qt,eqprec,nsub])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfom2.append([i,psig,epsq])
            if psig < pcut:
                break
    if opt == 1:
        for i in range(n):
            print("cycle:",i+1)
            if i==0:
                H.general_inc([0,0,0,1,1,0,0,0,0.0,0.95*qc,1,20,prec])
                H.und_stress_from_strain([1,qt+qc,eqprec,nsub*4])
            else:
                H.general_inc([0,0,0,1,1,0,0,0,0.0,0.95*qc,1,20,prec])
                H.und_stress_from_strain([1,qt+qc,eqprec,nsub])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfoh.append([i,psig,epsq])
            if epsq > eqcut:
                break
            H.und_stress_from_strain([1,qt+0.95*qc,-eqprec,nsub])
            H.general_inc([0,0,0,1,1,0,0,0,0.0,-0.95*qc,1,20,prec])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfom1.append([i,psig,epsq])
            if psig < pcut:
                break
            H.general_inc([0,0,0,1,1,0,0,0,0.0,-0.95*qc,1,20,prec])
            H.und_stress_from_strain([1,qt-qc,-eqprec,nsub])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfol.append([i,psig,epsq])
            if epsq > eqcut:
                break
            H.und_stress_from_strain([1,qt-0.95*qc,eqprec,nsub])
            H.general_inc([0,0,0,1,1,0,0,0,0.0,0.95*qc,1,20,prec])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfom2.append([i,psig,epsq])
            if psig < pcut:
                break
    if opt == 2:
        print('IIIIIIIIIIIIIIIIIIIIIIIINNNNNNNNNNNNNNN')
        for i in range(n):
            print("cycle:",i+1)
            if i==0:
                H.und_stress_from_strain([1,qt+qc,eqprec,nsub*4])
            else:
                H.general_inc([0,0,0,1,1,0,0,0,0.0,qc,1,20,prec])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfoh.append([i,psig,epsq])
            if epsq > eqcut:
                break
            H.general_inc([0,0,0,1,1,0,0,0,0.0,-qc,1,20,prec])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfom1.append([i,psig,epsq])
            if psig < pcut:
                break
            H.general_inc([0,0,0,1,1,0,0,0,0.0,-qc,1,20,prec])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfol.append([i,psig,epsq])
            if epsq > eqcut:
                break
            H.general_inc([0,0,0,1,1,0,0,0,0.0,qc,1,20,prec])
            lstate = H.returnstate()
            epsq,psig = lstate[2],lstate[3]
            testinfom2.append([i,psig,epsq])
            if psig < pcut:
                break
    H.end()
    reccu=H.returnrec()
    w = 1
    return(w,reccu,testinfoh,testinfom1,testinfom2,testinfol)

def TCUE_test(pt,qt,vt,eqc,n,pcut,prec,params,bug,firstd):
    nstp = 25
    H.const(params)
    H.start()
    start_t(pt,qt,vt,params,prec)
    if firstd == 1:
        H.general_inc([3,-1,0,0,0,0,0,1,0.0,eqc,1,nstp,prec])
        H.general_inc([3,-1,0,0,0,0,0,1,0.0,-2*eqc,1,2*nstp,prec])
        lstate = H.returnstate()
        qactual = lstate[4]
        H.general_inc([3,-1,0,1,0,0,0,0,0.0,qt-qactual,1,nstp,prec])
    paccu = []           
    lstate = H.returnstate()
    paccu.append(lstate[3])
    for i in range(n):
        print("cycle:",i+1)
        H.general_inc([0,0,0,0,1,0,0,1,0.0,eqc,1,nstp,prec])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,-eqc,1,nstp,prec])
        lstate = H.returnstate()
        paccu.append(lstate[3])
        if lstate[3] < pcut: break
        H.general_inc([0,0,0,0,1,0,0,1,0.0,-eqc,1,nstp,prec])
        H.general_inc([0,0,0,0,1,0,0,1,0.0,eqc,1,nstp,prec])
        lstate = H.returnstate()
        if lstate[3] < pcut: break
    H.end()
    reccu=H.returnrec()
    w = 1
    return(w,reccu,paccu)


def TMCD_test(pt,qt,vt,e1c,n,prec,params,eqlist,bug): #ini pc #cyclic e1 #cycles number #how precise #bug
    H.const(params)
    H.start()
    start_t(pt,qt,vt,params,prec)
    H.general_inc([3,-1,0,1,0,0,0,0,0.0,0.01,1,10,4])
    for i in range(len(eqlist)):
        print('cycle:',i+1)
        lstate = H.returnstate()
        eqstart = lstate[2]
        eqaim = eqlist[i]-eqstart
        H.general_inc([3,-1,0,0,0,0,0,1,0.0,eqaim,1,100,5*prec])
        lstate = H.returnstate()
        qsig = lstate[4]
        while qsig>0:
            H.general_inc([3,-1,0,0,0,0,0,1,0.0,-0.00001,1,1,5*prec])
            lstate = H.returnstate()
            qsig = lstate[4]
    H.end()
    tmcd=H.returnrec()
    w = 1
    return(w,tmcd)

# def TMCD_test2(pt,qt,vt,e1c,n,prec,params,eqlist,bug): #ini pc #cyclic e1 #cycles number #how precise #bug
#     H.const(params)
#     H.start()
#     start_t(pt,qt,vt,params,prec)
#     H.general_inc([3,-1,0,1,0,0,0,0,0.0,0.01,1,10,4])
#     for i in range(len(eqlist)):
#         print('cycle:',i+1)
#         lstate = H.returnstate()
#         eqstart = lstate[2]
#         H.general_inc([3,-1,0,0,0,0,1/3,1,0.0,(i+1)*e1c-eqstart,1,100,5*prec])
#         lstate = H.returnstate()
#         qsig = lstate[4]
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,-qsig,1,100,5*prec])
#     H.end()
#     tmcd=H.returnrec()
#     w = 1
#     return(w,tmcd)

# def TCUI_test(pt,qt,vt,qc,n,pcut,eqcut,prec,eqprec,params,bug,firstd):
#     nsub = 5
#     H.const(params)
#     H.start()
#     start_t(pt,qt,vt,params,prec)
#     if firstd == 1:
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,1.0,1,2,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc-1,1,18,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,-2*qc,1,20,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc,1,20,prec])
#     H.general_inc([0,0,0,1,1,0,0,0,0.0,qc,1,20,prec])
#     itern = H.returnitern()
#     if bug == 1:
#         if itern<81: #if bug
#             return (0,[])
#     for i in range(n):
#         print("cycle:",i+1)
#         H.und_stress_from_strain([1,0,-eqprec,nsub])
#         lstate = H.returnstate()
#         psig = lstate[3]
#         if psig < pcut:
#             break
#         H.und_stress_from_strain([1,-qc,-eqprec,nsub])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         if epsq > eqcut:
#             break
#         H.und_stress_from_strain([1,0,eqprec,nsub])
#         lstate = H.returnstate()
#         psig = lstate[3]
#         if psig < pcut:
#             break
#         H.und_stress_from_strain([1,qc,eqprec,nsub])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         if epsq > eqcut:
#             break
#     H.end()
#     reccu=H.returnrec()
#     w = 1
#     return(w,reccu)

# def TCUA_test(p0,q0,qc,n,pcut,eqcut,prec,eqprec,params,bug,firstd):
#     params[2] = p0-q0+0.1
#     H.const(params)
#     H.start()
#     if q0>p0:
#         print('q0>p0, need to code this')
#         return(0,[])
#     H.init_stress([p0-np.abs(q0),0])
#     H.general_inc([1,-1*np.sign(q0),0,1,0,0,0,0,0.0,q0,1,50,prec])
#     if firstd == 1:
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,1.0,1,2,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc-1,1,18,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,-2*qc,1,20,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc,1,20,prec])
#     # H.general_inc([0,0,0,1,1,0,0,0,0.0,qc,1,20,prec])
#     # itern = H.returnitern()
#     # if bug == 1:
#     #     if itern<81: #if bug
#     #         return (0,[])
#     for i in range(n):
#         print("cycle:",i+1)
#         H.und_stress_from_strain([1,q0+qc,eqprec,5])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         if epsq > eqcut:
#             break
#         H.und_stress_from_strain([1,q0,-eqprec,5])
#         lstate = H.returnstate()
#         psig = lstate[3]
#         if psig < pcut:
#             break
#         H.und_stress_from_strain([1,q0-qc,-eqprec,5])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         if epsq > eqcut:
#             break
#         H.und_stress_from_strain([1,q0,eqprec,5])
#         lstate = H.returnstate()
#         psig = lstate[3]
#         if psig < pcut:
#             break
        
#     H.end()
#     reccu=H.returnrec()
#     w = 1
#     return(w,reccu)

# def TCUA_test_new(p0,q0,qc,n,pcut,eqcut,prec,eqprec,params,bug,firstd): # 2024 addition
#     params[2] = 5
#     H.const(params)
#     H.start()
#     if q0>p0:
#         print('q0>p0, need to code this')
#         return(0,[])
#     H.init_stress([0.1,0])
#     H.general_inc([abs(q0)/p0,-1*np.sign(q0),0,1,0,0,0,0,0.0,q0,1,50,prec])
#     if firstd == 1:
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,1.0,1,2,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc-1,1,18,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,-2*qc,1,20,prec])
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,qc,1,20,prec])
#     # H.general_inc([0,0,0,1,1,0,0,0,0.0,qc,1,20,prec])
#     # itern = H.returnitern()
#     # if bug == 1:
#     #     if itern<81: #if bug
#     #         return (0,[])
#     testinfoh = []
#     testinfom1 = []
#     testinfom2 = []
#     testinfol = []
#     for i in range(n):
#         print("cycle:",i+1)
#         if i==0:
#             H.und_stress_from_strain([1,q0+qc,eqprec,20])
#         else:
#             H.und_stress_from_strain([1,q0+qc,eqprec,5])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         psig = abs(lstate[3])
#         testinfoh.append([i,psig,epsq])
#         if epsq > eqcut:
#             break
#         H.und_stress_from_strain([1,q0,-eqprec,5])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         psig = abs(lstate[3])
#         testinfom1.append([i,psig,epsq])
#         if psig < pcut:
#             break
#         H.und_stress_from_strain([1,q0-qc,-eqprec,5])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         psig = abs(lstate[3])
#         testinfol.append([i,psig,epsq])
#         if epsq > eqcut:
#             break
#         H.und_stress_from_strain([1,q0,eqprec,5])
#         lstate = H.returnstate()
#         epsq = abs(lstate[2])
#         psig = abs(lstate[3])
#         testinfom2.append([i,psig,epsq])
#         if psig < pcut:
#             break
        
#     H.end()
#     reccu=H.returnrec()
#     w = 1
#     return(w,reccu,testinfoh,testinfom1,testinfom2,testinfol)

# def TCUE_test(p0,q0,eqc,n,pcut,prec,params,bug,firstd):
    
#     if abs(q0)>p0:
#         params[2] = p0-np.abs(q0)/2+0.1
#         H.const(params)
#         H.start()
#         H.init_stress([p0-np.abs(q0)/2,0])
#         H.general_inc([2,-1,0,1,0,0,0,0,0.0,q0,1,50,prec])
#     else:
#         params[2] = p0-np.abs(q0)+0.1
#         H.const(params)
#         H.start()
#         H.init_stress([p0-np.abs(q0),0])
#         if q0>0:
#             H.general_inc([1,-1,0,1,0,0,0,0,0.0,q0,1,50,prec])
#         if q0<0:
#             H.general_inc([-1,-1,0,1,0,0,0,0,0.0,q0,1,50,prec])
#     if firstd == 1:
#         H.general_inc([3,-1,0,0,0,0,0,1,0.0,eqc,1,25,prec])
#         H.general_inc([3,-1,0,0,0,0,0,1,0.0,-2*eqc,1,25,prec])
#         lstate = H.returnstate()
#         qactual = lstate[4]
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,q0-qactual,1,25,prec])
#         # while qactual<q0:
#         #     H.general_inc([3,-1,0,0,0,0,0,1,0.0,0.00001,1,5,prec])
#         #     lstate = H.returnstate()
#         #     qactual = lstate[4]
#     paccu = []
#     H.general_inc([0,0,0,0,1,0,0,1,0.0,eqc,1,25,prec])
#     paccu.append(lstate[3])
#     for i in range(n):
#         print("cycle:",i+1)
#         H.general_inc([0,0,0,0,1,0,0,1,0.0,-2*eqc,1,25,prec])
#         lstate = H.returnstate()
#         psig = lstate[3]
#         if psig < pcut:
#             break
#         H.general_inc([0,0,0,0,1,0,0,1,0.0,2*eqc,1,25,prec])
#         lstate = H.returnstate()
#         paccu.append(lstate[3])
#         psig = lstate[3]
#         if psig < pcut:
#             break
#     H.end()
#     reccu=H.returnrec()
#     w = 1
#     return(w,reccu,paccu)



# def TMCD_test(p0,e1c,n,prec,params,bug): #ini pc #cyclic e1 #cycles number #how precise #bug
#     H.const(params)
#     H.start()
            
#     H.init_stress([p0,0])    
#     H.general_inc([3,-1,0,1,0,0,0,0,0.0,0.01,1,10,4])
#     for i in range(n):
#         print('cycle:',i+1)
#         e1aim = (i+1)*e1c
#         lstate = H.returnstate()
#         e1 = lstate[2]+1/3*lstate[1]
#         H.general_inc([3,-1,0,0,0,0,1/3,1,0.0,e1aim-e1,1,50,prec])
#         lstate = H.returnstate()
#         qsig = lstate[4]
#         H.general_inc([3,-1,0,1,0,0,0,0,0.0,-qsig,1,50,prec])
#     H.end()
#     tmcd=H.returnrec()
#     w = 1
#     return(w,tmcd)

