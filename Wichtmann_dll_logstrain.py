# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 16:11:33 2022

@author: scat7969

Here, conversion to log scale!

Line to run once charged:

[OElog,TMDlog,TMUlog,TMUMTlog,TMUAPlog,TCUIlog,TCUAlog,TCUElog,OEClog,ISOlog,TMCDlog,TMCUlog,TPlog]=Wicht_dl_log()
[OE,TMD,TMU,TMUMT,TMUAP,TCUI,TCUA,TCUE,OEC,ISO,TMCD,TMCU,TP]=Wicht_dl_log()

plotD(TMD[3]+TMD[13]+TMD[23]+TMDlog[3]+TMDlog[13]+TMDlog[23],[],['b:','g:','r:','b','g','r'],[900,-00,250,0,0.3,-0.0,0.03,-0.12],['linear loose','linear medium','linear dense','log loose','log medium','log dense'])

plotISO(ISO[4]+ISO[5]+ISO[6]+ISOlog[4]+ISOlog[5]+ISOlog[6],[],['b:','g:','r:','b','g','r'],[900,-00,850,0,0.3,-0.0,0.025,-0],['linear loose','linear medium','linear dense','log loose','log medium','log dense'])

"""
import numpy as np



def Wicht_dl_log():
    print("This program allows to put Wichtmann data into the Hyperdrive working format.")
    print("To use the data of one test, TMU[5] accesses the test data of TMU5.")
    print("Only exceptions are for TMU-MT and TMU-AP where TMUMT[] and TMUAP[] should be called.")
    print("This code takes under a minute to run.")
    #The tests imported are:
    #OE1 to OE12
    #TMD1 to TMD25
    #TMU1 to TMU12
    #TMU-MT1 to TMU-MT9
    #TMU-AP1 to TMU-AP3
    #TCUI1 to TCUI22
    #TCUA1 to TCUA19
    #TCUE1 to TCUE23
    #OEC1 to OEC4
    #ISO1 to ISO6
    #TMCD1 to TMCD7
    #TMCU1 to TMCU7
    #TP4,6,8,9,10,11
    
    #path to the directory with all of Wichtmann data
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    # path='C:\\Users\\scat7969\\OneDrive - Nexus365\\Research\\Wichtmann data\\'
    # C:\Users\lsimonin\OneDrive - UCL\Documents\05_Constitutive modelling
    path='C:\\Users\\lsimonin\\OneDrive - UCL\\Documents\\05_Constitutive modelling\\Wichtmann data\\'
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    
    #Define maximum number of points taken per test
    nb=30000
    #If number of points is inferior, it will take all points
    #This makes the data noisy for some tests, but the code can be modified to take
    
    
    #OE1 to OE12
    # WARNING # in this data I placed sigma1 in the p' column of rec
    # and e1 instead of eq
    OElog=[[]]
    for k in range(1,13):
        input_file_name = 'OE'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1.append(-np.log(1-float(textsplit[1])/100))
                S1.append(textsplit[0])
                e.append(textsplit[2])
                t.append(i-3)
    
        length=len(E1)
        interv=length//nb
        if interv>0:
            DATA_E1=np.array([float(E1[interv*i]) for i in range(nb)])
            DATA_E3=np.array([0 for i in range(nb)])
            DATA_P=np.array([float(S1[interv*i]) for i in range(nb)])
            DATA_Q=np.array([0 for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
            Data=[[]]
            for i in range(nb):
                Data[0].append(np.concatenate(([DATA_t[i]],[DATA_E1[i]],[DATA_E3[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
    
        else:
            DATA_E1=np.array([float(E1[i]) for i in range(length)])
            DATA_E3=np.array([0 for i in range(length)])
            DATA_P=np.array([float(S1[i]) for i in range(length)])
            DATA_Q=np.array([0 for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
            Data=[[]]
            for i in range(length):
                Data[0].append(np.concatenate(([DATA_t[i]],[DATA_E1[i]],[DATA_E3[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        OElog.append([dat])
    
    
    #TMD1 to TMD25
    TMDlog=[[]]
    for k in range(1,26):
        input_file_name = 'TMD'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                Ep.append( -np.log( 1 - float(textsplit[1])/100))
                Es.append(  np.log( ( 1 - float(textsplit[1])/100)**(1/3) / (1-float(textsplit[0])/100) ) )
                P.append(textsplit[6])
                Q.append(textsplit[5])
                e.append(textsplit[4]) 
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([i/nb],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TMDlog.append([dat])
    
    
    #TMU1 to TMU12
    TMUlog=[[]]
    for k in range(1,13):
        input_file_name = 'TMU'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[0])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([i/nb],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[0])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TMUlog.append([dat])
    
    
    #TMU-MT1 to TMU-MT9
    TMUMTlog=[[]]
    for k in range(1,10):
        input_file_name = 'TMU-MT'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[0])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([i/nb],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[0])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TMUMTlog.append([dat])
    
    #TMU-AP1 to TMU-AP3
    TMUAPlog=[[]]
    for k in range(1,4):
        input_file_name = 'TMU-AP'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[0])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([i/nb],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[0])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TMUAPlog.append([dat])
    
    
    #TCUI1 to TCUI22
    TCUIlog=[[]]
    for k in range(1,23):
        input_file_name = 'TCUI'+'%d'%k+'-test-data-vs-time'+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[8])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
                t.append(textsplit[0])
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([DATA_t[i]],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[0])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TCUIlog.append([dat])
    
    
    #TCUA1 to TCUA19
    TCUAlog=[[]]
    for k in range(1,20):
        input_file_name = 'TCUA'+'%d'%k+'-test-data-vs-time'+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[8])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
                e.append(0)
                t.append(textsplit[0])            
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([DATA_t[i]],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TCUAlog.append([dat])
    
    
    #TCUE1 to TCUE23
    TCUElog=[[]]
    for k in range(1,24):
        input_file_name = 'TCUE'+'%d'%k+'-test-data-vs-time'+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[8])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
                Ep.append(0)
                e.append(0)
                t.append(textsplit[0])            
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([DATA_t[i]],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TCUElog.append([dat])
    
    
    #OEC1 to OEC4
    # WARNING # in this data I placed sigma1 in the p' column of rec
    # Also it is e1 and not eq
    OEClog=[[]]
    for k in range(1,5):
        input_file_name = 'OEC'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1.append(-np.log(1-float(textsplit[1])/100))
                S1.append(textsplit[0])
                e.append(textsplit[2])
                t.append(i-3)
        length=len(E1)
        interv=length//nb
        if interv>0:
            DATA_E1=np.array([float(E1[interv*i]) for i in range(nb)])
            DATA_E3=np.array([0 for i in range(nb)])
            DATA_P=np.array([float(S1[interv*i]) for i in range(nb)])
            DATA_Q=np.array([0 for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
            Data=[[]]
            for i in range(nb):
                Data[0].append(np.concatenate(([DATA_t[i]],[DATA_E1[i]],[DATA_E3[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        else:
            DATA_E1=np.array([float(E1[i]) for i in range(length)])
            DATA_E3=np.array([0 for i in range(length)])
            DATA_P=np.array([float(S1[i]) for i in range(length)])
            DATA_Q=np.array([0 for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
            Data=[[]]
            for i in range(length):
                Data[0].append(np.concatenate(([DATA_t[i]],[DATA_E1[i]],[DATA_E3[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        OEClog.append([dat])
    
    
    #ISO1 to ISO6
    ISOlog=[[]]
    for k in range(1,7):
        input_file_name = 'ISO'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                Ep.append(-np.log( 1 - float(textsplit[1])/100))
                Es.append(0)
                P.append(textsplit[0])
                Q.append(0)
                e.append(textsplit[2])
                t.append(textsplit[0])            
        length=len(Ep)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([DATA_t[i]],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        ISOlog.append([dat])
    
    
    #TMCD1 to TMCD7
    TMCDlog=[[]]
    for k in range(1,8):
        input_file_name = 'TMCD'+'%d'%k+"-test-data-vs-time"+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                Ep.append(-np.log( 1 - float(textsplit[2])/100 ))
                Es.append(  np.log( ( 1 - float(textsplit[2])/100)**(1/3) / (1-float(textsplit[1])/100) ) )
                P.append(textsplit[4])
                Q.append(textsplit[5])
                e.append(textsplit[3])
                t.append(textsplit[0])            
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
            DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
            DATA_t=np.array([float(t[i]) for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([DATA_t[i]],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TMCDlog.append([dat])
    
    
    #TMCU1 to TMCU6
    TMCUlog=[[]]
    for k in range(1,7):
        input_file_name = 'TMCU'+'%d'%k+"-test-data-vs-time"+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                E1 = - np.log( 1 - float(textsplit[8])/100)
                Es.append(E1)
                Ep.append(0)
                P.append(textsplit[6])
                Q.append(textsplit[7])
                t.append(textsplit[0])
    
        length=len(Es)
        interv=length//nb
        # if interv>0:
        #     DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
        #     DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
        #     DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
        #     DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
        #     DATA_t=np.array([float(t[interv*i]) for i in range(nb)])
        # else:
        DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
        DATA_Es=np.array([float(Es[i]) for i in range(length)])
        DATA_P=np.array([float(P[i]) for i in range(length)])
        DATA_Q=np.array([float(Q[i]) for i in range(length)])
        DATA_t=np.array([float(t[i]) for i in range(length)])
        Data=[[]]
        # for i in range(min(nb,length)):
        for i in range(length):
            Data[0].append(np.concatenate(([DATA_t[i]],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[0])))
        dat=[]
        # for i in range(min(nb,length)):
        for i in range(length):
            dat.append(Data[0][i])
        TMCUlog.append([dat])
        
        
        
    #TP4,6,8,9,10,11
    TPlog=[[0] for i in range(12)]
    for k in [4,6,8,9,10,11]:
        input_file_name = 'TP'+'%d'%k+'.dat'
        input_file = open(path + input_file_name, 'r')
        t=[]
        E1=[]
        Ep=[]
        Es=[]
        S1=[]
        P=[]
        Q=[]
        e=[]
        L=[Ep,Es,P,Q]
        i=0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>10:
                Ep.append(-np.log( 1 - float(textsplit[1])/100 ))
                Es.append(  np.log( ( 1 - float(textsplit[1])/100)**(1/3) / (1-float(textsplit[0])/100) ) )
                P.append(textsplit[5])
                Q.append(textsplit[6])
                e.append(0)
    
        length=len(Es)
        interv=length//nb
        if interv>0:
            DATA_Ep=np.array([float(Ep[interv*i]) for i in range(nb)])
            DATA_Es=np.array([float(Es[interv*i]) for i in range(nb)])
            DATA_P=np.array([float(P[interv*i]) for i in range(nb)])
            DATA_Q=np.array([float(Q[interv*i]) for i in range(nb)])
            DATA_v=np.array([float(e[interv*i])+1 for i in range(nb)])
        else:
            DATA_Ep=np.array([float(Ep[i]) for i in range(length)])
            DATA_Es=np.array([float(Es[i]) for i in range(length)])
            DATA_P=np.array([float(P[i]) for i in range(length)])
            DATA_Q=np.array([float(Q[i]) for i in range(length)])
            DATA_v=np.array([float(e[i])+1 for i in range(length)])
        Data=[[]]
        for i in range(min(nb,length)):
            Data[0].append(np.concatenate(([i/nb],[DATA_Ep[i]],[DATA_Es[i]],[DATA_P[i]],[DATA_Q[i]],[0],[DATA_v[i]])))
        dat=[]
        for i in range(min(nb,length)):
            dat.append(Data[0][i])
        TPlog[k]=[dat]
        
        
    return [OElog,TMDlog,TMUlog,TMUMTlog,TMUAPlog,TCUIlog,TCUAlog,TCUElog,OEClog,ISOlog,TMCDlog,TMCUlog,TPlog]
    
    
    
    
    