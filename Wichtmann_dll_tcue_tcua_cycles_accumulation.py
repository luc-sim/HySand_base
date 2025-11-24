# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:12:12 2024

@author: Luc simonin

Looking at getting TCUE diminution of mean effective stress against number of cycles

USE:

[TCUEN,TCUAN,TCUANsimul]=Wicht_dl_log_cycles()
[TCUEN,TCUAN]=Wicht_dl_log_cycles()

"""
import numpy as np
import csv


def Wicht_dl_log_cycles():
    print("This program allows to put Wichtmann data into the Hyperdrive working format.")
    print("To use the data of one test, use TCUEN[i] or TCUAN[i] to access TCUE or TCUA number i.")

    #The tests imported are:
    #TCUE1 to 15 (only the small strains!)
    #TCUA1 to 19
    
    #path to the directory with all of Wichtmann data
    ######################## MODIFY TO YOUR CONVENIENCE ###########################
    path='C:\\Users\\lsimonin\\OneDrive - UCL\\Documents\\05_Constitutive modelling\\Wichtmann data\\'
    ######################## MODIFY TO YOUR CONVENIENCE ###########################

    #TCUE1 to TCUE14
    TCUEN=[[]]
    for k in range(1,15):
        input_file_name = 'TCUE'+'%d'%k+'-test-data-vs-N'+'.dat'
        input_file = open(path + input_file_name, 'r')
        Ncyc = []
        pav = []
        eqav = []
        i = 0
        for line in input_file:
            text = line.rstrip("\n\r")
            textsplit = text.split()
            i+=1
            if i>3:
                Ncyc.append(textsplit[0])
                pav.append(textsplit[2])
                E1 = - np.log( 1 - float(textsplit[4])/100)
                eqav.append(E1)
        nb = len(Ncyc)
        DATA_Ncyc=np.array([float(Ncyc[i]) for i in range(nb)])
        DATA_pav=np.array([float(pav[i]) for i in range(nb)])
        DATA_eqav=np.array([float(eqav[i]) for i in range(nb)])
        Data=[[]]
        for i in range(nb):
            Data[0].append(np.concatenate(([DATA_Ncyc[i]],[DATA_pav[i]],[DATA_eqav[i]])))
        dat=[]
        for i in range(nb):
            dat.append(Data[0][i])
        TCUEN.append([dat])
        
    #TCUA1 to TCUA19
    TCUAN=[[]]
    for k in range(1,20):
        input_file_name = 'TCUA'+'%d'%k+'-test-data-vs-N'+'.dat'
        input_file = open(path + input_file_name, 'r')
        Ncyc = []
        pav = []
        eqav = []
        i = 0
        for line in input_file:
            text = line.rstrip("\n\r\,")
            textsplit = text.split()
            i+=1
            if i>3:
                Ncyc.append(textsplit[0])
                pav.append(textsplit[2])
                E1 = - np.log( 1 - float(textsplit[4])/100)
                eqav.append(E1)
        nb = len(Ncyc)
        DATA_Ncyc=np.array([float(Ncyc[i]) for i in range(nb)])
        DATA_pav=np.array([float(pav[i]) for i in range(nb)])
        DATA_eqav=np.array([float(eqav[i]) for i in range(nb)])
        Data=[[]]
        for i in range(nb):
            Data[0].append(np.concatenate(([DATA_Ncyc[i]],[DATA_pav[i]],[DATA_eqav[i]])))
        dat=[]
        for i in range(nb):
            dat.append(Data[0][i])
        TCUAN.append([dat])
    
    # path='C:\\Users\\lsimonin\\OneDrive - UCL\\Documents\\05_Constitutive modelling\\2023_10_HySand articles\\2024_12_November revision by Guy\\code\\'
    # #TCUA1 to TCUA19
    # TCUANsimul=[[]]
    # for k in range(1,20):
    #     input_file_name = '2025_02_19_TCUANavg%i.csv'%k
    #     input_file = open(path + input_file_name, 'r')
    #     Ncyc = []
    #     pav = []
    #     eqav = []
    #     llength = [0,10,999,999,999,999,683,100,2237,999,999,687,999,1981,999,130,100,100,100,1781]
    #     i = 0
    #     reader = csv.reader(input_file)  # Uses ',' as the default delimiter
    #     for row in reader:
    #         if i>1:
    #             Ncyc.append(row[0])
    #             pav.append(row[1])
    #             eqav.append(row[2])
    #         i+=1
    #     nb = len(Ncyc)
    #     rnb = range(nb)
    #     DATA_Ncyc=np.array([float(Ncyc[i])-1 for i in rnb])
    #     DATA_pav=np.array([float(pav[i]) for i in rnb])
    #     DATA_eqav=np.array([float(eqav[i]) for i in rnb])
    #     Data=[[]]
    #     for i in range(nb):
    #         Data[0].append(np.concatenate(([DATA_Ncyc[i]],[DATA_pav[i]],[DATA_eqav[i]])))
    #     dat=[]
    #     for i in range(nb):
    #         dat.append(Data[0][i])
    #     TCUANsimul.append([dat])
    # return [TCUEN, TCUAN, TCUANsimul]
    return [TCUEN, TCUAN]
    