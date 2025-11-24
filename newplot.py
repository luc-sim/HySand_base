# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 15:08:13 2022

@author: Luc Simonin

Plotting routines
"""
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from datetime import datetime
#datetime.now().strftime("%Y_%m_%d")

# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 18}
font = {'family':'Tahoma','weight':'normal'}

plt.rc('font', **font)
cols=['k','grey','b','g','r','m','orange','c','lime','pink','yellow']
qmaj, qmin = 100, -100
pmaj, pmin = 220, 0
e2maj, e2min = 0.02, -0.02 #deviatoric
e1maj, e1min = 0.03, -0.08 #volumetric
q22, e22 = 120, 0.0001
limits0=[qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]
nr=2
nc=2
names = ['\u03b5'+'$_p$','\u03b5'+'$_q$',"p'",'q']
high=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
legend00 = ['','','','','','','','','','','','','','','','','','','']


def plotD(data,model,test_col=cols,limits=limits0,legend0=legend00,lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=5,xt2=5,yt2=5,ttl="",save=False,ori='v'):
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    if ori == 'v':
        fig, axes = plt.subplots(2, 1,figsize=(4,6))
    if ori == 'h':
        fig, axes = plt.subplots(1, 2,figsize=(10,3))  
    fig.tight_layout()
    
    for i in range(len(data)):
        recl = data[i]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]-recl[0][2]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        v = [item[6] for item in recl]
        axes[1].plot(e2, s2, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
        axes[0].plot(e2, e1, test_col[i],linewidth=lw,markersize=mw) 
    
    for j in range(len(model)):
        recl = model[j]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        i = j+len(data)
        axes[1].plot(e2, s2, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
        axes[0].plot(e2, e1, test_col[i],linewidth=lw,markersize=mw)  
    fonts = 12
    fig.suptitle(ttl,fontsize=15,font='Garamond',weight='semibold')
    if ori == 'v':
        fig.legend(ncol=ncl, loc='center', prop={'size': fonts-1}, bbox_to_anchor=(0.5, 0.075),facecolor='whitesmoke')
        plt.subplots_adjust(bottom=0.25,top=0.925,left=0.23,right=0.925,wspace=0.25, hspace=0.3)
    if ori == 'h':
        fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.92, 0.5),facecolor='whitesmoke')
        plt.subplots_adjust(bottom=0.2,top=0.87,left=0.12,right=0.85,wspace=0.25, hspace=0.1)
    axes[0].set_xlim([e2min,e2maj])
    axes[0].set_ylim([e1min,e1maj])
    axes[0].set_ylabel('\u03b5'+'$_p$',fontsize=fonts+2)
    axes[0].set_xlabel('\u03b5'+'$_q$',fontsize=fonts+2)
    axes[0].tick_params(axis='x',labelsize=fonts)
    axes[0].tick_params(axis='y',labelsize=fonts)
    axes[0].xaxis.set_major_locator(ticker.MaxNLocator(xt1))
    axes[0].yaxis.set_major_locator(ticker.MaxNLocator(yt1))
    axes[0].grid(linestyle=':',linewidth=1)
    axes[0].text(0.01*(e2maj-e2min)+e2min,0.89*(e1maj-e1min)+e1min,lab[0],horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[1].set_xlim([e2min,e2maj])
    axes[1].set_ylim([qmin,qmaj])
    axes[1].set_ylabel('$\it{q}$',fontsize=fonts+2)
    axes[1].set_xlabel('\u03b5'+'$_q$',fontsize=fonts+2)
    axes[1].tick_params(axis='x',labelsize=fonts)
    axes[1].tick_params(axis='y',labelsize=fonts)
    axes[1].xaxis.set_major_locator(ticker.MaxNLocator(xt2))
    axes[1].yaxis.set_major_locator(ticker.MaxNLocator(yt2))
    axes[1].grid(linestyle=':',linewidth=1)
    axes[1].text(0.01*(e2maj-e2min)+e2min,0.89*(qmaj-qmin)+qmin,lab[1],horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'TD.svg',format='svg',dpi=600) 
        

def plotU(data,model,test_col=cols,limits=limits0,legend0=legend00,lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=5,xt2=5,yt2=5,ttl="",save=False,ori='v'):
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    if ori == 'v':
        fig, axes = plt.subplots(2, 1,figsize=(4,6))
    if ori == 'h':
        fig, axes = plt.subplots(1, 2,figsize=(10,3))
    fig.tight_layout()
    for i in range(len(data)):
        recl = data[i]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]-recl[0][2]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        v = [item[6] for item in recl]
        axes[0].plot(s1, s2, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
        axes[1].plot(e2, s2, test_col[i],linewidth=lw,markersize=mw)
    for j in range(len(model)):
        recl = model[j]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]-recl[0][2]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        i = j+len(data)
        axes[0].plot(s1, s2, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
        axes[1].plot(e2, s2, test_col[i],linewidth=lw,markersize=mw)
    fonts = 12
    fig.suptitle(ttl,fontsize=15,font='Garamond',weight='semibold')
    if ori == 'v':
        fig.legend(ncol=ncl, loc='center', prop={'size': fonts-1}, bbox_to_anchor=(0.5, 0.075),facecolor='whitesmoke')
        plt.subplots_adjust(bottom=0.23,top=0.925,left=0.23,right=0.925,wspace=0.25, hspace=0.3)
    if ori == 'h':
        fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.92, 0.5),facecolor='whitesmoke')
        plt.subplots_adjust(bottom=0.2,top=0.87,left=0.12,right=0.85,wspace=0.25, hspace=0.1)
    axes[0].set_xlim([pmin,pmaj])
    axes[0].set_ylim([qmin,qmaj])
    axes[0].set_xlabel('$\it{p}$',fontsize=fonts+2)
    axes[0].set_ylabel('$\it{q}$',fontsize=fonts+2)
    axes[0].tick_params(axis='x',labelsize=fonts)
    axes[0].tick_params(axis='y',labelsize=fonts)
    axes[0].xaxis.set_major_locator(ticker.MaxNLocator(xt1))
    axes[0].yaxis.set_major_locator(ticker.MaxNLocator(yt1))
    axes[0].grid(linestyle=':',linewidth=1)
    axes[0].text(0.01*(pmaj-pmin)+pmin,0.89*(qmaj-qmin)+qmin,lab[0],horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[1].set_xlim([e2min,e2maj])
    axes[1].set_ylim([qmin,qmaj])
    axes[1].set_xlabel('\u03b5'+'$_q$',fontsize=fonts+2)
    axes[1].set_ylabel('$\it{q}$',fontsize=fonts+2)
    axes[1].tick_params(axis='x',labelsize=fonts)
    axes[1].tick_params(axis='y',labelsize=fonts)
    axes[1].xaxis.set_major_locator(ticker.MaxNLocator(xt2))
    axes[1].yaxis.set_major_locator(ticker.MaxNLocator(yt2))
    axes[1].grid(linestyle=':',linewidth=1)
    axes[1].text(0.01*(e2maj-e2min)+e2min,0.89*(qmaj-qmin)+qmin,lab[1],horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'TU.svg',format='svg',dpi=600) 

def plotISO(data,model,test_col=cols,limits=limits0,legend0=legend00,lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=5,xt2=5,yt2=5,ttl="",save=False):
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    fig, axes = plt.subplots(1, 1,figsize=(4,3))
    fig.tight_layout()
    for i in range(len(data)):
        recl = data[i]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]-recl[0][2]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        v = [item[6] for item in recl]
        axes.plot(s1, e1, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    for j in range(len(model)):
        recl = model[j]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        i = j+len(data)
        axes.plot(s1, e1, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    fonts = 12
    fig.suptitle(ttl,fontsize=15,font='Garamond',weight='semibold')
    fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.8, 0.4),facecolor='whitesmoke')
    plt.subplots_adjust(bottom=0.2,top=0.87,left=0.21,right=0.95,wspace=0.25, hspace=0.1)
    axes.set_xlim([pmin,pmaj])
    axes.set_ylim([e1min,e1maj])
    axes.set_ylabel('\u03b5'+'$_p$',fontsize=fonts+2)
    axes.set_xlabel('p',fontsize=fonts+2)
    axes.tick_params(axis='x',labelsize=fonts)
    axes.tick_params(axis='y',labelsize=fonts)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(xt1))
    axes.yaxis.set_major_locator(ticker.MaxNLocator(yt1))
    axes.grid(linestyle=':',linewidth=1)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'TI.svg',format='svg',dpi=600) 

def plotOE(data,model,v0s,test_col=cols,limits=[1000,0.1,1.05,0.68],legend0=legend00,lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=5,xt2=5,yt2=5,ttl="",save=False):
    [pmax,pmin,emax,emin]=limits
    fig, axes = plt.subplots(1, 1,figsize=(4,5))
    fig.tight_layout()
    for i in range(len(data)):
        recl = data[i]
        e1 = [item[1] for item in recl]
        s1 = [item[3] for item in recl]
        e = [item[-1]-1 for item in recl]
        v = [item[-1] for item in recl]
        axes.plot(s1, v, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    for j in range(len(model)):
        v0 = v0s[j]
        recl = model[j]
        ep0 = recl[0][1]
        cyc = [item[0]/200 for item in recl]
        e1 = [item[1]/3+item[2] for item in recl]-(recl[0][1]/3+recl[0][2])
        e3 = [item[1]/3-item[2]/2 for item in recl]-(recl[0][1]/3-recl[0][2]/2)
        s1 = [item[3]+2/3*item[4] for item in recl]
        s3 = [item[3]-1/3*item[4] for item in recl]
        p = [item[3] for item in recl]
        q = [item[4] for item in recl]
        e = [v0*np.exp(-(item[1]-ep0))-1 for item in recl]
        v = [v0*np.exp(-(item[1]-ep0)) for item in recl]
        i = j+len(data)
        axes.plot(s1, v, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    fonts = 12
    fig.suptitle(ttl,fontsize=15,font='Garamond',weight='semibold')
    # fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.875, 0.5),facecolor='whitesmoke')
    plt.subplots_adjust(bottom=0.1,top=0.92,left=0.18,right=0.93,wspace=0.25, hspace=0.1)
    axes.set_xscale('log')
    axes.set_xlim([pmin,pmax])
    # axes.set_ylim([v0-1-0.1,v0-1+0.02])
    axes.set_ylim([emin,emax])
    axes.set_ylabel('Specific volume '+'$\it{v}$',fontsize=fonts+2)
    axes.set_xlabel('$\sigma$'+'$_1$',fontsize=fonts+2)
    axes.tick_params(axis='x',labelsize=fonts)
    axes.tick_params(axis='y',labelsize=fonts)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(xt1))
    axes.xaxis.set_major_formatter(ticker.ScalarFormatter())
    axes.yaxis.set_major_locator(ticker.MaxNLocator(yt1))
    axes.grid(linestyle=':',linewidth=1)
    axes.xaxis.set_ticks([int(10*0.1*10**i)/10 for i in range(5)])
    # axes[0].text(0.01*(e2maj-e2min)+e2min,0.89*(e1maj-e1min)+e1min,'a)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'TO.svg',format='svg',dpi=600) 

def plotTCUAN(data,model,test_col=cols,limits=[500,0,400,0,0],legend0=legend00,lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=5,xt2=5,yt2=5,ttl="",save=False):
    [Ntot,pmin,pmaj,option,zeroset]=limits #option=0 to plot pav, option=1 to plot eqav; zeroset to initial value of eqav from post-drained start
    fig, axes = plt.subplots(1, 1,figsize=(4,3))
    fig.tight_layout()
    for i in range(len(data)):
        recl = data[i]
        N = [item[0] for item in recl]
        pav = [item[1] for item in recl]
        eqav = [item[2] for item in recl]
        if option == 0:
            axes.plot(N, pav, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
        if option == 1:
            axes.plot(N, eqav, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    for j in range(len(model)):
        recl = model[j]
        N = [item[0] for item in recl]
        pav = [item[1] for item in recl]
        eqav = [item[2] - zeroset for item in recl]
        i = j+len(data)
        if option == 0:
            axes.plot(N, pav, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
        if option == 1:
            axes.plot(N, eqav, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    fonts = 12
    fig.suptitle(ttl,fontsize=15,font='Garamond',weight='semibold')
    if option == 0:
        fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.8, 0.7),facecolor='whitesmoke')
    if option == 1:
        fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.37, 0.78),facecolor='whitesmoke')  
        # fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.37, 0.3),facecolor='whitesmoke')  
    plt.subplots_adjust(bottom=0.2,top=0.87,left=0.21,right=0.95,wspace=0.25, hspace=0.1)
    axes.set_xlim([0,Ntot])
    axes.set_ylim([pmin,pmaj])
    if option == 0:
        axes.set_ylabel('p'+'$_{av}$',fontsize=fonts+2)
    if option == 1:
        axes.set_ylabel('\u03b5'+'$_q$'+'$_{av}$',fontsize=fonts+2)
    axes.set_xlabel('Cycle number',fontsize=fonts+2)
    axes.tick_params(axis='x',labelsize=fonts)
    axes.tick_params(axis='y',labelsize=fonts)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(xt1))
    axes.yaxis.set_major_locator(ticker.MaxNLocator(yt1))
    axes.grid(linestyle=':',linewidth=1)
    if save==1:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'TCUAN.svg',format='svg',dpi=600) 
    
def plotTCUEN(data,model,test_col=cols,limits=[500,0,400],legend0=legend00,lab=['a)','b)'],lw=1.0,mw=1.0,ncl=3,xt1=5,yt1=5,xt2=5,yt2=5,ttl="",save=False):
    [Ntot,pmin,pmaj]=limits
    fig, axes = plt.subplots(1, 1,figsize=(4,3))
    fig.tight_layout()
    for i in range(len(data)):
        recl = data[i]
        N = [item[0] for item in recl]
        pav = [item[1] for item in recl]
        axes.plot(N, pav, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    for j in range(len(model)):
        recl = model[j]
        N = [item[0] for item in recl]
        pav = [item[1] for item in recl]
        i = j+len(data)
        axes.plot(N, pav, test_col[i],label=legend0[i],linewidth=lw,markersize=mw)
    fonts = 12
    fig.suptitle(ttl,fontsize=15,font='Garamond',weight='semibold')
    fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.8, 0.7),facecolor='whitesmoke')
    plt.subplots_adjust(bottom=0.2,top=0.87,left=0.21,right=0.95,wspace=0.25, hspace=0.1)
    axes.set_xlim([0,Ntot])
    axes.set_ylim([pmin,pmaj])
    axes.set_ylabel('p'+'$_{av}$',fontsize=fonts+2)
    axes.set_xlabel('Cycle number',fontsize=fonts+2)
    axes.tick_params(axis='x',labelsize=fonts)
    axes.tick_params(axis='y',labelsize=fonts)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(xt1))
    axes.yaxis.set_major_locator(ticker.MaxNLocator(yt1))
    axes.grid(linestyle=':',linewidth=1)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'TCUEN.svg',format='svg',dpi=600) 


# Two following functions used for calibration of m, kr and lambdas
def plotISOloglog(data,model,test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncl=5,save=False): #identify m
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    fig, axes = plt.subplots(1, 1,figsize=(6,2.5))
    fig.tight_layout()
    pr = 100 # Reference pressure
    for i in range(len(data)):
        recl = data[i]
        e1 = [item[1] for item in recl]-recl[0][1]
        lne1 = [np.log(e1[i]) for i in range(len(e1))]
        e2 = [item[2] for item in recl]-recl[0][2]
        s1 = [item[3] for item in recl]
        lns1 = [np.log(s1[i]/pr) for i in range(len(s1))]
        s2 = [item[4] for item in recl]
        v = [item[6] for item in recl]
        axes.plot(lns1, lne1, test_col[i],label=legend0[i],linewidth=lw)
    fonts = 12
    fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.9, 0.5),facecolor='whitesmoke')
    M = 5
    plt.subplots_adjust(bottom=0.22,top=0.95,left=0.15,right=0.8,wspace=0.25, hspace=0.1)
    # axes.set_xlim([pmin,pmaj])
    axes.set_ylim([-9,-3])
    axes.set_ylabel('ln('+'\u03b5'+'$_p$'+')',fontsize=fonts+2)
    axes.set_xlabel('ln(p/p'+'$_r$)',fontsize=fonts+2)
    axes.tick_params(axis='x',labelsize=fonts)
    axes.tick_params(axis='y',labelsize=fonts)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(M))
    axes.yaxis.set_major_locator(ticker.MaxNLocator(M))
    axes.grid(linestyle=':',linewidth=1)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'Isologlog.svg',format='svg',dpi=600) 

def plotISOpow(data,model,test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncl=5,save=False):
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    fig, axes = plt.subplots(1, 1,figsize=(6,2.5))
    fig.tight_layout()
    mm = float(input('Input the non-linear power m:'))
    for i in range(len(data)):
        recl = data[i]
        ep = [item[1] for item in recl]-recl[0][1]
        p = [item[3] for item in recl]
        pows1 = [1/(1-mm)*(p[i]/100)**(1-mm) for i in range(len(p))]
        axes.plot(pows1, ep, test_col[i],label=legend0[i],linewidth=lw)
    fonts = 12
    fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.9, 0.5),facecolor='whitesmoke')
    M = 5
    plt.subplots_adjust(bottom=0.25,top=0.95,left=0.18,right=0.8,wspace=0.25, hspace=0.1)    
    # axes.set_ylim([e1min,e1maj])
    axes.set_ylabel('-'+'\u03b5'+'$_p$',fontsize=fonts+2)
    axes.set_xlabel(r'$\frac{1}{1-m}$'+'('+r'$\frac{p}{pr}$'+')'+'$^{(1-m)}$',fontsize=fonts+2)
    axes.tick_params(axis='x',labelsize=fonts)
    axes.tick_params(axis='y',labelsize=fonts)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(M))
    axes.yaxis.set_major_locator(ticker.MaxNLocator(M))
    axes.grid(linestyle=':',linewidth=1)
    axes.set_ylim([0,0.025])
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'ISOpow.svg',format='svg',dpi=600) 

def plotall(data,model,test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncl=5):
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    fig, axes = plt.subplots(2, 2,figsize=(10,7),sharex='col',sharey='row')
    fig.tight_layout()
    for i in range(len(data)):
        recl = data[i]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]-recl[0][2]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        v = [item[6] for item in recl]
        axes[0,0].plot(s1, s2, test_col[i],label=legend0[i],linewidth=lw)
        axes[0,1].plot(e2, s2, test_col[i],linewidth=lw)
        axes[1,0].plot(s1, e1, test_col[i],linewidth=lw)
        axes[1,1].plot(e2, e1, test_col[i],linewidth=lw) 
    for j in range(len(model)):
        recl = model[j]
        e1 = [item[1] for item in recl]-recl[0][1]
        e2 = [item[2] for item in recl]
        s1 = [item[3] for item in recl]
        s2 = [item[4] for item in recl]
        i = j+len(data)
        axes[0,0].plot(s1, s2, test_col[i],label=legend0[i],linewidth=lw)
        axes[0,1].plot(e2, s2, test_col[i],linewidth=lw)
        axes[1,0].plot(s1, e1, test_col[i],linewidth=lw)
        axes[1,1].plot(e2, e1, test_col[i],linewidth=lw)   
    fonts = 12
    ncl0 = ncl
    if (len(data)+len(model))%ncl==0:
        nrow = (len(data)+len(model))//ncl
    else:
        nrow = (len(data)+len(model))//ncl+1
    fig.legend(ncl=ncl0, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.5, 0.975-(nrow-1)*0.025),facecolor='whitesmoke')
    M = 5
    ntickse1 = ticker.MaxNLocator(M)
    ntickse2 = ticker.MaxNLocator(M)
    nticksp = ticker.MaxNLocator(M)
    nticksq = ticker.MaxNLocator(M)
    plt.subplots_adjust(bottom=0.08,top=1.0-nrow*0.05,left=0.1,right=0.95,wspace=0.15, hspace=0.1)
    axes[0,0].set_xlim([pmin,pmaj])
    axes[0,0].set_ylim([qmin,qmaj])
    axes[0,0].set_ylabel('$\it{q}$',fontsize=fonts+2)
    axes[0,0].tick_params(axis='y',labelsize=fonts)
    axes[0,0].yaxis.set_major_locator(nticksq)
    axes[0,0].grid(linestyle=':',linewidth=1)
    axes[0,0].text(0.02*(pmaj-pmin)+pmin,0.92*(qmaj-qmin)+qmin,'a)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[0,1].set_xlim([e2min,e2maj])
    axes[0,1].set_ylim([qmin,qmaj])
    axes[0,1].grid(linestyle=':',linewidth=1)
    axes[0,1].text(0.02*(e2maj-e2min)+e2min,0.92*(qmaj-qmin)+qmin,'b)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[1,0].set_xlim([pmin,pmaj])
    axes[1,0].set_ylim([e1min,e1maj])
    axes[1,0].set_xlabel('$\it{p}$',fontsize=fonts+2)
    axes[1,0].set_ylabel('\u03b5'+'$_p$',fontsize=fonts+2)
    axes[1,0].tick_params(axis='x',labelsize=fonts)
    axes[1,0].tick_params(axis='y',labelsize=fonts)
    axes[1,0].xaxis.set_major_locator(nticksp)
    axes[1,0].yaxis.set_major_locator(ntickse1)
    axes[1,0].grid(linestyle=':',linewidth=1)
    axes[1,0].text(0.01*(pmaj-pmin)+pmin,0.92*(e1maj-e1min)+e1min,'c)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[1,1].set_xlim([e2min,e2maj])
    axes[1,1].set_ylim([e1min,e1maj])
    axes[1,1].set_xlabel('\u03b5'+'$_q$',fontsize=fonts+2)
    axes[1,1].tick_params(axis='x',labelsize=fonts)
    axes[1,1].xaxis.set_major_locator(ntickse2)
    axes[1,1].grid(linestyle=':',linewidth=1)
    axes[1,1].text(0.01*(e2maj-e2min)+e2min,0.92*(e1maj-e1min)+e1min,'d)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)

def plotpob(data,model,test_col=cols,limits=limits0,legend0=legend00,lw=1.0,ncol=5,save=False): # Plotting for small looping cycles
    [qmaj,qmin,pmaj,pmin,e2maj,e2min,e1maj,e1min]=limits
    fig, axes = plt.subplots(1, 3,figsize=(10,3))
    fig.tight_layout()
    for j in range(len(model)):
        recl = model[j]
        cyc = [item[0]/200 for item in recl]
        ep = [item[1] for item in recl]
        e1 = [item[1]/3+item[2] for item in recl]-(recl[0][1]/3+recl[0][2])
        e3 = [item[1]/3-item[2]/2 for item in recl]-(recl[0][1]/3-recl[0][2]/2)
        p = [item[3] for item in recl]
        q = [item[4] for item in recl]
        e = [1.857*np.exp(-(item[1]-ep[100]))-1 for item in recl]
        i = j+len(data)
        axes[0].plot(e1*100, q, 'b',label=legend0[i],linewidth=lw)
        axes[1].plot(cyc, e, 'b',linewidth=lw)
        axes[2].plot(cyc, e1*100, 'c',linewidth=lw)
        axes[2].plot(cyc, e3*100, 'b',linewidth=lw)
    fonts = 12
    M = 5
    # fig.legend(ncol=1, loc='center', prop={'size': fonts-2}, bbox_to_anchor=(0.92, 0.5),facecolor='whitesmoke')
    plt.subplots_adjust(bottom=0.2,top=0.87,left=0.065,right=0.99,wspace=0.32, hspace=0.1)
    axes[0].set_xlim([-0.05,0.4])
    axes[0].set_ylim([70,130])
    axes[0].set_ylabel('$\it{q}$',fontsize=fonts+2)
    axes[0].set_xlabel('\u03b5'+'$_1$'+'[%]',fontsize=fonts+2)
    axes[0].tick_params(axis='x',labelsize=fonts)
    axes[0].tick_params(axis='y',labelsize=fonts)
    axes[0].xaxis.set_major_locator(ticker.MaxNLocator(6))
    axes[0].yaxis.set_major_locator(ticker.MaxNLocator(7))
    # axes[0].grid(linestyle=':',linewidth=1)
    # axes[0].text(0.01*(e2maj-e2min)+e2min,0.89*(e1maj-e1min)+e1min,'a)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[1].set_xlim([0,30])
    axes[1].set_ylim([0.846,0.858])
    axes[1].set_ylabel('void ratio '+'$\it{e}$',fontsize=fonts+2)
    axes[1].set_xlabel('Number of cycles '+'$\it{N}$',fontsize=fonts+2)
    axes[1].tick_params(axis='x',labelsize=fonts)
    axes[1].tick_params(axis='y',labelsize=fonts)
    axes[1].xaxis.set_major_locator(ticker.MaxNLocator(3))
    axes[1].yaxis.set_major_locator(ticker.MaxNLocator(7))
    # axes[1].grid(linestyle=':',linewidth=1)
    # axes[1].text(0.01*(e2maj-e2min)+e2min,0.89*(qmaj-qmin)+qmin,'b)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    axes[2].set_xlim([0,30])
    axes[2].set_ylim([-0.3,1])
    axes[2].set_ylabel('\u03b5'+'$_1$'+', '+'\u03b5'+'$_3$'+'[%]',fontsize=fonts+2)
    axes[2].set_xlabel('Number of cycles '+'$\it{N}$',fontsize=fonts+2)
    axes[2].tick_params(axis='x',labelsize=fonts)
    axes[2].tick_params(axis='y',labelsize=fonts)
    axes[2].xaxis.set_major_locator(ticker.MaxNLocator(3))
    axes[2].yaxis.set_major_locator(ticker.MaxNLocator(7))
    # axes[2].grid(linestyle=':',linewidth=1)
    # axes[1].text(0.01*(e2maj-e2min)+e2min,0.89*(qmaj-qmin)+qmin,'b)',horizontalalignment='left',verticalalignment='bottom',fontsize=fonts+1)
    fig.text(0.9, 0.26,'\u03b5'+'$_3$'+'\n\n'+'\u03b5'+'$_1$',fontsize=14)
    if save:
        plt.savefig(datetime.now().strftime("%Y_%m_%d")+'pob.svg',format='svg',dpi=600)
