# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 17:24:37 2020

@author: AdminY
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def func_001():
    
    t = np.array([0.05, 0.1, 0.15, 0.2, 0.3])
    truncated = np.array([126886526.2874752, 149422953.72723293, 161568860.58943313, 187250307.79044265, 244426540.78683034])
    full = np.array([123160952.07183249, 143668023.38116586, 152430824.33375946, 168566148.0798896, 202251348.49787498])
    plt.plot(t,truncated,'r',label='truncated')
    plt.plot(t,full,'b',label='full')
    plt.legend(loc="upper left")
    plt.xlabel(r'ion concentration(M)')
    plt.ylabel(r'$D1(Å^{2}/s)$')
    plt.savefig(r'D:\temp\XPA\D1_full_truncated.png')
    plt.show()
    
def func_002():
    discard_lines = 2000
    with_DNA_path = r'D:\temp\XPA\XPA_Full_out_02'
    without_DNA_path = r'D:\temp\XPA\XPA_Full_out_03'
    w_str = 'Rg_xpa_dna'
    wo_str = 'Rg_xpa_no-dna'
    
    os.chdir(without_DNA_path)
    all_files = os.listdir('.')
    keys = [ f.split('_')[-2] for f in all_files ]
    keys = list(set(keys))
    
    i = 0
    for k in keys:
        t1 = []
        t2 = []
        i = i + 1
        for f in all_files:
            if k in f:
                os.chdir(without_DNA_path)
                fr = open(f,'r')
                print(f,' ...reading..')
                all_lines = fr.readlines()
                fr.close()
                t1.extend([ float(each) for each in all_lines[discard_lines:]] )
                
                os.chdir(with_DNA_path)
                fr = open(f.replace(wo_str,w_str),'r')
                all_lines = fr.readlines()
                fr.close()
                t2.extend( [ float(each) for each in all_lines[discard_lines:] ] )
                
        fig = plt.figure(i)
        ax = fig.add_subplot(2,1,1)
        n, bins, patches = ax.hist(np.array(t1), bins=100, normed=1,edgecolor='None',facecolor='red',label='wo DNA') 
        ax.set_title(str(k))
        ax.set_xlabel('Rg(A)')
        ax.set_ylabel('probability')
        plt.legend(loc='upper right')
        
        ax = fig.add_subplot(2,1,2)
        #ax.set_title('with DNA' + str(k))
        ax.set_xlabel('Rg(A)')
        ax.set_ylabel('probability')
        n, bins, patches = ax.hist(np.array(t2), bins=100, normed=1,edgecolor='None',facecolor='blue',label='with DNA')
        plt.legend(loc='upper right')
        
        plt.savefig(r'D:\\temp\\XPA\\' + 'Rg_' + str(k) + '.png')
        
    plt.show()

def func_003():
    t = np.array([0.1, 0.15, 0.2, 0.3])
    truncated = np.array([0.076923076923076927, 0.00041073384446878427, 0.00047326076668244201, 0.0089126559714795012])
    full = np.array([0.0010735816334957471, 0.010133578995854445, 0.0078301722637898041, 0.0011663597298956416])
    plt.plot(t,truncated,'r',label='truncated')
    plt.plot(t,full,'b',label='full')
    plt.legend(loc="upper left")
    plt.xlabel(r'ion concentration(M)')
    plt.ylabel(r'bound rate')
    plt.savefig(r'D:\temp\XPA\rate_full_truncated.png')
    plt.show()
    
#func_001()
func_003()
    
    
    