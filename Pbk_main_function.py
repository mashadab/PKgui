#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk main functionsolution
"""

from scipy.integrate import quad
from scipy.special import ellipk
from numpy import sin, cos, sqrt, array, pi, nan, inf, savetxt, linspace, isnan, concatenate, shape, any
from scipy.optimize import fsolve,least_squares
from matplotlib.pyplot import figure, xlabel, ylabel, plot, vlines, hlines, savefig, show, tight_layout
from warnings import filterwarnings
from matplotlib.pyplot import rcParams
from os import mkdir

rcParams.update({'font.size': 22})
rcParams.update({'font.family': 'Serif'})
filterwarnings("ignore")

def L_res_func(L,alpha, beta, C):
 return L - C*quad(lambda phi: ellipk(alpha + (beta - alpha)* sin(phi)**2)/sqrt(1 - alpha - (beta - alpha)* sin(phi)**2), 0, pi/2)[0]

def H_res_func(H,alpha, beta,C):
 return H - C*sqrt(alpha)*quad(lambda phi: (ellipk(alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

def H0_func(alpha, beta,C):
 return C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0]

def H1_res_func(H,alpha, beta,C):
 return 1 - H - H0_func(alpha, beta,C) - C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, pi/2)[0]

def x_func(L,alpha, beta,C,Psi):
 return L - C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

def z_func(H,alpha, beta,C,Psi):
 return H + H0_func(alpha, beta,C) + C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

def x_res_func(L,alpha, beta,C,Psi,x):
    return x - L + C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]#- x_func(L,alpha, beta,C,Psi)


def PbK_solution(H_full,L_full,H1,n,output_folder,unit):
    
    H = H_full/H1
    L = L_full/H1
    
    def full_equations(p):
        alpha, beta, C = p
        return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),H1_res_func(H,alpha, beta,C))
    
    res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,10)))
    
    print("======================")
    print("Given values")
    print("======================")
    print("1.) Lake level, H: \t", H1*H)
    print("2.) Aquifer length, L: \t", H1*L)
    print("======================")
    print("Output values")
    print("======================")
    print("1.) Seepage face height, H0: \t", H1*H0_func(res.x[0],res.x[1],res.x[2])) 
    print("2.) alpha: \t", res.x[0])
    print("3.) beta: \t", res.x[1])
    print("4.) C: \t \t", res.x[2] )
    
    
    
    
    print("======================")
    print("Free surface profiles")
    print("======================")
    print(f"Psi \t x {unit} \t z {unit}")
    Psi_array = linspace(0,10,n+1)
    x_array = []
    z_array = []
    xz_array = []
    for i in range(0,n):
        x  = x_func(L,res.x[0], res.x[1],res.x[2],Psi_array[i])
        z  = z_func(H,res.x[0], res.x[1],res.x[2],Psi_array[i]) 
        if x==0: z = 1
        if x==L: z = H + H0_func(res.x[0],res.x[1],res.x[2]) 
        xz_array.append([x,z])
        x_array.append(x)
        if i%(int(n/20)) == 0 and x != inf and x!=-inf and x>=0:    
            print(Psi_array[i],'\t',H1*x,'\t','\t','\t','\t', H1*z)
    x_array  = array(x_array)
    xz_array = array(xz_array)
    xz_array[x_array<0,:] = nan
    xz_array[x_array>L,:] = nan
    xz_array = H1*xz_array
    fig = figure(figsize=(6,6) , dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(H1*xz_array[:,0],H1*xz_array[:,1],'b-')
    ax.vlines(0,0,H1,colors='blue') 
    ax.vlines(H1*L,0,H1*H+H1*H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
    ax.vlines(H1*L,0,H1*H,colors='blue') 
    ax.hlines(H1*H,H1*L,H1*1.1*L,colors='blue')   
    ax.hlines(0,0,H1*1.1*L,colors='blue')   
    ax.set_xlabel(f'x [{unit}]')
    ax.set_ylabel(f'z [{unit}]')
    tight_layout(pad=1, w_pad=0.8, h_pad=1)
    #fig.show()
    
    H0 = H0_func(res.x[0],res.x[1],res.x[2]) 
    
    mkdir(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}")
    fig.savefig(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profile.pdf")
    
    fig = figure(figsize=(10,10), dpi=50)
    ax = fig.add_subplot(111)
    ax.plot(xz_array[:,0],xz_array[:,1],'b-')
    ax.vlines(0,0,H1,colors='blue') 
    ax.vlines(H1*L,0,H1*H+H1*H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
    ax.vlines(H1*L,0,H1*H,colors='blue') 
    ax.hlines(H1*H,H1*L,H1*1.1*L,colors='blue')   
    ax.hlines(0,0,H1*1.1*L,colors='blue')   
    ax.set_xlabel(f'x [{unit}]')
    ax.set_ylabel(f'z [{unit}]')
    tight_layout(pad=1, w_pad=0.8, h_pad=1)
    #fig.show()
    
    H0 = H0_func(res.x[0],res.x[1],res.x[2]) 
    fig.savefig(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profile.png")

    names = ['Unit','Dam length L', 'Lower lake level H', 'Upper lake level H1', 'Seepage face height H0', 'alpha', 'beta', 'C' ]
    scores = [unit, L_full, H_full, H1,  H0*H1, res.x[0],res.x[1],res.x[2] ]
    
    xz_array = xz_array[~isnan(xz_array).any(axis=1),:]
    
    savetxt(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/details.csv", [p for p in zip(names, scores)], delimiter=',', fmt='%s')
    savetxt(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profiles_XandZ.csv", xz_array, delimiter=",")
    
    return H0, res.x, xz_array

