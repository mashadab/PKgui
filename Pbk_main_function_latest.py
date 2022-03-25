#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk main functionsolution
"""

from scipy.integrate import quad
from scipy.special import ellipk
from numpy import sin, cos, sqrt, array, pi, nan, isnan, inf, savetxt, linspace, isnan, concatenate, shape, any, argsort
from scipy.optimize import fsolve,least_squares
from matplotlib.pyplot import figure, xlabel, ylabel, plot, vlines, hlines, savefig, show, tight_layout
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from warnings import filterwarnings
from matplotlib.pyplot import rcParams
from os import makedirs
from os.path import exists
from shutil import rmtree

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
    return x - L + C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

def QbyK_func(alpha, beta,C):
 return C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0] + C*sqrt(alpha)*quad(lambda phi: (ellipk(1 - alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

def QbyK_res_func(QbyK,alpha, beta,C):
 return QbyK - C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0] - C*sqrt(alpha)*quad(lambda phi: (ellipk(1 - alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

def H1_func(H,alpha, beta,C):
 return H + H0_func(alpha, beta,C) + C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, pi/2)[0]

def H_func(alpha, beta,C):
 return C*sqrt(alpha)*quad(lambda phi: (ellipk(alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

def QH0byQH_func(alpha, beta):
 return quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0] / (sqrt(alpha)*quad(lambda phi: (ellipk(1 - alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0])



def PbK_solution_full(H0,H_full,L_full,H1,n,output_folder,Q,K,unit,Tunit):
    if ~isnan(H_full) and ~isnan(H1) and isnan(H0):
        H_scale = H1
        H = H_full/H_scale
        L = L_full/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),H1_res_func(H,alpha, beta,C))
        a = max(H,L)
        
        res = least_squares(full_equations, (0.0001, 0.1,0.1*a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        
        print("======================")
        print("Given values")
        print("======================")
        print("1.) Lake level, H: \t", H_scale*H, f"{unit}")
        print("2.) Aquifer length, L: \t", H_scale*L, f"{unit}")
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2]), f"{unit}") 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )      
        if isnan(Q) and ~isnan(K): 
            Q = K*H_scale*QbyK_func(res.x[0],res.x[1],res.x[2]) 
            QbyK = H_scale*QbyK_func(res.x[0],res.x[1],res.x[2]) 
            print("5.) Specific discharge, Q: \t \t", Q, f'{unit}^2/{Tunit}') 
        elif ~isnan(Q) and isnan(K): 
            K = Q/(H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])) 
            QbyK = H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])
            print("5.) Hydraulic conductivity, K: \t \t", K ,f'{unit}/{Tunit}') 
        elif isnan(Q) and isnan(K): 
            QbyK = H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])     
            print("5.) Neither K or Q given") 

    elif ~isnan(H_full) and isnan(H1) and isnan(H0):
        H_scale = H_full
        H = H_full/H_scale
        L = L_full/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,10)))
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Higher lake level, H1: \t", H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]),f"{unit}")   
        
    elif isnan(H_full) and ~isnan(H1) and isnan(H0):        
        H_scale = H1
        H1 = H1/H_scale
        L = L_full/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H1_res_func(H_func(alpha, beta,C),alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,10)))
        H   = H_func(res.x[0],res.x[1],res.x[2])
        H_full= H_func(res.x[0],res.x[1],res.x[2])
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Lower lake level, H: \t", H_scale*H_func(res.x[0],res.x[1],res.x[2]),f"{unit}")   
        
    print(f"6.) Specific discharge over hydraulic conductivity, q/K {unit}: \t \t", H_scale*QbyK_func(res.x[0], res.x[1],res.x[2]))          
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
        if x==0: z = H1_func(H,res.x[0],res.x[1],res.x[2])
        if x==L: z = H + H0_func(res.x[0],res.x[1],res.x[2]) 
        xz_array.append([x,z])
        x_array.append(x)
        if i%(int(n/20)) == 0 and x != inf and x!=-inf and x>=0:    
            print(Psi_array[i],'\t',H_scale*x,'\t','\t','\t','\t', H_scale*z)
    x_array  = array(x_array)
    xz_array = array(xz_array)
    xz_array[x_array<0,:] = nan
    xz_array[x_array>L,:] = nan
    xz_array = H_scale*xz_array
    xz_array = xz_array[~isnan(xz_array).any(axis=1),:] 
    
    p = argsort(xz_array[:,0])
    
    xz_array = xz_array[p,:]
    
    
    
    fig = figure(figsize=(6,6) , dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(xz_array[:,0],xz_array[:,1],'b-')
    ax.vlines(0,0,H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]),colors='blue') 
    ax.vlines(H_scale*L,0,H_scale*H+H_scale*H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
    ax.vlines(H_scale*L,0,H_scale*H,colors='blue') 
    ax.hlines(H_scale*H,H_scale*L,H_scale*1.1*L,colors='blue')   
    ax.hlines(0,0,H_scale*1.1*L,colors='blue')   
    ax.set_xlabel(f'x [{unit}]')
    ax.set_ylabel(f'z [{unit}]')
    tight_layout(pad=1, w_pad=0.8, h_pad=1)
    #fig.show()

    H0 = H_scale*H0_func(res.x[0],res.x[1],res.x[2]) 
    H1 = H_scale*H1_func(H,res.x[0],res.x[1],res.x[2])

    path = f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}"
    if not exists(path):
        makedirs(path)
    else:
        rmtree(path)           # Removes all the subdirectories!
        makedirs(path) 
    
    if not output_folder =='/tmp':
        fig.savefig(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profile.pdf")
    
    fig = figure(figsize=(10,10), dpi=50)
    ax = fig.add_subplot(111)
    ax.plot(xz_array[:,0],xz_array[:,1],'b-')
    ax.vlines(0,0,H_scale,colors='blue') 
    ax.vlines(H_scale*L,0,H_scale*H+H_scale*H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
    ax.vlines(H_scale*L,0,H_scale*H,colors='blue') 
    ax.hlines(H_scale*H,H_scale*L,H_scale*1.1*L,colors='blue')   
    ax.hlines(0,0,H_scale*1.1*L,colors='blue')   
    ax.set_xlabel(f'x [{unit}]')
    ax.set_ylabel(f'z [{unit}]')
    tight_layout(pad=1, w_pad=0.8, h_pad=1)
    #fig.show()

    fig.savefig(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profile.png")

    QbyK = QbyK_func(res.x[0],res.x[1],res.x[2]) 
    QH0byQH = QH0byQH_func(res.x[0],res.x[1])
    
    names = ['Length Unit','Time Unit','Dam length L', 'Lower lake level H', 'Upper lake level H1', 'Seepage face height H0', 'Specific discharge Q', 'Hydraulic conductivity K' ,'QbyK', 'alpha', 'beta', 'C' ,'Q_H0/Q_H']
    scores= [unit, Tunit, L_full, H_full, H1,  H0, Q, K,QbyK, res.x[0],res.x[1],res.x[2]*H_scale,QH0byQH]


    if not output_folder =='/tmp':
        savetxt(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/details.csv", [p for p in zip(names, scores)], delimiter=',', fmt='%s')
        savetxt(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profiles_XandZ.csv", xz_array, delimiter=",")

    res.x[2] = res.x[2]*H_scale

    output_folder = f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}"

    return H0, H_full, L_full, res.x, xz_array,Q,K, H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]/H_scale),QH0byQH,QbyK,output_folder

