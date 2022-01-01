#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk solution
"""

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
from scipy.optimize import fsolve,least_squares
import math
import matplotlib.pyplot as plt

def L_res_func(L,alpha, beta, C):
 return L - C*integrate.quad(lambda phi: special.ellipk(alpha + (beta - alpha)* np.sin(phi)**2)/np.sqrt(1 - alpha - (beta - alpha)* np.sin(phi)**2), 0, np.pi/2)[0]

def H_res_func(H,alpha, beta,C):
 return H - C*np.sqrt(alpha)*integrate.quad(lambda phi: (special.ellipk(alpha * np.sin(phi)**2)*np.sin(phi))/np.sqrt((1 - alpha* np.sin(phi)**2)*(beta - alpha* np.sin(phi)**2)), 0, np.pi/2)[0]

def H0_func(alpha, beta,C):
 return C*integrate.quad(lambda phi: (special.ellipk(np.cos(phi)**2)*np.sin(phi)*np.cos(phi))/np.sqrt((1 - (1-alpha)* np.sin(phi)**2)*(1 - (1-beta)* np.sin(phi)**2)), 0, np.pi/2)[0]

def H1_res_func(H,alpha, beta,C):
 return 1 - H - H0_func(alpha, beta,C) - C*integrate.quad(lambda phi: (special.ellipk(np.cos(phi)**2)*np.sin(phi))/np.sqrt((1 - alpha* np.sin(phi)**2)*(1 - beta* np.sin(phi)**2)), 0, np.pi/2)[0]

def full_equations(p):
    alpha, beta, C = p
    return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),H1_res_func(H,alpha, beta,C))

def x_func(L,alpha, beta,C,Psi):
 return L - C*integrate.quad(lambda phi: (special.ellipk(np.sin(phi)**2)*np.sin(phi))/np.sqrt((1 - alpha* np.sin(phi)**2)*(1 - beta* np.sin(phi)**2)), 0, Psi)[0]

def z_func(H,alpha, beta,C,Psi):
 return H + H0_func(alpha, beta,C) + C*integrate.quad(lambda phi: (special.ellipk(np.cos(phi)**2)*np.sin(phi))/np.sqrt((1 - alpha* np.sin(phi)**2)*(1 - beta* np.sin(phi)**2)), 0, Psi)[0]

def x_res_func(L,alpha, beta,C,Psi,x):
    return x - L + C*integrate.quad(lambda phi: (special.ellipk(np.sin(phi)**2)*np.sin(phi))/np.sqrt((1 - alpha* np.sin(phi)**2)*(1 - beta* np.sin(phi)**2)), 0, Psi)[0]#- x_func(L,alpha, beta,C,Psi)

##########################################
#Enter H and L values
H = 0.0001     #lake level
L = 1     #length of reservoir
n = 2001    #number of points for free surface profile
##########################################

res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,10)))

print("======================")
print("Given values")
print("======================")
print("1.) Lake level, H: \t", H)
print("2.) Aquifer length, L: \t", L)
print("======================")
print("Output values")
print("======================")
print("1.) Seepage face height, H0: \t", H0_func(res.x[0],res.x[1],res.x[2])) 
print("2.) alpha: \t", res.x[0])
print("3.) beta: \t", res.x[1])
print("4.) C: \t \t", res.x[2] )

print("======================")
print("Free surface profiles")
print("======================")
print("Psi \t \t \t \t x \t \t \t \t \t z")
'''
x_array = np.linspace(0,L,n)
for x in x_array:
    def full_equations_new(p):
        Psi = p
        return (x_res_func(L,res.x[0], res.x[1],res.x[2],Psi,x))
    
    res_new = least_squares(full_equations_new, (1e-2), bounds = ((0), (1)))

    print(res_new.x[0],x,z_func(H,res.x[0], res.x[1],res.x[2],res_new.x))
'''
Psi_array = np.linspace(0,10,n)
x_array = []
z_array = []
xz_array = []
for i in range(0,n):
    x  = x_func(L,res.x[0], res.x[1],res.x[2],Psi_array[i])
    z  = z_func(H,res.x[0], res.x[1],res.x[2],Psi_array[i])    
    xz_array.append([x,z])
    x_array.append(x)
    if i%(int(n/10)) == 0:    
        print(Psi_array[i],'\t','\t','\t','\t',x,'\t','\t','\t','\t', z)
x_array  = np.array(x_array)
xz_array = np.array(xz_array)
xz_array[x_array<0,:] = np.nan
xz_array[x_array>L,:] = np.nan
#xz_array = xz_array[~np.isnan(xz_array).any(axis=1),:]
#xz_array = xz_array[~np.isinf(xz_array).any(axis=1),:]

fig = plt.figure(figsize=(8,8) , dpi=100)
plt.plot(xz_array[:,0],xz_array[:,1],'b-')
plt.vlines(0,0,1,colors='blue') 
plt.vlines(L,0,H+H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
plt.vlines(L,0,H,colors='blue') 
plt.hlines(H,L,1.1*L,colors='blue')   
plt.hlines(0,0,1.1*L,colors='blue')   

H0 = H0_func(res.x[0],res.x[1],res.x[2]) 

np.savetxt(f"H{H}_L{L}_H0_{H0}.csv", xz_array, delimiter=",")


