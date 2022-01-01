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

#Enter H and L values
H = 0.1
L = 1.0

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

res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,10)))

print("======================")
print("Given values")
print("======================")
print("1.) H : \t", H)
print("2.) L : \t", H)
print("\n")
print("======================")
print("Output values")
print("======================")
print("1.) Seepage face height, H0:", H0_func(res.x[0],res.x[1],res.x[2]),"\n 2.) alpha: \t", res.x[0],"\n 3.) beta: \t", res.x[1],"\n 4.) C: \t", res.x[2] )



