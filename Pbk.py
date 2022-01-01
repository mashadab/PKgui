#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk solution
"""

from scipy.integrate import quad
from scipy.special import ellipk
from numpy import sin, cos, sqrt, array, pi, nan, inf, savetxt, linspace, isnan
from scipy.optimize import fsolve,least_squares
from matplotlib.pyplot import figure, xlabel, ylabel, plot, vlines, hlines, savefig

def L_res_func(L,alpha, beta, C):
 return L - C*quad(lambda phi: ellipk(alpha + (beta - alpha)* sin(phi)**2)/sqrt(1 - alpha - (beta - alpha)* sin(phi)**2), 0, pi/2)[0]

def H_res_func(H,alpha, beta,C):
 return H - C*sqrt(alpha)*quad(lambda phi: (ellipk(alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

def H0_func(alpha, beta,C):
 return C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0]

def H1_res_func(H,alpha, beta,C):
 return 1 - H - H0_func(alpha, beta,C) - C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, pi/2)[0]

def full_equations(p):
    alpha, beta, C = p
    return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),H1_res_func(H,alpha, beta,C))

def x_func(L,alpha, beta,C,Psi):
 return L - C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

def z_func(H,alpha, beta,C,Psi):
 return H + H0_func(alpha, beta,C) + C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

def x_res_func(L,alpha, beta,C,Psi,x):
    return x - L + C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]#- x_func(L,alpha, beta,C,Psi)

##########################################
#Enter H and L values
H = 0.5     #lake level
L = 0.5     #length of reservoir
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
print("Psi \t \t x \t \t z")
'''
x_array = linspace(0,L,n)
for x in x_array:
    def full_equations_new(p):
        Psi = p
        return (x_res_func(L,res.x[0], res.x[1],res.x[2],Psi,x))
    
    res_new = least_squares(full_equations_new, (1e-2), bounds = ((0), (1)))

    print(res_new.x[0],x,z_func(H,res.x[0], res.x[1],res.x[2],res_new.x))
'''
Psi_array = linspace(0,10,n)
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
    if i%(int(n/20)) == 0 and x != inf and x!=-inf:    
        print(Psi_array[i],'\t',x,'\t','\t','\t','\t', z)
x_array  = array(x_array)
xz_array = array(xz_array)
xz_array[x_array<0,:] = nan
xz_array[x_array>L,:] = nan
#xz_array = xz_array[~isnan(xz_array).any(axis=1)]

fig = figure(figsize=(8,8) , dpi=100)
plot(xz_array[:,0],xz_array[:,1],'b-')
vlines(0,0,1,colors='blue') 
vlines(L,0,H+H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
vlines(L,0,H,colors='blue') 
hlines(H,L,1.1*L,colors='blue')   
hlines(0,0,1.1*L,colors='blue')   
xlabel('x')
ylabel('z')
H0 = H0_func(res.x[0],res.x[1],res.x[2]) 


savefig(f"H{H}_L{L}_H0_{H0}.pdf")
savetxt(f"H{H}_L{L}_H0_{H0}.csv", xz_array, delimiter=",")


