#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PK Non-linear equation solution 
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
License: MIT

_______________________________________________
Main Variables

#H: Lower lake height
#H1: Upper lake height
#H0: Seepage lake height
#L: Length of the dam
#N: Setting resolution for free surface height
#U: length units
#U2: time units
#Q: total flow rate 
#K: hydraulic cond.
#FOLDER: output folder directory
_______________________________________________

"""

#importing required libraries
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

#setting default font and ignore warnings
rcParams.update({'font.size': 22})
rcParams.update({'font.family': 'Serif'})
filterwarnings("ignore")


##defining variables as functions (_func) and their residual (_res_func) from the nonlinear functions given in Appendix of the manuscript

#length residual
def L_res_func(L,alpha, beta, C):
 return L - C*quad(lambda phi: ellipk(alpha + (beta - alpha)* sin(phi)**2)/sqrt(1 - alpha - (beta - alpha)* sin(phi)**2), 0, pi/2)[0]

#length 
def L_func(alpha, beta, C):
 return C*quad(lambda phi: ellipk(alpha + (beta - alpha)* sin(phi)**2)/sqrt(1 - alpha - (beta - alpha)* sin(phi)**2), 0, pi/2)[0]

#lower lake residual
def H_res_func(H,alpha, beta,C):
 return H - C*sqrt(alpha)*quad(lambda phi: (ellipk(alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

#lower lake 
def H_func(alpha, beta,C):
 return C*sqrt(alpha)*quad(lambda phi: (ellipk(alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

#seepage face residual
def H0_res_func(H0,alpha, beta,C):
 return H0 - C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0]

#seepage face residual
def H0_func(alpha, beta,C):
 return C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0]

#upper lake residual
def H1_res_func(H1,H,alpha, beta,C):
 return H1 - H - H0_func(alpha, beta,C) - C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, pi/2)[0]

#upper lake
def H1_func(H,alpha, beta,C):
 return H + H0_func(alpha, beta,C) + C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, pi/2)[0]

#horizontal distance from seepage face
def x_func(L,alpha, beta,C,Psi):
 return L - C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

#vertical height of free surface
def z_func(H,alpha, beta,C,Psi):
 return H + H0_func(alpha, beta,C) + C*quad(lambda phi: (ellipk(cos(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

#horizontal distance from seepage face residual
def x_res_func(L,alpha, beta,C,Psi,x):
    return x - L + C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(1 - beta* sin(phi)**2)), 0, Psi)[0]

#Ratio of total flow rate to hydraulic conductivity
def QbyK_func(alpha, beta,C):
 return C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0] \
  + C*sqrt(alpha)*quad(lambda phi: (ellipk(1 - alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

#Ratio of total flow rate to hydraulic conductivity residual
def QbyK_res_func(QbyK,alpha, beta,C):
 return QbyK - C*quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0]\
 - C*sqrt(alpha)*quad(lambda phi: (ellipk(1 - alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0]

#Fraction of fluid flowing out through seepage face
def QH0byQH_func(alpha, beta):
 return quad(lambda phi: (ellipk(sin(phi)**2)*sin(phi)*cos(phi))/sqrt((1 - (1-alpha)* sin(phi)**2)*(1 - (1-beta)* sin(phi)**2)), 0, pi/2)[0] / (sqrt(alpha)*quad(lambda phi: (ellipk(1 - alpha * sin(phi)**2)*sin(phi))/sqrt((1 - alpha* sin(phi)**2)*(beta - alpha* sin(phi)**2)), 0, pi/2)[0])

#horizontal distance from seepage face residual for high aspect ratio
def x_func_high_aspect_ratio(L,C,Psi): #Psi is m in PbK book
 return  L -0.5*C*quad(lambda m: (ellipk(m))/(1-m), 0, Psi)[0]

#vertical elevation for high aspect ratio
def z_func_high_aspect_ratio(H0,C,Psi): #Psi is m in PbK book
 return H0+0.5*C*quad(lambda m: (ellipk(1-m))/(1-m), 0, Psi)[0]

##main solution function 
def PbK_solution_full(H0,H_full,L_full,H1,n,output_folder,Q,K,unit,Tunit):
    ##Different cases depending on the number of inputs
    
    if ~isnan(H_full) and ~isnan(H1) and isnan(H0) and ~isnan(L_full): #H1, H, L given; _full means dimensional; H0 not given
        
        #Scaling wrt upper lake
        H_scale = H1
        H = H_full/H_scale
        L = L_full/H_scale
        H1= H1/H_scale
        
        #Making a function of nonlinear system of equations to be solved
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),H1_res_func(H1,H,alpha, beta,C))
        a = max(H,L)
        
        #Solving the equations by least squares
        res = least_squares(full_equations, (0.0001, 0.1,0.1*a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        
        #Printing output in terminal
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
        
        #Solving the equations by least squares
        if isnan(Q) and ~isnan(K): 
            Q = K*H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])
            print('hahaprint Q', Q,'hahaprint QbyK',QbyK_func(res.x[0],res.x[1],res.x[2]),'hahaprint Hscale',H_scale)
            QbyK = H_scale*QbyK_func(res.x[0],res.x[1],res.x[2]) 
            print("5.) Specific discharge, Q: \t \t", Q, f'{unit}^2/{Tunit}') 
        elif ~isnan(Q) and isnan(K): 
            K = Q/(H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])) 
            QbyK = H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])
            print("5.) Hydraulic conductivity, K: \t \t", K ,f'{unit}/{Tunit}') 
        elif isnan(Q) and isnan(K): 
            QbyK = H_scale*QbyK_func(res.x[0],res.x[1],res.x[2])     
            print("5.) Neither K or Q given") 

    elif ~isnan(H_full) and isnan(H1) and ~isnan(H0) and ~isnan(L_full):  #H0, H, L given; H1 not given
    
        #Scaling again
        H_scale = L_full
        H = H_full/H_scale
        L = L_full/H_scale
        H0 = H0/H_scale
        
        #Making a function of nonlinear system of equations to be solved
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),H0_res_func(H0,alpha, beta,C))
        a = max(H,L)
        
        #Solving the equations by least squares
        res = least_squares(full_equations, (0.0001, 0.1,0.1*a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        
        #Printing output in terminal
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
            
            
    elif ~isnan(H_full) and ~isnan(H1) and ~isnan(H0) and isnan(L_full): #H0, H, H1 given; L not given
        H_scale = H1
        H = H_full/H_scale
        H0 = H0/H_scale
        H1= H1/H_scale

        #Making a function of nonlinear system of equations to be solved
        def full_equations(p):
            alpha, beta, C = p
            return (H1_res_func(H1,H,alpha, beta,C), H_res_func(H,alpha, beta,C),H0_res_func(H0,alpha, beta,C))
        a = max(H1,H0)
        
        #Solving the equations by least squares
        res = least_squares(full_equations, (0.0001, 0.1,a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        
        L = L_func(res.x[0],res.x[1],res.x[2]) 
        L_full = L*H_scale
        
        #Printing output in terminal
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
 
    elif isnan(H_full) and ~isnan(H1) and ~isnan(H0) and ~isnan(L_full): #H0, L, H1 given; H not given
        H_scale = H1
        L = L_full/H_scale
        H0 = H0/H_scale
        H1= H1/H_scale
        def full_equations(p):
            alpha, beta, C = p
            return (H1_res_func(H1,H_func(alpha, beta,C),alpha, beta,C), L_res_func(L,alpha, beta,C),H0_res_func(H0,alpha, beta,C))
        a = max(H1,L)
        
        res = least_squares(full_equations, (0.0001, 0.1,0.1*a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        
        H   = H_func(res.x[0],res.x[1],res.x[2])  
        H_full = H_scale * H
        
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

    ##Adding the cases with given Q/K now            
    elif isnan(H_full) and isnan(H1) and ~isnan(H0) and ~isnan(L_full) and ~isnan(Q/K): #H0, L, Q/K given; H1, H not given     
        H_scale = L_full
        L = L_full/H_scale
        H0 = H0/H_scale
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H0_res_func(H0,alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        a = max(H0,L)
        res = least_squares(full_equations, (0.0001, 0.1,a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        H   = H_func(res.x[0],res.x[1],res.x[2])
        H_full= H_func(res.x[0],res.x[1],res.x[2])*H_scale
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Lower lake level, H: \t", H_scale*H_func(res.x[0],res.x[1],res.x[2]),f"{unit}")   
    
        
    elif isnan(H_full) and ~isnan(H1) and isnan(H0) and ~isnan(L_full) and ~isnan(Q/K):  #H1, L, Q/K given; H0, H not given       
        H_scale = H1
        H1 = H1/H_scale
        L = L_full/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H1_res_func(H1,H0_func(alpha, beta,C),alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        a = max(H1,L,Q/K)
        
        res = least_squares(full_equations, (0.0001, 0.1,0.1*a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        H   = H_func(res.x[0],res.x[1],res.x[2])
        H_full= H_func(res.x[0],res.x[1],res.x[2])*H_scale
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Lower lake level, H: \t", H_scale*H_func(res.x[0],res.x[1],res.x[2]),f"{unit}")   
    

    elif isnan(H_full) and ~isnan(H1) and ~isnan(H0) and isnan(L_full) and ~isnan(Q/K):   #H1, H0, Q/K given; L, H not given        
        H_scale = H1
        H1 = H1/H_scale
        H0 = H0/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (H0_res_func(H0,alpha, beta,C), H1_res_func(H1,H_func(alpha, beta,C),alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        #Calculate other functions
        res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,inf)))
        H   = H_func(res.x[0],res.x[1],res.x[2])
        H_full= H_func(res.x[0],res.x[1],res.x[2])*H_scale
        L   = L_func(res.x[0],res.x[1],res.x[2])
        L_full = L_func(res.x[0],res.x[1],res.x[2])*H_scale
        
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2]),ftol=1e-12) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Lower lake level, H: \t", H_scale*H_func(res.x[0],res.x[1],res.x[2]),f"{unit}")   

    
    elif ~isnan(H_full) and isnan(H1) and isnan(H0) and ~isnan(L_full) and ~isnan(Q/K):   #H, L, Q/K given; H0, H1 not given 
        H_scale = L_full
        H = H_full/H_scale
        L = L_full/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (L_res_func(L,alpha, beta,C), H_res_func(H,alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-12)
        
        L   = L_func(res.x[0],res.x[1],res.x[2])
        L_full = L_func(res.x[0],res.x[1],res.x[2])*H_scale    
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Higher lake level, H1: \t", H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]),f"{unit}")   
        
 
    elif ~isnan(H_full) and isnan(H1) and ~isnan(H0) and isnan(L_full) and ~isnan(Q/K):
        if H_full != 0 and H0 != 0:
            H_scale = max(H_full,H0)
        elif H_full == 0 and H0 != 0:
            H_scale = H0
        elif H_full != 0 and H0 == 0:
            H_scale = H_full
        else:
            H_scale = Q/K
            
        H = H_full/H_scale
        H0= H0/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (H0_res_func(H0,alpha, beta,C), H_res_func(H,alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        a = min(H0,H_full)
        res = least_squares(full_equations, (0.0001, 0.1,a), bounds = ((0, 0,0), (1,1,inf)),ftol=1e-14)
        L   = L_func(res.x[0],res.x[1],res.x[2])
        L_full = L_func(res.x[0],res.x[1],res.x[2])*H_scale
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Higher lake level, H1: \t", H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]),f"{unit}")   
        
    elif ~isnan(H_full) and ~isnan(H1) and isnan(H0) and isnan(L_full) and ~isnan(Q/K):   #H1, H, Q/K given; H0, L not given 
        H_scale = H1
        H = H_full/H_scale
        H1= H1/H_scale
        
        def full_equations(p):
            alpha, beta, C = p
            return (H1_res_func(H1,H_func(alpha, beta,C),alpha, beta,C), H_res_func(H,alpha, beta,C),QbyK_res_func(Q/(K*H_scale),alpha, beta,C))
        
        res = least_squares(full_equations, (0.0001, 0.1,1), bounds = ((0, 0,0), (1,1,inf)))

        L   = L_func(res.x[0],res.x[1],res.x[2])
        L_full = L_func(res.x[0],res.x[1],res.x[2])*H_scale
        
        print("======================")
        print("Output values")
        print("======================")
        print("1.) Seepage face height, H0: \t", H_scale*H0_func(res.x[0],res.x[1],res.x[2])) 
        print("2.) alpha: \t", res.x[0])
        print("3.) beta: \t", res.x[1])
        print("4.) C: \t \t", res.x[2] )  
        print("5.) Higher lake level, H1: \t", H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]),f"{unit}")   
        
       
    print(f"6.) Specific discharge over hydraulic conductivity, q/K {unit}: \t \t", H_scale*QbyK_func(res.x[0], res.x[1],res.x[2]))          
    
    ##Calculating the free surface profile
    print("======================")
    print("Free surface profiles")
    print("======================")
    print(f"Psi \t x {unit} \t z {unit}")
    
    
    ##Making an array of parameter Psi ref to Equation 13 and 14 of the MS  
    Psi_array = linspace(0,10,n+1)
    
    if res.x[0]<1e-3 and res.x[1]>1-1e-3:
        
        #The least square rarely works so we do it directly from line 471 onwards
        #horizontal length residual in terms of parameter psi (Eq 13)
        def x_func_high_aspect_ratio_res(Psi): #Psi is m in PbK book
         return  L -0.5*res.x[2]*quad(lambda m: (ellipk(m))/(1-m), 0, Psi)[0]       
        print(L,res.x[2])
        
        Psi_max = least_squares(x_func_high_aspect_ratio_res, (1e-5), bounds = (-inf,inf))
        print('Psi max is',Psi_max.x,x_func_high_aspect_ratio_res(Psi_max.x))

        #vertical elevation of GW residual in terms of parameter psi (Eq 14)
        def x_func_high_aspect_ratio_res(Psi): #Psi is m in PbK book
         return  0.5*res.x[2]*quad(lambda m: (ellipk(m))/(1-m), 0, Psi)[0]

        Psi_min = least_squares(x_func_high_aspect_ratio_res, (1e-5), bounds = (-inf,inf))
        print('Psi min is',Psi_min.x,x_func_high_aspect_ratio_res(Psi_min.x))
        
        Psi_array = linspace(0,1,n)
    
    #Calculating x and z from a predefined Psi
    x_array = []
    z_array = []
    xz_array = []       
    x_newarray = linspace(0,L,n+1)
    for i in range(0,n):
        if res.x[0]<1e-3 and res.x[1]>1-1e-3:
            x  = x_func_high_aspect_ratio(L,res.x[2],Psi_array[i])
            z  = z_func_high_aspect_ratio(H0_func(res.x[0],res.x[1],res.x[2]),res.x[2],Psi_array[i])         
        else:
            x  = x_func(L,res.x[0], res.x[1],res.x[2],Psi_array[i])
            z  = z_func(H,res.x[0], res.x[1],res.x[2],Psi_array[i]) 
        if x==0: z = H1_func(H,res.x[0],res.x[1],res.x[2])
        if x==L: z = H + H0_func(res.x[0],res.x[1],res.x[2]) 
        xz_array.append([x,z])
        x_array.append(x)
        if i%(int(n/20)) == 0 and x != inf and x!=-inf and x>=0:    
            print(Psi_array[i],'\t',H_scale*x,'\t','\t','\t','\t', H_scale*z)
            
    #Omitting out the x values which go less than zero and beyond L  
    x_array  = array(x_array)
    xz_array = array(xz_array)
    xz_array[x_array<0,:] = nan
    xz_array[x_array>L,:] = nan
    xz_array = H_scale*xz_array
    xz_array = xz_array[~isnan(xz_array).any(axis=1),:] 
    
    p = argsort(xz_array[:,0])
    
    xz_array = xz_array[p,:]
    
    #Plotting the free surface along with the lakes
    fig = figure(figsize=(6,6) , dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(xz_array[:,0],xz_array[:,1],'b-')
    ax.vlines(0,0,H_scale*H1_func(H_func(res.x[0],res.x[1],res.x[2]),res.x[0],res.x[1],res.x[2]),colors='blue') 
    ax.vlines(H_scale*L,0,H_scale*H_func(res.x[0],res.x[1],res.x[2])+H_scale*H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
    ax.vlines(H_scale*L,0,H_scale*H_func(res.x[0],res.x[1],res.x[2]),colors='blue') 
    ax.hlines(H_scale*H_func(res.x[0],res.x[1],res.x[2]),H_scale*L,H_scale*1.1*L,colors='blue')   
    ax.hlines(0,0,H_scale*1.1*L,colors='blue')   
    ax.set_xlabel(f'x [{unit}]')
    ax.set_ylabel(f'z [{unit}]')
    tight_layout(pad=1, w_pad=0.8, h_pad=1)
    #fig.show()

    H0 = H_scale*H0_func(res.x[0],res.x[1],res.x[2]) 
    H1 = H_scale*H1_func(H,res.x[0],res.x[1],res.x[2])

    #Output folder location
    path = f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}" #path of the folder
    
    #making a subdirectory if it does not exist already
    if not exists(path):
        makedirs(path)
    else:
        rmtree(path)           # Removes all the subdirectories!
        makedirs(path) 
    
    if not output_folder =='/tmp': #save figure as pdf
        fig.savefig(f"{output_folder}/L{L_full}{unit}_H{H_full}{unit}_H1_{H1}{unit}_N{n}/free-surface-profile.pdf")
    
    #plotting as png to show
    fig = figure(figsize=(10,10), dpi=50)
    ax = fig.add_subplot(111)
    ax.plot(xz_array[:,0],xz_array[:,1],'b-')
    ax.vlines(0,0,H_scale*H1_func(H_func(res.x[0],res.x[1],res.x[2]),res.x[0],res.x[1],res.x[2]),colors='blue') 
    ax.vlines(H_scale*L,0,H_scale*H_func(res.x[0],res.x[1],res.x[2])+H_scale*H0_func(res.x[0],res.x[1],res.x[2]),colors='blue')  
    ax.vlines(H_scale*L,0,H_scale*H_func(res.x[0],res.x[1],res.x[2]),colors='blue') 
    ax.hlines(H_scale*H_func(res.x[0],res.x[1],res.x[2]),H_scale*L,H_scale*1.1*L,colors='blue')   
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

    return H0, H_full, L_full, res.x, xz_array,Q,K, H_scale*H1_func(H,res.x[0],res.x[1],res.x[2]/H_scale),QH0byQH,QbyK*H_scale,output_folder

