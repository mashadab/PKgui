#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PK GUI Driver Code for Mac
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

#import libraries
import PySimpleGUI as sg
from Pbk_main_function_latest import *

sg.theme('Dark Green 7') #picking theme

#window layout of two columns using PySimpleGUI
file_list_column = [  [sg.Txt('Polubarinova-Kochina Solution',font = ("Serif", 28,"bold"))],
           [sg.Txt('  ',font = ("Serif", 10))],
           [sg.Txt('Developed by:',font = ("Serif", 15,'bold')),sg.Txt('M.A. Shadab*, E. Hiatt & M.A. Hesse',font = ("Serif", 15))],
           [sg.Txt('\t           The University of Texas at Austin',font = ("Serif", 15,'italic'))],
           [sg.Txt('\t           * Email: mashadab@utexas.edu',font = ("Serif", 15))],
           [sg.Txt('',font = ("Serif", 5))],
           [sg.Txt('Notes: - All units must be consistent. Fill only numbers.',font = ("Serif", 15))],
           [sg.Txt('            - For low-aspect ratio dams.',font = ("Serif", 15))],
           [sg.Txt('            - Fill 3 variables when not picking Q & K.',font = ("Serif", 15))],
           [sg.Txt('            - Fill 4 variables when picking Q or K or both.',font = ("Serif", 15))],
           [sg.Txt('_________________________________________________________',font = ("Serif", 18))],
           [sg.Txt('',font = ("Serif", 5))],
           [sg.Txt('Input variables',font = ("Serif", 18,'bold','italic')),
            sg.Txt('(Check relevant boxes)',font = ("Serif", 16,'bold'))],
           [sg.Txt('Units:   Length [ L ]',font = ("Serif", 15)), 
            sg.In(size=(4,1), key='-U-',font = ("Serif", 15)),
            sg.Txt(',',font = ("Serif", 20)),
            sg.Txt('Time [ T ]',font = ("Serif", 15)),
            sg.In(size=(4,1), key='-U2-',font = ("Serif", 15))],
           [sg.Checkbox('Dam length, L [ L ]',font = ("Serif", 15), default=False, key="-CheckL-"), sg.Txt('  \t                 ',font = ("Serif", 10)),sg.In(size=(8,1), key='-L-',font = ("Serif", 15))],
           [sg.Checkbox('Lower lake height, H  [ L ]',font = ("Serif", 15), default=False, key="-CheckH-"),sg.Txt('                   ',font = ("Serif", 10)),sg.In(size=(8,1), key='-H-',font = ("Serif", 15))],
           [sg.Checkbox('Upper lake height, H1  [ L ]',font = ("Serif", 15), default=False, key="-CheckH1-"),sg.Txt('                ',font = ("Serif", 10)),sg.In(size=(8,1), key='-H1-',font = ("Serif", 15))],
           [sg.Checkbox('Flow rate by width, Q_H1 [ L^2/T ]',font = ("Serif", 15), default=False, key="-CheckQ-"),sg.Txt('',font = ("Serif", 10)),sg.In(size=(8,1), key='-Q-',font = ("Serif", 15))],
           [sg.Checkbox('Seepage face height, H0  [ L ]',font = ("Serif", 15), default=False, key="-CheckH0-"),sg.Txt('         ',font = ("Serif", 10)),sg.In(size=(8,1), key='-H0-',font = ("Serif", 15))],
           [sg.Checkbox('Hydraulic conductivity, K  [ L/T ]',font = ("Serif", 15), default=False, key="-CheckK-"),sg.Txt('    ',font = ("Serif", 10)),sg.In(size=(8,1), key='-K-',font = ("Serif", 15))],
           [sg.Txt('Free surface resolution:',font = ("Serif", 15)),
            sg.Radio('Low',"RADIO1",font = ("Serif", 15), default=True, key="-CheckLowRes-"),
            sg.Radio('High',"RADIO1",font = ("Serif", 15), default=False, key="-CheckHighRes-"),
            sg.Radio('Very high',"RADIO1",font = ("Serif", 15), default=False, key="-CheckVHRes-")],
           [sg.Txt(' ',font = ("Serif", 2))], 
           [sg.Button('\t        Calculate \t          ', bind_return_key=True,font = ("Serif", 20,'bold','italic'))],
           [sg.Text("Output Folder  [e.g. /Users/admin/Desktop]",font = ("Serif", 15))],
           [
            sg.In(size=(24,1), key="-FOLDER-",font = ("Serif", 15)),
            sg.FolderBrowse(font = ("Serif", 20)),
            sg.Txt(' ',font = ("Serif", 15)),
            sg.Button('  Save  ', bind_return_key=True,font = ("Serif", 20,'italic','bold'))
            ]
           ]

#output variables and photo
image_viewer_column = [
           [sg.Txt('Output details',font = ("Serif", 18,'bold','italic'))],
           [sg.Txt('   L :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT5-') ,
            sg.Txt('    Q :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT6-')],           
           [sg.Txt('  H :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT4-') ,
            sg.Txt('  H0 :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT-')],
           [sg.Txt('H1 :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT7-') ,
            sg.Txt('    K :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT8-')],
           [sg.Txt('Q/K:',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT9-') ,
            sg.Txt('Q_H0/Q_H:',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT11-')] ,
            [sg.Txt('Alpha :',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT1-') ,
             sg.Txt('Beta :',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT2-')] ,
            [sg.Txt('   C :',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT3-')],
            [sg.Txt('',font = ("Serif", 20)),sg.Txt(size=(50,1), key='-OUTPUT10-')],
           [sg.Txt(' ',font = ("Serif", 2))],
            [sg.Image(size=(550,490),key="-IMAGE-")], 
]

layout = [
    [
         sg.Column(file_list_column),
         sg.VSeparator(),
         sg.Column(image_viewer_column),
     ]    
]

window = sg.Window('Polubarinova-Kochina Solution', layout)

#Running the program
while True:
    event, values = window.read()

    if event != sg.WIN_CLOSED:  
        if values['-CheckH-'] == False: #H: Lower lake height
            H = nan
        else: 
            H = float(values['-H-'])
        if values['-CheckH1-'] == False: #H1: Upper lake height
            H1 = nan                    
        else: 
            H1 = float(values['-H1-'])  
            H1_input = H1
        if values['-CheckH0-'] == False: #H0: Seepage lake height
            H0 = nan
        else: 
            H0 = float(values['-H0-'])
            H0_input = H0
        if values['-CheckL-'] == False:  #L: Length of the dam
            L = nan
        else: 
            L = float(values['-L-'])

        if values['-CheckLowRes-'] == True: #N: Setting resolution for free surface ht
            N = 100
        elif values['-CheckHighRes-'] == True: 
            N = 2000
        else:
            N = 50000
            
        if values['-U-'] == '': #U: length units
            unit = 'L'
        else: 
            unit = str(values['-U-'])

        if values['-U2-'] == '': #U2: time units
            Tunit = 'T'
        else: 
            Tunit = str(values['-U2-'])        

        if values['-CheckQ-'] == False:   #Q: total flow rate  
            Q = nan
        else: 
            Q = float(values['-Q-'])

        if values['-CheckK-'] == False:  #K: hydraulic cond.
            K = nan
        else:
            K = float(values['-K-'])
        
        if values['-FOLDER-'] == '':  #FOLDER: output folder directory
            output_folder = '/tmp'
        else: 
            output_folder = str(values['-FOLDER-'])
                  
        window['-OUTPUT10-'].update('',font = ("Serif", 15)),sg.Txt(size=(50,1))
        window["-IMAGE-"].update('',size=(550,490))

        #indicator of nan
        if (isnan(H) and ~isnan(H1) and isnan(H0) and ~isnan(L) and ~isnan(Q/K)) or (~isnan(H) and isnan(H1) and ~isnan(H0) and isnan(L) and ~isnan(Q/K)):
            buzzer = 1
        else: 
            buzzer = 0

        #solving the equations nonlinearly           
        H0, H, L, res, xz_array,Q,K, H1,QH0byQH,QbyK, output_folder_full = PbK_solution_full(H0,H,L,H1,N,output_folder,Q,K,unit,Tunit)
        calc = f'{H0:.7f} [{unit}]'
        calc1 = f'{res[0]:.7f}'
        calc2 = f'{res[1]:.7f}'        
        calc3 = f'{res[2]:.7f} [{unit}]'
        calc4 = f'{H:.7f} [{unit}]'
        calc5 = f'{L:.7f} [{unit}]'
        calc6 = f'{Q:.7f} [{unit}^2/{Tunit}]'       
        calc7 = f'{H1:.7f} [{unit}]'
        calc8 = f'{K:.7f} [{unit}/{Tunit}]'
        calc9 = f'{QbyK:.7f} [{unit}]'
        calc10= f'{QH0byQH:.7f}'
        
        print('Worked')
        
        #Showing output on a case by case basis
        if isnan(H0) and isnan(H1) and isnan(res):
            print('Did not work')
            calc = 'Invalid H0'
            calc1 = 'Invalid Alpha'
            calc2 = 'Invalid Beta'
            calc3 = 'Invalid C'
            calc4 = 'Invalid H'
            calc5 = 'Invalid L'
            calc6 = 'Invalid Q'
            calc7 = 'Invalid H1'
            calc8 = 'Invalid K'
            calc9 = 'Invalid QbyK'
            calc10= 'Invalid Q_H0/Q_H'

        window['-OUTPUT-'].update(calc,font = ("Serif", 15))
        window['-OUTPUT1-'].update(calc1,font = ("Serif", 15))
        window['-OUTPUT2-'].update(calc2,font = ("Serif", 15))
        window['-OUTPUT3-'].update(calc3,font = ("Serif", 15))
        window['-OUTPUT4-'].update(calc4,font = ("Serif", 15))
        window['-OUTPUT5-'].update(calc5,font = ("Serif", 15))
        window['-OUTPUT6-'].update(calc6,font = ("Serif", 15))
        window['-OUTPUT7-'].update(calc7,font = ("Serif", 15))
        window['-OUTPUT8-'].update(calc8,font = ("Serif", 15))
        window['-OUTPUT9-'].update(calc9,font = ("Serif", 15))
        window['-OUTPUT11-'].update(calc10,font = ("Serif", 15))   
        
        try:
            H1_input
        except NameError:
            H1_input = H1
            
            print('L/H1 ratio is',L/H1_input)

        if L/H1_input>=3.5 or (H1_input**2-H**2)/L**2<0.1: #Checking if the PK solution is valid
        
            if res[0]<1e-3 and res[1]>1-1e-3:
                if not N ==50000:
                    filename = f"Caution: Limiting case!"  
                    window['-OUTPUT10-'].update(filename,font = ("Serif", 15),text_color='Red')
                filename = f"{output_folder_full}/free-surface-profile.png"
                window["-IMAGE-"].update(filename=filename)

            else: #Checking if the PK solution is valid: Throw an error
                filename = f"Error: The aspect ratio is high!"
                window['-OUTPUT10-'].update(filename,font = ("Serif", 20),text_color='Red')
            
        else: 
            if not N ==50000: #Checking if higher resolution is needed 
                filename = f"Caution: Try higher resolution if free surface is disconnected!"
                window['-OUTPUT10-'].update(filename,font = ("Serif", 15),text_color='Red')
               

            if buzzer==1: #Checking if accuracy is low
                print('Hello there')
                filename = f"Caution: Accuracy may be low! Pick other variables."
                window['-OUTPUT10-'].update(filename,font = ("Serif", 15),text_color='Red')                
                
            filename = f"{output_folder_full}/free-surface-profile.png" #Saving the image
            window["-IMAGE-"].update(filename=filename)
        
    else:
        break

window.close()