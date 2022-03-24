#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk GUI
"""

import PySimpleGUI as sg
from Pbk_main_function_latest import *

sg.theme('Dark Green 7')


#window layout of two columns
file_list_column = [  [sg.Txt('Polubarinova-Kochina Solution',font = ("Serif", 28,"bold"))],
           [sg.Txt('  ',font = ("Serif", 10))],
           [sg.Txt('Developed by:',font = ("Serif", 15,'bold')),sg.Txt('M.A. Shadab*, E. Hiatt & M.A. Hesse',font = ("Serif", 15))],
           [sg.Txt('The University of Texas at Austin    \ ',font = ("Serif", 15,'italic')),sg.Txt('* mashadab@utexas.edu',font = ("Serif", 15))],
           [sg.Txt('',font = ("Serif", 2))],
           [sg.Txt('Notes: - All units must be consistent.',font = ("Serif", 15))],
           [sg.Txt('            - For low-aspect ratio dams.',font = ("Serif", 15))],
           [sg.Txt('            - Fill 3 variables when not picking Q & K.',font = ("Serif", 15))],
           [sg.Txt('            - Fill 4 variables when picking Q or K or both.',font = ("Serif", 15))],
           [sg.Txt('_________________________________________________________',font = ("Serif", 18))],
           [sg.Txt('Input variables',font = ("Serif", 18,'bold','italic')),
            sg.Txt('(Check relevant boxes)',font = ("Serif", 16,'bold'))],
           [sg.Txt('Units:   Length [e.g. m]',font = ("Serif", 15)), 
            sg.In(size=(4,1), key='-U-',font = ("Serif", 15)),
            sg.Txt(',',font = ("Serif", 20)),
            sg.Txt('Time [e.g. day]',font = ("Serif", 15)),
            sg.In(size=(4,1), key='-U2-',font = ("Serif", 15))],
           [sg.Checkbox('Dam length, L [e.g. 100]',font = ("Serif", 15), default=False, key="-CheckL-"), sg.Txt('                     ',font = ("Serif", 10)),sg.In(size=(8,1), key='-L-',font = ("Serif", 15))],
           [sg.Checkbox('Lower lake height, H  [e.g. 10]',font = ("Serif", 15), default=False, key="-CheckH-"),sg.Txt('        ',font = ("Serif", 10)),sg.In(size=(8,1), key='-H-',font = ("Serif", 15))],
           [sg.Checkbox('Upper lake height, H1  [e.g. 110]',font = ("Serif", 15), default=False, key="-CheckH1-"),sg.Txt(' ',font = ("Serif", 10)),sg.In(size=(8,1), key='-H1-',font = ("Serif", 15))],
           [sg.Checkbox('Specific discharge, Q  [e.g. 1]',font = ("Serif", 15), default=False, key="-CheckQ-"),sg.Txt('         ',font = ("Serif", 10)),sg.In(size=(8,1), key='-Q-',font = ("Serif", 15))],
           [sg.Checkbox('Seepage face height, H0  [e.g. 1]',font = ("Serif", 15), default=False, key="-CheckH0-"),sg.Txt(' ',font = ("Serif", 10)),sg.In(size=(8,1), key='-H0-',font = ("Serif", 15))],
           [sg.Checkbox('Hydraulic conductivity, K  [e.g. 2]',font = ("Serif", 15), default=False, key="-CheckK-"),sg.Txt('',font = ("Serif", 10)),sg.In(size=(8,1), key='-K-',font = ("Serif", 15))],
           [sg.Txt('Free surface resolution:',font = ("Serif", 15)),
            sg.Radio('Low',"RADIO1",font = ("Serif", 15), default=True, key="-CheckLowRes-"),
            sg.Radio('High',"RADIO1",font = ("Serif", 15), default=False, key="-CheckHighRes-"),
            sg.Radio('Very high',"RADIO1",font = ("Serif", 15), default=False, key="-CheckVHRes-")],
           [sg.Txt(' ',font = ("Serif", 2))], 
           [sg.Button('\t   Calculate \t \t ', bind_return_key=True,font = ("Serif", 20,'bold','italic'))],
           [sg.Text("Output Folder  [e.g. /Users/admin/Desktop]",font = ("Serif", 15))],
           [
            sg.In(size=(24,1), key="-FOLDER-",font = ("Serif", 15)),
            sg.FolderBrowse(font = ("Serif", 20)),
            sg.Txt(' ',font = ("Serif", 15)),
            sg.Button('  Save  ', bind_return_key=True,font = ("Serif", 20,'italic','bold'))
            ]
           ]

#for now will only show the name of the chosen file
image_viewer_column = [
           [sg.Txt('Output details',font = ("Serif", 18,'bold','italic'))],
           [sg.Text("Location:",font = ("Serif", 15)),
            sg.Text(size=(50,1), key="-TOUT-",font = ("Serif", 15))],
           [sg.Txt('   L :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT5-') ,
            sg.Txt('    Q :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT6-')],           
           [sg.Txt('  H :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT4-') ,
            sg.Txt('  H0 :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT-')],
           [sg.Txt('H1 :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT7-') ,
            sg.Txt('    K :',font = ("Serif", 15)), sg.Txt(size=(20,1), key='-OUTPUT8-')],
           [sg.Txt('Q/K:',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT9-') ,
            sg.Txt('Alpha :',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT1-')  ],
           [sg.Txt('Beta :',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT2-') ,
            sg.Txt('   C :',font = ("Serif", 15)),sg.Txt(size=(20,1), key='-OUTPUT3-')  ],
             [sg.Txt('',font = ("Serif", 15)),sg.Txt(size=(40,1), key='-OUTPUT10-')],
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

while True:
    event, values = window.read()

    if event != sg.WIN_CLOSED:  
        if values['-CheckH-'] == False:
            H = nan
        else: 
            H = float(values['-H-'])
        if values['-CheckH1-'] == False:
            H1 = nan
        else: 
            H1 = float(values['-H1-'])
            H1_input = H1
        if values['-CheckL-'] == False:
            L = nan
        else: 
            L = float(values['-L-'])

        if values['-CheckLowRes-'] == True:
            N = 100
        elif values['-CheckHighRes-'] == True: 
            N = 1000
        else:
            N = 5000
            
        if values['-U-'] == False:
            unit = 'Length-unit'
        else: 
            unit = str(values['-U-'])

        if values['-U2-'] == '':
            Tunit = 'Time-unit'
        else: 
            Tunit = str(values['-U2-'])        

        if values['-CheckQ-'] == False:
            Q = nan
        else: 
            Q = float(values['-Q-'])

        if values['-CheckK-'] == False:
            K = nan
        else:
            K = float(values['-K-'])
        
        if values['-FOLDER-'] == '':
            output_folder = '/tmp'
        else: 
            output_folder = str(values['-FOLDER-'])
                  
        window['-OUTPUT10-'].update('',font = ("Serif", 15)),sg.Txt(size=(0,1))
        window["-IMAGE-"].update('',size=(550,490))


        H0, H, L, res, xz_array,Q,K, H1,QbyK = PbK_solution_full(nan,H,L,H1,N,output_folder,Q,K,unit,Tunit)
        calc = f'{H0} [{unit}]'
        calc1 = res[0]
        calc2 = res[1]        
        calc3 = f'{res[2]} [{unit}]'
        calc4 = f'{H} [{unit}]'
        calc5 = f'{L} [{unit}]'
        calc6 = f'{Q} [{unit}/{Tunit}]'       
        calc7 = f'{H1} [{unit}]'
        calc8 = f'{K} [{unit}/{Tunit}]'
        calc9 = f'{QbyK}'
        
        print('Worked')
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
        
        try:
            H1_input
        except NameError:
            H1_input = H1
            
            print('L/H1 ratio is',L/H1_input)

        if L/H1_input>3.5:# or (H1_input**2-H**2)/L**2<0.1:
            filename = f"Error: The aspect ratio is low!"
            window['-OUTPUT10-'].update(filename,font = ("Serif", 15))
            
        else:
            filename = f"{output_folder}/L{L}{unit}_H{H}{unit}_H1_{H1}{unit}_N{N}/free-surface-profile.png"
            window["-IMAGE-"].update(filename=filename)
            
            if not output_folder =='/tmp':
                filename = f"{output_folder}/L{L}{unit}_H{H}{unit}_H1_{H1}{unit}_N{N}"
                window["-TOUT-"].update(filename,font = ("Serif", 15))
        
    else:
        break

window.close()