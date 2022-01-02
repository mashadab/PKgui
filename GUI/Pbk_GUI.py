#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk GUI
"""

import PySimpleGUI as sg
from Pbk_main_function import *

sg.theme('Dark Green 7')


#window layout of two columns
file_list_column = [  [sg.Txt('Polubarinova-Kochina Solutions',font = ("Serif", 22,"bold"))],
           [sg.Txt('===========Input variables ============',font = ("Serif", 18))],
           [sg.Txt('Length of the aquifer, L',font = ("Serif", 13))],
           [sg.In(size=(8,1), key='-L-',font = ("Serif", 13))],
           [sg.Txt('Lake level, H',font = ("Serif", 13))],
           [sg.In(size=(8,1), key='-H-',font = ("Serif", 13))],
           [sg.Txt('Rough number of points at free-surface, N [e.g. 1001]',font = ("Serif", 13))],
           [sg.In(size=(8,1), key='-N-',font = ("Serif", 13))],
           [sg.Text("Output Folder",font = ("Serif", 13))],
           [
            sg.In(size=(30,1), key="-FOLDER-",font = ("Serif", 13)),
            sg.FolderBrowse(font = ("Serif", 13)),
            ],
           [sg.Txt('',font = ("Serif", 5))], 
           [sg.Button('Calculate', bind_return_key=True,font = ("Serif", 15))],
           [sg.Txt('===========Output variables ===========',font = ("Serif", 18))],
           [sg.Txt('Seepage face height, H0 :',font = ("Serif", 13))],
           [sg.Txt(size=(30,1), key='-OUTPUT-')  ],
           [sg.Txt('Alpha :',font = ("Serif", 13))],
           [sg.Txt(size=(30,1), key='-OUTPUT1-')  ],
           [sg.Txt('Beta :',font = ("Serif", 13))],
           [sg.Txt(size=(30,1), key='-OUTPUT2-')  ],
           [sg.Txt('C :',font = ("Serif", 13))],
           [sg.Txt(size=(30,1), key='-OUTPUT3-')  ],
           [sg.Txt('===============About Us ===============',font = ("Serif", 18))],
           [sg.Txt('M.A. Shadab, E. Hiatt & M.A. Hesse',font = ("Serif", 15))],
           [sg.Txt('The University of Texas at Austin',font = ("Serif", 15))],
           [sg.Txt('Contact: mashadab@utexas.edu',font = ("Arial", 15))]]

#for now will only show the name of the chosen file
image_viewer_column = [
    
    [sg.Text("Output will be stored as .csv and .pdf in the same folder as exe file. Example of png filename:",font = ("Serif", 13))],    
    [sg.Text(size=(50,1), key="-TOUT-")],
    [sg.Image(size=(600,600),key="-IMAGE-")],
]

layout = [
    [
         sg.Column(file_list_column),
         sg.VSeparator(),
         sg.Column(image_viewer_column),
     ]    
]

window = sg.Window('Polubarinova-Kochina Solutions', layout)

while True:
    event, values = window.read()

    if event != sg.WIN_CLOSED:
        try:
            H = float(values['-H-'])
            L = float(values['-L-'])
            N = int(values['-N-'])
            output_folder = str(values['-FOLDER-'])            
            H0, res, xz_array = PbK_solution(H,L,N,output_folder)
            calc = H0
            calc1 = res[0]
            calc2 = res[1]        
            calc3 = res[2]
        except:
            calc = 'Invalid H0'
            calc1 = 'Invalid Alpha'
            calc2 = 'Invalid Beta'
            calc3 = 'Invalid C'

        window['-OUTPUT-'].update(calc,font = ("Serif", 13))
        window['-OUTPUT1-'].update(calc1,font = ("Serif", 13))
        window['-OUTPUT2-'].update(calc2,font = ("Serif", 13))
        window['-OUTPUT3-'].update(calc3,font = ("Serif", 13))
        
        filename = f"{output_folder}/H{H}_L{L}_H0_{calc}.png"
        window["-TOUT-"].update(filename,font = ("Serif", 13))
        window["-IMAGE-"].update(filename=filename)
    else:
        break