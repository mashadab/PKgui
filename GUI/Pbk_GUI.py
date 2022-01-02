#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk GUI
"""

import PySimpleGUI as sg
from Pbk_main_function import *

sg.theme('Dark Green 7')


#window layout of two columns
file_list_column = [ [sg.Txt('===========Input variables ============',font = ("Serif", 13))],
           [sg.Txt('Length of the aquifer, L',font = ("Serif", 13))],
           [sg.In(size=(8,4), key='-L-',font = ("Serif", 13))],
           [sg.Txt('Lake level, H',font = ("Serif", 13))],
           [sg.In(size=(8,4), key='-H-',font = ("Serif", 13))],
           [sg.Txt('===========Output variables ===========',font = ("Serif", 13))],
           [sg.Txt('Seepage face height, H0 :',font = ("Serif", 13))],
           [sg.Txt(size=(30,4), key='-OUTPUT-')  ],
           [sg.Txt('Alpha :',font = ("Serif", 13))],
           [sg.Txt(size=(30,4), key='-OUTPUT1-')  ],
           [sg.Txt('Beta :',font = ("Serif", 13))],
           [sg.Txt(size=(30,4), key='-OUTPUT2-')  ],
           [sg.Txt('C :',font = ("Serif", 13))],
           [sg.Txt(size=(30,4), key='-OUTPUT3-')  ],
           [sg.Button('Calculate', bind_return_key=True,font = ("Serif", 13))]]

#for now will only show the name of the chosen file
image_viewer_column = [
    
    [sg.Text("Output will be stored as .csv and .pdf in the same folder as exe file.",font = ("Serif", 13))],    
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

window = sg.Window('Poluborina-Kochina Solution', layout)

while True:
    event, values = window.read()

    if event != sg.WIN_CLOSED:
        try:
            H = float(values['-H-'])
            L = float(values['-L-'])
            H0, res, xz_array = PbK_solution(H,L)
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
        
        filename = f"H{H}_L{L}_H0_{calc}.png"
        window["-TOUT-"].update(filename,font = ("Serif", 13))
        window["-IMAGE-"].update(filename=filename)
    else:
        break