#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk GUI
"""

import PySimpleGUI as sg
from Pbk_main_function import *

sg.theme('Dark Green 7')


#window layout of two columns
file_list_column = [ [sg.Txt('Length of the aquifer, L')],
           [sg.In(size=(8,4), key='-L-')],
           [sg.Txt('Lake level, H')],
           [sg.In(size=(8,4), key='-H-')],
           [sg.Txt('Seepage face height, H0 :')],
           [sg.Txt(size=(20,4), key='-OUTPUT-')  ],
           [sg.Txt('Alpha :')],
           [sg.Txt(size=(20,4), key='-OUTPUT1-')  ],
           [sg.Txt('Beta :')],
           [sg.Txt(size=(20,4), key='-OUTPUT2-')  ],
           [sg.Txt('C :')],
           [sg.Txt(size=(20,4), key='-OUTPUT3-')  ],
           [sg.Button('Calculate', bind_return_key=True)]]

#for now will only show the name of the chosen file
image_viewer_column = [
    
    [sg.Text("Output figure will appear here.")],    
    [sg.Text(size=(40,1), key="-TOUT-")],
    [sg.Image(key="-IMAGE-")],
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

        window['-OUTPUT-'].update(calc)
        window['-OUTPUT1-'].update(calc1)
        window['-OUTPUT2-'].update(calc2)
        window['-OUTPUT3-'].update(calc3)
    else:
        break