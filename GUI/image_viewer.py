#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pbk GUI
"""

import PySimpleGUI as sg
import os

#window layout of two columns
file_list_column = [
    [
    sg.Text("Image Folder"),
    sg.In(size=(25,1), enable_events=True, key="-FOLDER-"),
    sg.FolderBrowse(),
    ],
    [
        sg.Listbox(
             values=[], enable_events=True,size=(40,20),
             key="-FILE LIST-"
         )
     ],
]

#for now will only show the name of the chosen file
image_viewer_column = [
    
    [sg.Text("Choose an image from the list on the lest:")],    
    [sg.Text(size=(40,1), key="-TOUT-")],
    [sg.Image(key="-IMAGE-")],
]

# ----- full layout-----------
layout = [
    [
         sg.Column(file_list_column),
         sg.VSeparator(),
         sg.Column(image_viewer_column),
     ]    
]

#making a window
window = sg.Window("Image Viewer", layout)

#event loop
while True:
    event, values = window.read()
    if event == "Exit" or event == sg.WIN_CLOSED:
        break
    
    #folder name was filled in, so make a list of files
    if event == "-FOLDER-":
        folder = values["-FOLDER-"]
        try:
            #get list of files in folder
            file_list = os.listdir(folder)
        except:
            file_list = []
        fnames= [
            f
            for f in file_list
            if os.path.isfile(os.path.join(folder, f))
            and f.lower().endswith((".png",".gif"))
        ]
        window["-FILE LIST-"].update(fnames)
    elif event == "-FILE LIST": #a file was chosen from the list
        try:
            filename = os.path.join(
                values["-FOLDER"], values["-FILE LIST-"][0]    
            )
            window["-TOUT-"].update(filename)
            window["-IMAGE-"].update(filename=filename)
        except:
            pass
#window.close()



'''
# image_viewer.py

import io
import os
import PySimpleGUI as sg
from PIL import Image


file_types = [("JPEG (*.jpg)", "*.jpg"),
              ("All files (*.*)", "*.*")]

def main():
    layout = [
        [sg.Image(key="-IMAGE-")],
        [
            sg.Text("Image File"),
            sg.Input(size=(25, 1), key="-FILE-"),
            sg.FileBrowse(file_types=file_types),
            sg.Button("Load Image"),
        ],
    ]

    window = sg.Window("Image Viewer", layout)

    while True:
        event, values = window.read()
        if event == "Exit" or event == sg.WIN_CLOSED:
            break
        if event == "Load Image":
            filename = values["-FILE-"]
            if os.path.exists(filename):
                image = Image.open(values["-FILE-"])
                image.thumbnail((400, 400))
                bio = io.BytesIO()
                image.save(bio, format="PNG")
                window["-IMAGE-"].update(data=bio.getvalue())

    window.close()


if __name__ == "__main__":
    main()
'''