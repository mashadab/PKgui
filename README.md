# Polubarinova-Kochina-solutions
## ([Standalone GUI executables for Windows and Mac](https://drive.google.com/drive/u/0/folders/184aby8uWy1ZTMGidqhQwq9rfznjocVY1))
Oden Institute for Computational Engineering and Sciences / Jackson School of Geosciences / University of Texas Institute for Geophysics
The University of Texas at Austin


## Getting Started

### Overview

This codes solves the Polubarina-Kochina equations for low aspect ratio dam problems where upper lake level is close to the length of the dam. The GUI standalone executables are compatible with Mac and Windows (available [here](https://drive.google.com/drive/u/0/folders/184aby8uWy1ZTMGidqhQwq9rfznjocVY1)). The results have been validated with [1].

![cover](/cover/cover.png?raw=true)

### Variables [Dimension] (refer to paper [1])

L: Length of dam [L],   H: Lower lake height [L],   H1: Upper lake height [L],   H0: Seepage face height [L],   C: Pbk parameter [L], \
Alpha: Pbk parameter [-],   Beta:  Pbk parameter [-],   z: Free surface height [L],   x: horizontal dimension [L],   Q: Total flow rate by width in third dimension [L^2/T], Q_H0/Q_H: Ratio of flow rate from seepage face vs from lake [-]

### Relevant file paths
- [Standalone executable GUI](https://drive.google.com/drive/u/0/folders/184aby8uWy1ZTMGidqhQwq9rfznjocVY1)
- [Python with GUI](https://github.com/mashadab/Polubarinova-Kochina-solutions/blob/main/Pbk_GUI_latest.py)

## Quick Usage
### Executable GUI (Windows and Mac)

1. Download the executables from [Google drive](https://drive.google.com/drive/u/0/folders/184aby8uWy1ZTMGidqhQwq9rfznjocVY1).

2. Run the exe file (takes a few seconds to boot).:

    a. Mac: - If it does not recognize the file, then:
            - Change permissions of the file using ```chmod 777 Pbk_Mac```
            - Choose Apple menu > System Preferences, click Security & Privacy, then click General. You can grant an exception for a blocked app by clicking the “Open Anyway” button in the General pane (takes about 2-10 minutes to appear).
            
    b. Windows: Run as Administrator

### Python with GUI (Windows, Mac and Linux)
Dependences for Python program
- [Python](https://www.python.org/) version 3.5+
- [Numpy](http://www.numpy.org/) >= 1.16
- [scipy](https://www.scipy.org/) >=1.5
- [PySimpleGUI](https://pypi.org/project/PySimpleGUI/) >= 4.55.1
- [Matplotlib](https://matplotlib.org/) >= 3.5.0


1. Install the dependencies in a "Spyder" environment:

    i. For windows, remove scipy and numpy installed in conda and install using pip 
```
conda remove --force numpy, scipy
conda install python=3.7
pip install numpy
pip install scipy
pip install matplotlib
pip install PySimpleGUI
```

2. Download the github repository and unzip the package contents or clone the repository.
```
git clone https://github.com/mashadab/Polubarinova-Kochina-solutions
```

3. Run the python program "Pbk_GUI.py" in Spyder environment


## Output (Windows and Mac)
The output files include:

1. details.csv: Information regading the input and output variables
2. free-surface-profiles_XandZ.csv: Free surface profile (X vs Z)
3. free-surface-profile.pdf: High-res image of the output figure
4. free-surface-profile.png: Low-res image of the output figure

## Authors
- Mohammad Afzal Shadab ([mashadab@utexas.edu](mailto:mashadab@utexas.edu))
- Eric Hiatt ([eric.hiatt@utexas.edu](mailto:eric.hiatt@utexas.edu))
- Marc Andre Hesse ([mhesse@jsg.utexas.edu](mailto:mhesse@jsg.utexas.edu))

<!--- Cite the code: [![DOI](https://zenodo.org/badge/373661080.svg)](https://zenodo.org/badge/latestdoi/373661080) -->


## References / Related publications
[1] Hornung, U. and Krueger, T., 1985. Evaluation of the Polubarinova‐Kochina formula for the dam problem. Water Resources Research, 21(3), pp.395-398.
[2] Polubarinova-Koch, P.I., 2015. Theory of ground water movement. Princeton university press.
