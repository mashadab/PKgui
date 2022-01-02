# Polubarinova-Kochina-solutions
Oden Institute for Computational Engineering and Sciences / Jackson School of Geosciences / University of Texas Institute for Geophysics
The University of Texas at Austin


## Getting Started

### Overview

This codes solves the Polubarina-Kochina equations for low aspect ratio dam problems where upper lake level is close to the length of the dam. The GUI standalone executables are compatible with Mac and Windows (available [here](https://drive.google.com/drive/u/0/folders/184aby8uWy1ZTMGidqhQwq9rfznjocVY1)). The results have been validated with [1].

![cover](/cover/cover.png?raw=true)

### GUI Standalone executables (available [here](https://drive.google.com/drive/u/0/folders/184aby8uWy1ZTMGidqhQwq9rfznjocVY1))
These files may require permissions to write:
- Mac: Change permissions using chmod 777 Pbk_Windows
- Windows: Run as Administrator

### Dependences

Polubarinova-Kochina-solutions python code requires the following packages to function:
- [Python](https://www.python.org/) version 3.5+
- [Numpy](http://www.numpy.org/) >= 1.16
- [scipy](https://www.scipy.org/) >=1.5
- [PySimpleGUI](https://pypi.org/project/PySimpleGUI/) >= 4.55.1
- [Matplotlib](https://matplotlib.org/) >= 3.5.0


### Variables (refer to paper [1])
L: Length of dam [L] \
H: Lower lake height [L] \
H1: Upper lake height [L]\
H0: Seepage face height [L]\
C: Pbk parameter \
Alpha: Pbk parameter \
Beta:  Pbk parameter \
z: Free surface height [L] \
x: horizontal dimension [L] \

### Running seepagePINN


### Quick Usage (MacOS)


1. Install the dependencies in a "Conda environment":

    i. Create an environment: conda create **environment name**\
    ii. Activate the environment: conda activate **environment name**\
    iii. Install the dependent libraries (given in dependencies): conda install **library name**
```
conda create -n seepage python=3.7
conda activate seepage 
conda install tensorflow==1.14
conda install matplotlib pandas scipy h5py
```
<!--
```
conda create -n seepage -c uvilla -c conda-forge fenics==2019.1.0 matplotlib scipy jupyter python=3.7
conda activate seepage
conda install -c conda-forge tensorflow==1.13.2
conda install -c conda-forge numpy=1.16.6 -y
conda install -c conda-forge pandas -y
conda install -c anaconda scipy=1.5.3 -y
conda install -c anaconda h5py=3.3.0 -y
```
-->
2. Download the github repository and unzip the package contents or clone the repository.
```
git clone https://github.com/dc-luo/seepagePINN.git
```
3. Move to the specific folder on steady results
```
cd seepagePINN/src/steady/paper/
```
4. Run the python program in Mac terminal using experimental_all.py [-h] [-c CASE] [-n N_EPOCH] [-m {dinucci,dupuit}] [-r]
for example:
```
python experimental_all.py -c 1mm -n 20000
```
and to visualize the training results
```
python viz_exp.py -c 1mm -u --show --legend
```
## Authors
- Mohammad Afzal Shadab
- Eric Hiatt
- Marc Andre Hesse

<!--- Cite the code: [![DOI](https://zenodo.org/badge/373661080.svg)](https://zenodo.org/badge/latestdoi/373661080) -->


## References / Related publications
[1] Hornung, U. and Krueger, T., 1985. Evaluation of the Polubarinova‐Kochina formula for the dam problem. Water Resources Research, 21(3), pp.395-398.
