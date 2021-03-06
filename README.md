# ASTEC

This file contains instructions for installing and running the ASTEC reconstruction algorithm.
ASTEC is fully described in the following article: [Contact-dependent cell communications drive morphological invariance during ascidian embryogenesis](https://www.biorxiv.org/content/early/2017/12/24/238741)

## I - CONTENTS OF THE  REPOSITORY
Once uncompressed, the folder contains the following elements:
  - ASTEC :  main source code for the ASTEC package
  - LICENCE   : the licence terms you accept by using the workflow
  - index.ipynb : a notebook example which shows how to run the different functions of the ASTEC segmentation and tracking algorithm
  - README.TXT    : this file

## II - INSTALLATION AND SOFTWARE REQUIREMENTS

ASTEC was fully tested on Unbuntu 14.04 64 bits.
### II.i - REQUIREMENTS
In order to be able to compile the different source code you need to install several packages, mostly from the terminal application.

to run the different codes:
  * python 2.7 or higher  
    - Installation : Should be installed, use the command line `python -V` to get the installed version or visit https://www.python.org/
  * and a c compiler
    - Installation : Should be installed, use the command line `dpkg --list | grep compiler` to get the list of avaiblable C compiler versions or install gcc with `sudo apt-get install gcc`

the following libraries are necessary:
  * pip, an installer for python (https://pypi.python.org/pypi/pip)
    - Installation : run command line `sudo apt-get install python-pip python-dev build-essential` 
  * cmake , a software to build source codes (http://www.cmake.org)
    - Installation : run command line `sudo apt-get install cmake` 

the following libraries are necessary:
  * numpy,scipy,matplotlib, different scientific packages for python  (http://www.numpy.org, http://www.scipy.org, http://matplotlib.org)
    - Installation : run command line `sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose`
  * jupyter notebook, a open-source framework to share scientific code (http://jupyter.org)
    - Installation : run command line `sudo pip install jupyter`
  * zlib, a compression library	(http://www.zlib.net) :
    - Installation : run command line `sudo apt-get install zlib1g-dev`
  *libhdf5-dev,cython,h5py  a library to read/write hdf5 images (https://www.hdfgroup.org/HDF5/)
    - Installation : run command line `sudo apt-get install libhdf5-dev`, 
                then run `sudo pip install cython`, 
                then run `sudo pip install h5py`
  *pylibtiff, a library to read/write tiff/tif images (https://pypi.python.org/pypi/libtiff/)
    - Installation : run `cd path/to/Package/ASTEC/CommunFunctions/libtiff-0.4.0/; sudo python setup.py install`

### II.ii - INSTALLATION
In order to compile the necessary c folders:
in a terminal:
 - run `cd path/to/Package/ASTEC/CommunFunctions/`
 - run `./compil.sh`


### II.iii - TROUBELSHOOTING
The installation of jupyter may failed according of the version of your ipython : to correct this issue run `sudo pip install IPython==5.0 `



## III - RUNNING ASTEC
The code of ASTEC, MARS and the FUSION algorithm are written in python and can easly be run with a web interface using a Notebook (http://jupyter.org)  

You can find a user interface tutorial here (http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb)

To start the notebook, open a terminal and  run the command `jupyter notebook path/to/Package/index.ipynb`

