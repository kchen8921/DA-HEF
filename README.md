[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kchen8921/DA-HEF/master)

# Data Assimilation 

[![N|Solid](https://upload.wikimedia.org/wikipedia/en/thumb/1/17/Pacific_Northwest_National_Laboratory_logo.svg/200px-Pacific_Northwest_National_Laboratory_logo.svg.png)](https://www.pnnl.gov/)

In the project, data assimilation methods such as EnKF and ES-MDA are used to estimate hydrologic exchange flux between surface water and groundwater by assimilating observed field data at Hanford site (e.g. temperature, hydaulic heads...). This repository provides the entire workflow for the implementation of DA. Objectives include but not limited to:

  - Estimate permeability and thermal conductivity in the riverbed
  - Infer subdaily hydrologic exchange flux

# Contents
# 
| Workflow | Contents |
| ------ | ------ |
| 1. [Pre-Processing](https://github.com/kchen8921/SFA-DA/blob/master/pre-processing.ipynb) | 1.1 Preprocess raw voltage data |
|                   | 1.2 Convert voltage to temperature |

| 2. [Data Assimilation](https://github.com/kchen8921/SFA-DA/blob/master/Data%20Assimilation.ipynb) | 2.1 Data Assimilation Framework |
|                      | 2.2 Configuration of forward simulator |
|                      | 2.3 Configuration of observation data |
|                      | 2.4 Input setting and execution |
|                      | 2.5 Submit batch job on cluster |



# Installation and Configuration

1. The workflow is written in [Jupyter notebook](http://jupyter.org/) which supports both Python and R. A recommended distribution of Jupyter notebook is [Anaconda](https://www.anaconda.com/download/).
  (1) To start Jupyter notebook on Mac/Linux after installation of Anaconda, typing the following command in terminal:
    ```sh
    jupyter notebook
    ```
    (2) To start Jupyter notebook on Windows, just click the desktop icon. 
2. The numerical model for subsurface flow at Hanford site is built using [PFLOTRAN](http://www.pflotran.org/), a massively parallel reactive flow and transport model for describing surface and subsurface processes.
3. The workflow is adapted to the supercomputers at [National Energy Research Scientific Computing Center (NERSC)](http://www.nersc.gov/).

![Estimated flux](https://github.com/kchen8921/SFA-DA/blob/2799ad89472f6a627e0631474edd6be1562277b1/model/TH1D/figure/ENKF-ESMDA.png)
![Estimated perm](https://github.com/kchen8921/SFA-DA/blob/master/model/TH1D/figure/flux%20estimation.png)

