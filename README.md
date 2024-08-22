# KirasFM-NN: KirasFM + Feedforawrd Neural Network

This repository consists of a modified version of the Maxwell Solver KirasFM, which is based on a fixed-point domain decomposition method and a feedforward neural network-enhanced that is used to compute an approximation of the surface operator occurring in the fixed point formulation of the Maxwell equation.

## Citation
Please use the *Cite this repository* button in the *About* section of this repository.

## Installation - Dependencies
The dependencies required are deal.II, PyTorch, and Jupyter Lab. In the following, it is assumed that a recent C++ compiler (e.g., GCC) and a recent Python version are installed.

### Obtaining deal.II 
This is a CMake script that installs dea.lII along with its dependencies. And it should be the easiest way to install deal.II with all its dependencies. For more details, see https://github.com/kinnewig/dcs2.

tl;dr:
1. Step: Download dcs2:
```
git clone git@github.com:kinnewig/dcs2.git
cd dcs2
```
2. Step: Run the install script:
```
./dcs.sh  -b </path/to/build> -p </path/to/install>
```

Remember to replace `</path/to/build>` with the path where you would like to store the temporary files created while installing deal.II (the folder can be deleted once you successfully installed deal.II).

Also, remember to replace `</path/to/install>` with the folder where you would like to install deal.II.

If you have any problems feel free to open an issue on: https://github.com/kinnewig/dcs2/issues

### (Optional) Create a venv:
Before you begin with installing the required Python packages, it is recommended to create a venv (virtual environment)
```
python -m venv </path/to/new/virtual/environment>
```
Remember to replace `</path/to/new/virtual/environment>` with the path where you would like to store the venv.

After creating the venv, you can enter the venv by
```
. .</path/to/new/virtual/environment>/bin/activate
```

For more details on venv, see: https://docs.python.org/3/library/venv.html

### Install PyTourch
To install PyTorch with the CPU backend, use
```
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

For more details on PyTorch see https://pytorch.org/get-started/locally/ 

### Install Jupyter Lab 
To install Jupyter Lab, use
```
pip3 install jupyterlab
```

Jupyter Lab can be started by
```
jupyter lab
```
After that, Jupyter Lab should open in the browser.
