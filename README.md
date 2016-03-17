# UQ-mesa
Uncertainty Quantification (UQ) in Modules and Experiments for Stellar Astrophysics (MESA) 
The scripts, doc, and otherwise found in this repo are related to my UQ study of the stellar evolution code, MESA.\n
Some of my analyzing scripts utilise the mesa module for python from NuGridPy (nugridpy.phys.uvic.ca) and should be noted in the comments of the code. If you're interested in this project, you'll also find D. Willcox's work handy, and can be found at https://github.com/dwillcox/nuq\n
- isthecodedone.py is a script that runs through subdirectories of a a MESA run folder and prints out which folders are done and what their termination code was.\n
- scatter_w_inputs.py uses the mesa module to extract final parameter values and extracts inlist information. User inputs desired x axis, y axis, and colormap. \n
- read_n_scatter.py is the basis for scatter_w_inputs except it only plots Blocker v. Mass with Reimers colormap.