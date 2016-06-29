# UQ-mesa
Uncertainty Quantification (UQ) in Modules and Experiments for Stellar Astrophysics (MESA) 
The scripts, doc, and otherwise found in this repo are related to my UQ study of the stellar evolution code, MESA.
Some of my analyzing scripts utilise the mesa module for python from NuGridPy (nugridpy.phys.uvic.ca) and should be noted in the comments of the code. If you're interested in this project, you'll also find D. Willcox's work handy, and can be found at https://github.com/dwillcox/nuq
- isthecodedone.py is a script that runs through subdirectories of a a MESA run folder and prints out which folders are done and what their termination code was.
- scatter_w_inputs.py uses the mesa module to extract final parameter values and extracts inlist information. User inputs desired x axis, y axis, and colormap. 
- read_n_scatter.py is the basis for scatter_w_inputs except it only plots Blocker v. Mass with Reimers colormap.
- plot_w_inputs.py creates line plots for Blocker v. Mass (and organizes by Reimers in the process).
- make_hr.py makes either single HR plot for multiple folders, or combines multiple folder's HR info to one plot.
- quickplot.py reads .out files created by initial scatter_w_inputs.py to make plots quickly
- RandFolders.py makes a directory of folders where Reimers and Blocker values are random & unique to each folder