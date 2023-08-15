# AdsorptionSampling
The codes here expedite the running of AIMD liquid metal (LM) sampling of adsorption. In particular, they help with setting up the "shadow" molecular 
dynamics run that removes the adsorbate from LM surfaces and calculates single point energies of just the LM alone. This gives the reference state for 
calculation of adsorption energies. 

# setupShadowMD.py
This code takes a series of images (e.g. Trajectory, or concatenated OUTCAR), removes the adsorbate from them, then copies the VASP input files from a template 
calculation (in the same directory) to set up a series of single-point energy evaluations that can be run using a for-loop. The code is currently set to treat
anything with "C" "O" or "H" elements as an adsorbate, but this can be modified as desired. By default, the code will start from image zero in the trajectory
and run until the last image. If this is not desired, the start and end point and spacing can all be set.

# collateResultsShadowMD.py
This script will extract the energies from a shadow MD run, and pack them all into one file. Note that it uses the Henkelman group's VTST scripts.

# plotShadowMD.py
This script makes a plot of the adsorption energy over sampling-time. These are the plots visible in the Ruffman et al's main manuscript. The code has a translation
energy which should be manually set, in order to ensure a consistent reference in multi-step pathways (see manuscript).

# StepEnergies.py and EnergyHistogram.py
Small scripts for selecting images/energies at specific intervals from an AIMD trajectory and plotting these in a histogram. These scripts can be used for 
convenient data processing, and to calculate the translation energy.
