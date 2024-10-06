# Models, animations and Signatures of metastability

### 
*Code to accompany Metastability demystified June 2024*

<table><tr><td>Fran Hancock
fran.hancock@kcl.ac.uk
January 2024</td></tr></table>

# Animations

### Setup your local path

addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PostDoc/Projects/META_MODELLING'))

### If you just wish to run the animations without re-running the simulation
### You will still need to download the Brain Dynamics Toolbox here https://bdtoolbox.org/

In the modelling/models folder, load a model into MATLAB, e.g.

load HMM.mat sys
gui = bdGUI(sys)

### Kuramoto Shanahan model

#### to run the Kuramoto Shanahan model in the GUI:

Kij = random_connections; sys = KuramotoShan(Kij); gui = bdGUI(sys);

#### save the model from the GUI in models/

#### In the animations folder, run the simulation plot and save the video

Plot_phases_Shanahan_extra(4000,10000,10,1)


###  Hansel-Mato-Meurnier model

#### run the Hansel-Mato-Meunier model in the GUI

n=5; kij=rand(n,n); Kij=kij-diag(diag(kij)) + diag(1); sys=HMM(Kij); gui=bdGUI(sys);

#### In the animations folder, run the simulation plot and save the video

Plot_phases_HMM(1,2000,2,1)


### Extended Haken-Kelso-Bunz model

#### run the HKB model in the GUI

n=8;Kij=ones(n);Kij=Kij-diag(diag(Kij));
Kij=0.105*Kij;sys=HKB(Kij);gui=bdGUI(sys);

#### In the animations folder, run the simulation plot and save the video

Plot_phases_HKB(7500,10000,10,1)

# Signatures of metastability  

## Instructions

#### The code for the signatures is in the Signatures folder

AAL116 parcellated data for 20 Healthy controls and 20 patients with Schizophrenia are supplied
Source - COBRE - See Metastability as a candidate neuromechanistic marker for schizophrenia, 10.1371/journal.pone.0282707


Modify the **loadParameters.m** if you are running this code on your own data.

Otherwise

Run **get_signatures.m**

This file will calculate the 4 signatures in the manuscript and output the results in a csv file for statistical anaysis in R Studio

In the folder _statistics_, run **01_Global_metrics.Rmd** in RStudio

The statistics will be saved in excel under the folder _Results_

The figure will be saved in the folder _Figures_

For readability, the signatures have the following names in the figure:

std-KOP = Metastability

mean-VAR = Dynamical complexity

std-SPECT = Dynamical agility

std-IGNITE = Ignition flexibility 


#### NOTE: For MAC users, depending on which OS version you are running, the export to Excel may not work.

If all runs well, you should have the following figure in the Figures folder 
![Boxplots Healthy Controls versus Schizophrenia](kruskal_Global_Metrics.pdf)












