# corsing-wavelet-fourier-adr

Matlab code used to generate the figures in "Wavelet-Fourier CORSING techniques for multi-dimensional advection-diffusion-reaction equations" by Simone Brugiapaglia, Stefano Micheletti, Fabio Nobile, Simona Perotto, 2018.

## What do I need to run this code?     

In order to run this code, you need the following toolboxes:

1. Matlab Symbolic toolbox

2. OMP-Box (v10 or later)
   You can download this package from 
   http://www.cs.technion.ac.il/~ronrubin/software.html

## How to run this code?                                

Running this code is extremly simple! Any Figure X in the paper can be reproduced by running the script FigureX.m. In order to make this efficient, some of the figures are produced from data previously computed and stored in the folder data in FigureX_article.mat. If you want to run those hevier experiments again, you can run the script Compute_FigureX.m and then load the new data file FigureX_new.mat in the corresponding script FigureX.m. 
To print all figures, run the script Print_all_figures.m.
