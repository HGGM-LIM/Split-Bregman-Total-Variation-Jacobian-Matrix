# Split-Bregman-Total-Variation-Jacobian-Matrix
This repository contains demos that show how to solve  2D and 3D constrained Total Variaton image reconstruction problems using the Split Bregman formulation. 

TV_SB_2D.m and TV_SB_3D.m solve the constrained total variation problem

![TV](https://github.com/HGGM-LIM/Split-Bregman-Total-Variation-Jacobian-Matrix/blob/master/TVFormula.jpg)

where A is a linear operator (a matrix) that projects the image u to the data f. The code works for general linear inverse problems. It currently expects A to be a matrix; it can be easily modified to use A as a MATLAB function by changing A and A' for functions that compute forward and adjoint operations.

These demos solve the compressed sensing problem for magnetic resonance imaging as an exemplar. A is the Fourier transform provided as a matrix operator, f is undersampled data, and u is a 2D image. 

The repository contains the following files:

- **Demo_TV_SB_2D.m:** Demo for 2D TV reconstruction

- **Demo_TV_SB_3D.m:** Demo for 3D TV reconstruction

- **TV_SB_2D.m:** 2D TV reconstruction

- **TV_SB_3D.m:** 3D TV reconstruction

If you use this code, please, cite the following paper: Abascal JF, Chamorro-Servent J, Aguirre J, Arridge S, Correia T, Ripoll
J, Vaquero JJ, Desco M. Fluorescence diffuse optical tomography using the split Bregman method. Med Phys. 38(11):6275-84, 2011. DOI: http://dx.doi.org/10.1118/1.3656063

If you need to contact the author, please do so at 
juanabascal78@gmail.com, juchamser@gmail.com, desco@hggm.es
