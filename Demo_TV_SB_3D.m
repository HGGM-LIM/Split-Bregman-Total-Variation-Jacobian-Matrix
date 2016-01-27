% Demo_TV_SB_3D.m 
%
% Demo to solve the 3D constrained total variation image reconstruction
% problem using the Split Bregman formulation. 
%
% The Split Bregman separates L2- and L1-norm functionals in such a way
% that they can be solved analytically in two alternating steps. In the
% first step a quadratic problem is solved using a Gauss-Newton Krylov
% method, and the second step is solved using a shrinkage formula. 
%
% TV_SB_3D.m solves the constrained total variation problem 
% min_u ||grad_x,y u||_1 st. ||Au-f||^2 < delta, where A is a linear
% operator (a matrix) the projects the image u to the data f. The code
% works for general linear inverse problems. It currently expects A to be a
% matrix; it can be easily modified to use A as a MATLAB function by
% changing A and A' for functions that compute forward and adjoint
% operators.
%
% This demo solves the compressed sensing problem for magnetic resonance
% imaging as an exemplar. A is the Fourier transform provided as a matrix
% operator, f is undersampled data, and u is a 3D image.
%
% If you use this code, please, cite the following paper: 
% Abascal JF, Chamorro-Servent J, Aguirre J, Arridge S, Correia T, Ripoll
% J, Vaquero JJ, Desco M. Fluorescence diffuse optical tomography using the
% split Bregman method. Med Phys. 38(11):6275-84, 2011.   
% DOI: http://dx.doi.org/10.1118/1.3656063
%
% Code downloaded from repository: 
% https://github.com/HGGM-LIM/Split-Bregman-Total-Variation-Jacobian-Matrix
% 
% Juan FPJ Abascal, Judit Chamorro-Servent
% Departamento de Bioingeniería e Ingeniería Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, juchamser@gmail.com, desco@hggm.es

% This method is a modification of Goldstein'n code mrics.m downloaded from 
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323–343. 

N       = [20 20 10]; 
sparsity = .85; % use only 30% on the K-Space data for CS 

% build an image of a square
image   = zeros(N);
image(N(1)/4:3*N(1)/4,N(2)/4:3*N(2)/4,N(3)/4:3*N(3)/4)=255;
figure; imagesc(abs(image(:,:,5))); 

% Fourier transform as matrix operator, slice by slice FFT2 
% Change this by a linear operator or Jacobian 
A       = zeros(prod(N),prod(N));
for ip = 1:N(3)
    indThis     = 1+prod(N(1:2))*(ip-1):prod(N(1:2))*(ip);
    A(indThis,indThis)  = kron(fft(eye(N(1:2))),fft(eye(N(1:2)))); 
end

% Simulate data. For an operator substitute for A*image
data    = A*image(:);

% build the sampling matrix, R
% Without undersampling, make R=ones(size(data)) or do not pass R to TV_SB_3D.m
rand('state',0);
R       = rand(N);
R       = R<sparsity;
R(1,1)  = 1;
R       = R(:); 

% Simulate the CS undersampled data
data    = R.*data;

% Recover the image
mu      = 1;
lambda  = 1;
gamma   = 1;
nInner  = 1;
nBreg   = 100;
target  = image;
[u,errAll] = TV_SB_3D(A,data,N,mu,lambda,gamma,nInner,nBreg,target,R);

% Display results
figure;
subplot(2,2,1);
imagesc(abs(image(:,:,5))); colormap('gray'); colorbar;
title('Original');
subplot(2,2,2);
imagesc(abs(u(:,:,5))); colormap('gray'); colorbar;
title('Recovered');
subplot(2,2,3);
plot(errAll); axis tight; title(['Sol. error' ]);

%