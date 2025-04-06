%
% Matrix spectral factorization 
% for 2d FCC lattice of  graphene lattice Aluminum
%
% Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
%
%
% Xiantao Li, Penn State U, March 2025
%
	   
clc;

%% PARAMETERS
m = 12;       % size of matrices P_{jk} (here, 8x8)
r = 2;       % number of squares (adjust as needed)
deg = 2;     % maximum degree in z and w for each Q_l(z,w)


% the input matrix coefficients P_{jk} for j,k=-1,0,1.
% computed from embedded atom potential for aluminum 

load PmatAl

% find Q
Q=optimize_trivariate_spectral_factorization_deg3(Pmat, r);

save Qmatr2 Q r
