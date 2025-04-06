%
% Matrix spectral factorization 
% for 2d graphene lattice
%
% Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
%
%
% Xiantao Li, Penn State U, March 2025
%

clear all;

%% PARAMETERS
m = 8;       % size of matrices P_{jk} (here, 8x8)
mdim=m; 
r = 2;       % number of squares in Fejer-Rietz (adjust as needed)
deg = 2;     % maximum degree in z1 and z2 for each Q_l(z1,z2)
tol = 1e-1;  % tolerance for approximate coefficient matching
max_degree=2;% degree of Q 

%% Define your input matrix coefficients P_{jk} for j,k=-1,0,1.
% For demonstration purposes, we create a cell array Pmat (3x3)
% with each entry an m-by-m matrix.
% Replace these with your actual data.

load fc  % load the force constant matrix

% arrange the force constant matrix P_{i,j}, i,j =-1,0,1 into 3 by 3 cells		  
Pmat = cell(3,3);
for i = 1:3
    for j = 1:3
        Pmat{i,j} = ( squeeze( F(i,j,:,:) ) + squeeze( F(4-i,4-j,:,:) )' )/2;
    end
end

% find the factors Q, stored as a cell Q{l,j,k}, l=1,2,..,r 
Q = optimize_bivariate_spectral_factorization_deg2(Pmat, r);

save Qmatr2 Q r

% now check whether the factorization is correct

%first create the matrix D for a 10x10 system

Nx=10; 
Ny=10;

D= zeros(Nx*Ny*mdim,Nx*Ny*mdim);
for i=1:Nx
    for j=1:Ny

        nij= grid_2d_index(i, j, Nx, Ny);

        for k=-1:1
            for l=-1:1
                ik= i + k;
                jl= j + l;

                if ik >= 1 && ik <= Nx && jl >= 1 && jl <= Ny

                    nkl= grid_2d_index(ik, jl, Nx, Ny);

                    D( (nij-1)*mdim+1:nij*mdim, (nkl-1)*mdim+1:nkl*mdim )= Pmat{k+2,l+2}'; %squeeze( F(k+2,l+2,:,:) );

                end
            end
        end
    end
end

% now assemble the matrix Q= B^T
Nx1= Nx + max_degree;
Ny1= Ny + max_degree; 
B= zeros(Nx*Ny*mdim,Nx1*Ny1*mdim*r);

for i=1:Nx
    for j=1:Ny

        nij= grid_2d_index(i, j, Nx, Ny);

        for k=0:max_degree
            for l=0:max_degree
                ik= i + k;
                jl= j + l;

                nkl= grid_2d_index(ik, jl, Nx1, Ny1);
                    
                for a=1:r
                    nstart= (a-1)*Nx1*Ny1*mdim;
                    B( (nij-1)*mdim+1:nij*mdim, nstart+(nkl-1)*mdim+1:nstart+nkl*mdim )= Q{a,k+1,l+1}' ;
                end
            end
        end
    end
end
Q= B';
		
fprintf('The relative error of D-Q^T*Q:  %f\n', norm(D-Q'*Q, inf)/norm(D, inf) )


function n = grid_2d_index(i, j, Nx, Ny)
    % GRID_COLUMN_INDEX Maps 2D grid indices (i,j) to a linear index n when 
    % ordering points by columns
    %
    % Inputs:
    %   i - row index (1 to Nx)
    %   j - column index (1 to Ny)
    %   Nx - number of rows in the grid
    %   Ny - number of columns in the grid
    %
    % Output:
    %   n - linear index when ordering points column-wise
    
    % Validate input ranges
    if i < 1 || i > Nx || j < 1 || j > Ny
        error('Indices out of grid bounds');
    end
    
    % Calculate linear index by columns
    % The formula maps (i,j) to a sequential index moving down columns first
    n = (j - 1) * Nx + i;
end
		
		
