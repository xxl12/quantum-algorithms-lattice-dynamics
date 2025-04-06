%
% simulation of wave packet using quantum simulation algorithms 
% for 2d graphene lattice  
%
% Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
%
%
% Xiantao Li, Penn State U, March 2025
%

clear all

load fc 

load Qmatr2  % matrix coefficients of Q(z)  

m=size(Q{2,2},1); % dim of the matrix 
mdim=m;
max_degree=2; % degree of Q

% size of the lattice system 
Nx=128; 
Ny=36;

% find the matrix Q, Q=B^T

Nx1= Nx + max_degree;
Ny1= Ny + max_degree; 
B= sparse(Nx*Ny*mdim,Nx1*Ny1*mdim*r);

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


% initialize the displacement and velocity

u= zeros(Nx*Ny*mdim,1);
v= u;

% lattice structure of graphene 
a0=  1.46060011237d0;
b = sqrt(3d0)/2d0;
c = 1/2d0;

L= a0 * 2d0 * b;
H= a0 * 2d0 * (1d0 + c);

x0=[0 0;
    0 1;
    b 1+c;
    b 2+c]*a0;

% width of the wave packet
sigma=L*6;

% determine the phonon frequency and polarization vector

load fc.mat

Pmat = cell(3,3);
for i = 1:3
    for j = 1:3
        Pmat{i,j} = ( squeeze( F(i,j,:,:) ) + squeeze( F(4-i,4-j,:,:) )' )/2;
    end
end


xi = pi/4;
Dxi=0;
for i=-1:1
    for j=-1:1
        Dxi= Dxi + Pmat{i+2,j+2}*exp(-1i*xi*i); 
    end
end
Dxi= (Dxi+Dxi')/2;
[V, D]= eig(Dxi); omega=sqrt(abs(D(1,1))); w=(V(:,1)); 



% initialize the displacement and velocity and grid points
xgrid=zeros(Nx*Ny*mdim/2,1);;
ygrid=zeros(Nx*Ny*mdim/2,1);
for i=1:Nx
    for j=1:Ny

        nij= grid_2d_index(i, j, Nx, Ny);


        xgrid( (nij-1)*4+1:nij*4 )= (i-1)*L + x0(:,1);
        ygrid( (nij-1)*4+1:nij*4 )= (j-1)*H + x0(:,2);

        xc= L*Nx/4;
        yc= H*Ny/2;

        g= exp(- ((L*(i-1)-xc)^2+((j-1)*H-yc)^2)/(2*sigma^2)  );
        u( (nij-1)*mdim+1:nij*mdim )= g*cos( (i-1) *xi )*w;
        v( (nij-1)*mdim+1:nij*mdim )=g*sin( (i-1) *xi )*w*omega;

    end
end

% create a void 
n_remove=[];
g_remove=[];
for i=1:Nx
    for j=1:Ny
        if i>Nx-Nx/3 && i<5/6*Nx && abs(j-Ny/2)<Ny/4

            nij= grid_2d_index(i, j, Nx, Ny);
            
            n_remove=[ n_remove (nij-1)*mdim+1:nij*mdim ];

            g_remove=[ g_remove (nij-1)*4+1:nij*4];
        
            
        end
    end
end

B(n_remove,:)=[];

u(n_remove)= [];
v(n_remove)= [];

xgrid(g_remove)= [];
ygrid(g_remove)= [];

tri=delaunay(xgrid,ygrid);

% update the mesh
centroids = zeros(size(tri,1), 2);
for i = 1:size(tri,1)
    centroids(i,:) = mean([xgrid(tri(i,:)) ygrid(tri(i,:))], 1);
end

inHole = centroids(:,1) >= Nx*2/3*L & centroids(:,1) <= 5/6*Nx*L & ...
         centroids(:,2) >= 1/4*Ny*H & centroids(:,2) <= 3/4*Ny*H;
         
% Keep only triangles outside the hole
tri = tri(~inHole, :);


% build the Hamiltonian
[nr nc]= size(B);
H= sparse(nr+nc,nr+nc);
H(1:nr,nr+1:end) = -B;
H = H + H'; 
psi = [v; 1i*B'*u]; 

dt=1/32;
nsteps= 1024;

V= [];


for n=1:nsteps
    
    %    psi = U*psi; using RK4

    k1 = -1i * H *  psi;
    k2 = -1i * H * (psi + dt*k1/2);
    k3 = -1i * H * (psi + dt*k2/2);
    k4 = -1i * H * (psi + dt*k3);
        
    psi = psi + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    if mod(n,64)==0 
        v= psi(1:nr);  V=[V v];
    end

end

save wave2dqa.mat tri xgrid ygrid V 

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



