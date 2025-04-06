%
% matrix spectral factorization 
% 1d diatoms interactions 
%
% Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
%
%
% Xiantao Li, Penn State U, March 2025
%

clear all 

% number of unit cells 
N=16;  

% 2x2 force constant matrices; K_{-1}=K_1^T  
K1= [0  0; -1 0];
K0= [2 -1;-1 2];

% the matrix D
D= zeros(2*N,2*N);
for n=1:2:2*N-1
    D(n:n+1,n:n+1)= K0;
    if n<2*N-1
        D(n:n+1,n+2:n+3)= K1; 
    end
    if n>1
        D(n:n+1,n-2:n-1)= K1';
    end
end

% coefficient matrices from Q(z) from the Freje-Rietz theorem 
Q0= [1  0; 1 -1]; 
Q1= [0 -1; 0  0];

B=zeros(2*N,2*N+2);
for n=1:2:2*N-1
    B(n:n+1,n:n+1)  = Q0'; 
    B(n:n+1,n+2:n+3)= Q1'; 
end
Q= B'; 

% check results 
fprintf('The norm of D-Q^T*Q:  %f\n', norm(D-Q'*Q) )
