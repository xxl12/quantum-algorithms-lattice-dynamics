%
% Matrix spectral factorization 
% for 1d nearest and next nearest neighbor interactions 
%
% Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
%
%
% Xiantao Li, Penn State U, March 2025
%


% The Toeplitz matrix D constructed from the force constants  
N=16;
D=5/3*diag(ones(N,1)) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1) ...
    + 1/6*diag(ones(N-2,1),2) + 1/6*diag(ones(N-2,1),-2);

% compute the matrix Q from the Freje-Rietz theorem 
B=zeros(N,N+2);
for n=1:N
    B(n,n:n+2)= [ (1+1/sqrt(3))/2 -1  (1-1/sqrt(3))/2  ];
end
Q = B';

fprintf('The norm of D-Q^T*Q:  %f\n', norm(D-Q'*Q) )




