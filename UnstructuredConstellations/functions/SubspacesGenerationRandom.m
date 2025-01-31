function V = SubspacesGenerationRandom(T,M,L)
%SUBSPACES_GENERATION generates L random subspaces of dimension T x M 
% and returns their bases.
%
% INPUTS
% Name: T
% Type: positive integer scalar
% Description: subspaces row dimension 
% 
% Name: M
% Type: positive integer scalar 
% Description: subspaces column dimension
% 
% Name: L
% Type:positive integer scalar
% Description: number of subspaces 
% 
% OUTPUT
% Name: V
% Type: TxMxL tensor
% Description: subspaces bases
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

V = zeros(T,M,L);
for i = 1:L
    X = 1/sqrt(2) * (randn(T,M) + 1i * randn(T,M)); % X(i,j) ~ CN(0,1)
    
    [Q,~] = qr(X); % compute QR decomposition
    
    v = Q(:,1:M); % select subspace basis
    
    V(:,:,i) = v; % store subspace basis
end

end