function V = subspaces_generation_random(T,M,N)
%SUBSPACES_GENERATION
%This function generates N random subspaces of dimension T x M and returns 
%their bases.
%
% INPUTS
% T: subspaces row dimension
% M: subspaces column dimension
% N: number of subspaces 
% 
% OUTPUT
% V: subspaces bases 

V = zeros(T,M,N);
for i = 1:N
    X = 1/sqrt(2) * (randn(T,M) + 1i * randn(T,M)); % X(i,j) ~ CN(0,1)
    
    [Q,~] = qr(X); % compute QR decomposition
    
    v = Q(:,1:M); % select subspace basis
    
    V(:,:,i) = v; % store subspace basis
   
end

end