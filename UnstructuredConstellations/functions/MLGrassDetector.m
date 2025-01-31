function codewordIndex = MLGrassDetector(GrassCodebook,Y)

% MLGRASSDETECTOR Description: This function performs maximum-likelihood 
% detection for non-coherent communications based on unstructured 
% Grassmannian constellations.
%
%   MLGRASSDETECTOR(GRASSCODEBOOK,Y) returns the index of the estimated
%   received codeword that belongs to codebook "GrassCodebook" from the
%   received signal "Y".
%
%   INPUTS:
%   Name: GrassCodebook
%   Type: TxMxNumCodewords tensor
%   Description: The elements are orthonormal bases of a Grassmann
%   constellation
%
%   Name: Y
%   Type: TxN double precision matrix (either real or complex)
%   Description: received signal
%
%   OUTPUT:
%   Name: codewordIndex
%   Type: positive integer scalar
%   Description: index of the ML detected codeword 
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

[~,~,NumCodewords] = size(GrassCodebook);
loglikelihood = zeros(1,NumCodewords);

for nn = 1:NumCodewords 
    loglikelihood(nn) = real(trace(Y'*GrassCodebook(:,:,nn)*GrassCodebook(:,:,nn)'*Y));
end

[~,codewordIndex] = max(loglikelihood);

end
