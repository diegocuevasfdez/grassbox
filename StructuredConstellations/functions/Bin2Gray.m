function g = Bin2Gray(bits, B)

%BIN2GRAY Gray encoding of a binary sequence.  
%
%   G = BIN2GRAY(BITS, B) applies Gray encoding to a sequence of bits
%   partitioned in groups of B bits.
%
%   INPUTS:
%   Name: bits
%   Type: 2*M*(T-M)*Bx1 binary vector
%   Description: sequence of bits to be Gray encoded.
%
%   Name: B
%   Type: positive integer scalar
%   Description: size of the partition of bits to be Gray encoded each time.
%
%   OUTPUT:
%   Name: g
%   Type: 2*M*(T-M)*Bx1 binary vector
%   Description: Gray encoded bits
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

Nbits = length(bits); % This should be 2*M*(T-M)*B
baux = reshape(bits,B,round(Nbits/B)); % bit sequences as rows
g = [baux(1,:); xor(baux(1:end-1,:),baux(2:end,:))];
g = g(:);

end
