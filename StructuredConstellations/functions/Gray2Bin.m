function bits = Gray2Bin(g,B)

%GRAY2BIN Gray to binary conversion of a bit sequence.  
%
%   BITS = GRAY2BIN(G, B) applies Gray decoding to a sequence of bits
%   partitioned in groups of B bits.
%
%   INPUTS:
%   Name: g
%   Type: 2*M*(T-M)*Bx1 binary vector
%   Description: sequence of bits to be Gray decoded.
%
%   Name: B
%   Type: positive integer scalar
%   Description: size of the partition of bits to be Gray decoded each time.
%
%   OUTPUT
%   Name: bits
%   Type: 2*M*(T-M)*Bx1 binary vector
%   Description: binary-coded bits
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

Nbits = length(g); % This should be 2*M*(T-M)*B
gaux = reshape(g,B,round(Nbits/B)); % bit sequences as rows
bits = zeros(size(gaux));
bits(1,:) = gaux(1,:);
for i = 2 : B
    bits(i,:) = xor(bits(i-1,:), gaux(i,:));
end
bits = bits(:);

end