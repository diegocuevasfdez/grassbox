function X = CubeSplitEncoding(B, cell_index, tx_symbols)

%CUBESPLITENCODING Encoding for the Cube-Split constellation. This
%constellation design was proposed in the following paper:
% 
% [1] K. -H. Ngo, A. Decurninge, M. Guillaud and S. Yang, "Cube-Split: A
% Structured Grassmannian Constellation for Non-Coherent SIMO
% Communications," in IEEE Transactions on Wireless Communications, vol.
% 19, no. 3, pp. 1948-1964, Mar. 2020, doi: 10.1109/TWC.2019.2959781.
%
%   X = CUBESPLITENCODING(B, CELL_INDEX, TX_SYMBOLS) generates the 
%   Cube-Split encoded codeword for a certain number of bits "B" per 
%   lattice and for the desired "cell_index" and "tx_symbols".
%
%   INPUTS:
%   Name: B
%   Type: positive integer scalar
%   Description: number of bits used to codify each grid. If the code uses
%   B bits per coordinate, each grid will have P = 2^B points. We encode a
%   total of 2(T-1) coordinates, where T is the coherence time.
%
%   Name: cell_index
%   Type: positive integer scalar
%   Description: index of the cell onto which we will apply the mapping.
%   It takes values between 0 and T-1.
%
%   Name: tx_symbols
%   Type: 2(T-1)x1 positive integer vector
%   Description: column vector with 2(T-1) coordinates. Each element 
%   represents the constellation point for one coordinate.
%
%   OUTPUT:
%   Name: X
%   Type: Tx1 double precision vector
%   Description: encoded point in Gr(1,C^T). 
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

% Grid in the Euclidean space [1, Eq. (9)]
Aj = 1/2^(B+1):1/2^B:(2^(B+1)-1)/2^(B+1); % grids
Aj_chosen = Aj(tx_symbols); % chosen symbols on the grid

% CUBE-SPLIT MAPPING

% Inverse transform sampling to generate N(0,1)
% [1, Eqs. 12-13]
real_a = Aj_chosen(1:2:end); real_a = real_a(:);
imag_a = Aj_chosen(2:2:end); imag_a = imag_a(:);

real_an = norminv(real_a);
imag_an = norminv(imag_a);
w = real_an + 1i*imag_an;

t = sqrt((1-exp(-abs(w).^2/2))./(1+exp(-abs(w).^2/2))) .* w ./ abs(w);

% Apply mapping g_i
% [1, Eq. (11)]
array1 = t(1:cell_index-1);
array2 = t(cell_index:end);
t_1 = [array1; 1; array2];

gi_a = 1 / sqrt(1+sum(abs(t).^2)) * t_1;

X = gi_a; % TX codeword

end