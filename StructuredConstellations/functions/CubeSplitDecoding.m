function [rx_cell_index, rx_symbols] = CubeSplitDecoding(B, Y)

%CUBESPLITDECODING Decoding for the Cube-Split constellation. This
%constellation design was proposed in the following paper:
% 
% [1] K. -H. Ngo, A. Decurninge, M. Guillaud and S. Yang, "Cube-Split: A
% Structured Grassmannian Constellation for Non-Coherent SIMO
% Communications," in IEEE Transactions on Wireless Communications, vol.
% 19, no. 3, pp. 1948-1964, Mar. 2020, doi: 10.1109/TWC.2019.2959781.
%
%   [RX_CELL_INDEX, RX_SYMBOLS] = CUBESPLITDECODING(B, Y) estimates the 
%   transmitted cell index and the transmitted symbols from the received 
%   signal "Y" and for the chosen value of bits per grid "B".
%
%   INPUTS:
%   Name: B
%   Type: positive integer scalar
%   Description: number of bits used to codify each grid. If the code uses
%   B bits per coordinate, each grid will have P = 2^B points. We decode a
%   total of 2(T-1) coordinates, where T is the coherence time.
%
%   Name: Y
%   Type: TxN double precision vector
%   Description: TxN received signal, where T is the coherence time and N
%   the number of receive antennas. 
%
%   OUTPUT:
%   Name: rx_cell_index
%   Type: positive integer scalar
%   Description: estimated cell index onto which the mapping was applied.
%   It takes values between 0 and T-1.
%
%   Name: rx_symbols
%   Type: 2(T-1)x1 positive integer vector
%   Description: column vector with 2(T-1) coordinates. Each element 
%   represents the constellation point decided for each coordinate.
%   
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

[T,~] = size(Y); % T --> coherence time
Aj = 1/2^(B+1):1/2^B:(2^(B+1)-1)/2^(B+1); % grids
num_gridpoints = length(Aj); % number of points in each grid

% DENOISING
% [1, Sec. IV A) Step 1)]
[U,~,~] = svd(Y,'econ');
u = U(:,1); % left eigenvector corresponding to largest eigenvalue

% CELL AND LOCAL COORDINATES ESTIMATION
% [1, Sec. IV A) Step 2)]

% Cell estimation
[~,rx_cell_index] = max(abs(u));

% Inverse mapping
t = [u(1:rx_cell_index-1);u(rx_cell_index+1:end)] / u(rx_cell_index); % [1, Eq. (15)]
w = sqrt(2*log((1+abs(t).^2)./(1-abs(t).^2))) .* t ./ abs(t); % [1, Eq. (14)]

real_w = real(w);
imag_w = imag(w);

real_a = normcdf(real_w); 
imag_a = normcdf(imag_w);

a_aux = [real_a';imag_a'];
a = a_aux(:); % received points on each grid

% Symbols estimation
A = repmat(a,1,num_gridpoints)-repmat(Aj,2*(T-1),1);
A = A.^2;
[~,mininda] = min(A,[],2);
rx_symbols = mininda;

end