function X = ExpMapEncoding(M, Q, tx_symbols)

%EXPMAPENCODING Encoding for the Exp-Map constellation. This
%constellation design was proposed in the following paper:
% 
% [1] I. Kammoun, A. M. Cipriano and J. -C. Belfiore, "Non-Coherent Codes
% over the Grassmannian," in IEEE Transactions on Wireless Communications, 
% vol. 6, no. 10, pp. 3657-3667, Oct. 2007, doi: 10.1109/TWC.2007.06059
%
%   X = EXPMAPENCODING(M, Q, TX_SYMBOLS) generates the 
%   Exp-Map encoded codeword for "M" transmit antennas and for the desired 
%   "tx_symbols" and modulation order "Q".
%
%   INPUTS:
%   Name: M
%   Type: positive integer scalar
%   Description: number of transmit antennas.
%
%   Name: Q
%   Type: positive integer scalar
%   Description: QAM modulation order of coherent symbols.
%
%   Name: tx_symbols
%   Type: 2M(T-M)x1 positive integer vector
%   Description: column vector with 2M(T-M) coordinates. Each element 
%   represents the constellation point for one coordinate.
%
%   OUTPUT:
%   Name: X
%   Type: TxM double precision vector
%   Description: encoded point in Gr(M,C^T). 
%   
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

T = length(tx_symbols)/M + M; % coherence time

if M == 1

    % Q-QAM modulator
    s = qammod(tx_symbols,Q);
    
    % Homothetic factor [1, Eq. (25)]
    alpha = pi / (2 * (sqrt(2*(T-1)) * (sqrt(Q)-1) + 1));

    % Coherent codeword
    v = alpha * s;
    
    % Non-coherent codeword [1, Eq. (8)]
    rho = norm(v);
    X = [cos(rho) ; -sin(rho)/rho * v]; % TX signal

elseif M == 2

    % Q-QAM modulator
    s = qammod(tx_symbols,Q,'UnitAveragePower',true);

    % Homothetic factor [1, Sec. V C)]
    if T == 4

        alpha = 0.3;

    elseif T == 6

        alpha = 0.23;

    else

        error('Coherence time not supported')

    end
    
    % Coherent codeword
    s_aux = reshape(s,T-M,M).';
    V = alpha * s_aux;
    
    % Non-coherent codeword
    [P,Sigma,Q] = svd(V,'econ');
    lambdas = diag(Sigma);
    C = diag(cos(lambdas));
    S = diag(sin(lambdas));
    X = [P*C*P'; -Q*S*P']; % TX signal

else

    error('Number of transmit antennas not supported');

end

end