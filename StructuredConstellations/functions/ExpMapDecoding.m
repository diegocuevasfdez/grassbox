function rx_symbols = ExpMapDecoding(M, Q, Y)

%EXPMAPDECODING Decoding for the Exp-Map constellation. This
%constellation design was proposed in the following paper:
% 
% [1] I. Kammoun, A. M. Cipriano and J. -C. Belfiore, "Non-Coherent Codes
% over the Grassmannian," in IEEE Transactions on Wireless Communications, 
% vol. 6, no. 10, pp. 3657-3667, Oct. 2007, doi: 10.1109/TWC.2007.06059
%
%   rx_symbols = EXPMAPDECODING(M, Q, Y) estimates the transmitted symbols 
%   belonging to a QAM constellation of modulation order "Q" from the 
%   received signal "Y" and for the chosen value of transmit antennas "M".
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
%   Name: Y
%   Type: TxN double precision vector
%   Description: TxN received signal, where T is the coherence time and N
%   the number of receive antennas.
%
%   OUTPUT: 
%   Name: rx_symbols
%   Type: 2M(T-M)x1 positive integer vector
%   Description: column vector with 2M(T-M) coordinates. Each element 
%   represents the constellation point decided for each coordinate.
%   
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

[T,N] = size(Y); % T --> coherence time, N --> number of receive antennas

if M == 1

    % Homothetic factor [1, Eq. (25)]
    alpha = pi / (2 * (sqrt(2*(T-1)) * (sqrt(Q)-1) + 1));
    
    % SIMPLIFIED DECODING [1, Sec. IV B)]

    % Compute correlation matrix
    Ry = 1/N * (Y * Y');

    % Extract eigenvector corresponding to greatest eigenvalue of Ry
    [V,D] = eig(Ry);
    [~,ind] = sort(diag(D),'descend');
    xc = V(:,ind(1));

    % Obtain normalized received signal in CS form
    zc = xc / norm(xc);
    z = exp(-1j * angle(zc(1))) * zc;
    
    % Apply inverse exponential map [1, Eq. (9)]
    vz = -acos(z(1))/sqrt(1-z(1)^2) * z(2:end);
    vz_alpha = vz / alpha;
    
    % Estimate received symbols
    rx_symbols = qamdemod(vz_alpha,Q);

elseif M == 2

    % Homothetic factor [1, Sec. V C)]

    if T == 4

        alpha = 0.3;

    elseif T == 6

        alpha = 0.23;

    else

        error('Coherence time not supported')

    end

    % SIMPLIFIED DECODING [1, Sec. V B)]
    
    % Obtain subspace rough estimation of the sent codeword
    [U,~,~] = svd(Y,'econ');
    Yc = U(:,1:M);
    Yc1 = Yc(1:M,1:M);

    % Obtain from basis Yc1 the corresponding basis Xe in CS form
    [Pe, ~, Q1] = svd(Yc1,'econ');
    Xe = Yc*Q1*Pe'; % estimated RX codeword

    % Channel estimation
    He = Xe'*Y;
    
    % Estimate coherent codewords
    YHe = Y/He; % Y * inv(He)
    Ve = -YHe(M+1:end,:)';
    Ve_t = Ve.';
    s_rx = 1/alpha * Ve_t(:); % received QAM symbols

    % Estimate QAM symbols
    rx_symbols = qamdemod(s_rx,Q,'UnitAveragePower',true);

else

    error('Number of transmit antennas/coherence time not supported');

end

end