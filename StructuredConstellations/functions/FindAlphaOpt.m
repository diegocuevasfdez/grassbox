function alpha_opt = FindAlphaOpt(T,M,N,B,varargin)

%COMPUTE_ALPHA_OPT Computation of best alpha parameter in terms of symbol
%error rate (SER) or bit error rate (BER) for structured Grassmannian 
%constellation design method 'Grass-Lattice'. This constellation design was
%proposed in the following papers:
%
% [1] D. Cuevas, J. Alvarez-Vizoso, C. Beltran, I. Santamaria, V. Tucek and
% G. Peters, "A Measure Preserving Mapping for Structured Grassmannian 
% Constellations in SIMO Channels," IEEE Global Communications Conference
% (GLOBECOM), Rio de Janeiro, Brazil, Dec. 2022, 
% doi: 10.1109/GLOBECOM48099.2022.10000914
%
% [2] D. Cuevas, J. Alvarez-Vizoso, C. Beltran, I. Santamaria, V. Tucek and
% G. Peters, "Constellations on the Sphere with Efficient Encoding-Decoding
% for Noncoherent Communications," in IEEE Transactions on Wireless
% Communications, vol. 23, no. 3, pp. 1886-1898, Mar. 2024, 
% doi: 10.1109/TWC.2023.3292935
%
% [3] D. Cuevas, C. Beltran, M. Gutierrez, I. Santamaria and V. Tucek, 
% "Structured Multi-Antenna Grassmannian Constellations for Noncoherent
% Communications," IEEE 13rd Sensor Array and Multichannel Signal
% Processing Workshop (SAM), Corvallis, OR, USA, Jul. 2024,
% doi: 10.1109/SAM60225.2024.10636457
%
%   ALPHA_OPT = COMPUTE_ALPHA_OPT(T,M,N,B) computes best alpha parameter 
%   in terms of SER or BER for a given coherence time T, number of transmit 
%   antennas M, number of receive antennas N and number of bits used to 
%   codify each real component in Grass-Lattice method.
%
%   [ ... ] = MANIFOLD_OPTIMIZATION(..., params) specifies
%   optional parameter name/value pairs to control the SER/BER simulation
%   performed by COMPUTE_ALPHA_OPT. Parameters are:
%
%   'SNR' - Signal-to-Noise Ratio for which we compute the SER/BER. 
%   Default is 20 dB.
%
%   'NumSim' - Number of iterations for each alpha value in the SER/BER 
%   simulation. Default is 1e5.
%
%   'Alpha' - Range of values of the alpha parameter. Default is 
%   [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45] (values must be in the 
%   interval (0,0.5)).
%
%   'PlotResults' - If 'True', plots the curve SER/BER vs. Alpha.
%   Default is 'False'.
%
%   'ErrorRate' - If 'SER', use symbol error rate as figure of merit. If
%   'BER', use bit error rate as figure of merit. Default is 'BER'.
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

% Default values
params = struct();
params.SNR = 15; % in dB
params.NumSim = 1e5;
params.Alpha = 0.01:0.01:0.25;
params.PlotResults = 0;
params.ErrorRate = 'BER';

if nargin < 4
    error(message('TooFewInputs'));
elseif nargin == 5
    params_aux = varargin{1};
    for arg = fieldnames(params_aux)'
        parameter = arg{1};
        param_value = params_aux.(parameter);
        switch parameter
            case 'SNR'
                params.SNR = param_value;
            case 'NumSim'
                params.NumSim = param_value;
            case 'Alpha'
                params.Alpha = param_value;
            case 'PlotResults'
                params.PlotResults = param_value;
            case 'ErrorRate'
                params.ErrorRate = param_value;
        end
    end
elseif nargin > 5
    error('TooManyInputs');
end

% Simulation parameters
eff = 2*M*(T-M)*B; % Spectral efficiency
VarNoise = (M/T)*10.^(-params.SNR/10); % Noise variance

% Grass-Lattice parameters
P = 2^B; % P = 2^B points equispaced btw [alpha, 1-alpha]

% Simulations
SER_alpha = zeros(size(params.Alpha,1));
BER_alpha = zeros(size(params.Alpha,1));

for iter_alpha = 1:length(params.Alpha)

    alpha_i = params.Alpha(iter_alpha);
    latticeGL = alpha_i + (0:P-1)*(1-2*alpha_i)/(P-1);

    fprintf("Alpha = %.2f",alpha_i)
    fprintf("\n")

    counter_symbols_GL = 0;
    counter_bits_GL = 0;

    for ss = 1:params.NumSim

        if ~mod(ss,floor(params.NumSim/10)), fprintf('.'); end

        % Encoder and Mapper
        TX_bits_GL = randi([0,1],[2*M*(T-M)*B,1]); % bit string to transmit (Gray encoded)
        TX_bitsbin = Gray2Bin(TX_bits_GL,B);
        TX_symbols_GL = bit2int(TX_bitsbin,B)+1; % bit-to-symbol mapper
        X_GL = GrassLatticeEncoding(M,latticeGL,TX_symbols_GL);

        % Rayleigh MIMO channel (fixed during the coherence block)
        H = randn(M,N)/sqrt(2) + 1i*randn(M,N)/sqrt(2);

        % AWGN
        Noise = sqrt(VarNoise/2)*randn(T,N) + 1i*sqrt(VarNoise/2)*randn(T,N);

        % RX Signal
        Y_GL = X_GL*H + Noise;

        % Decoder and Demapper
        RX_symbols_GL = GrassLatticeDecoding(M,latticeGL,Y_GL);
        RX_bitsbin = int2bit(RX_symbols_GL-1,B);
        RX_bits_GL = Bin2Gray(RX_bitsbin, B);
        counter_symbols_GL = counter_symbols_GL + (1-all(RX_symbols_GL==TX_symbols_GL));
        counter_bits_GL = counter_bits_GL + sum(TX_bits_GL~=RX_bits_GL);
    end

    SER_alpha(iter_alpha) = counter_symbols_GL / params.NumSim;
    BER_alpha(iter_alpha) = counter_bits_GL / (eff * params.NumSim);
    
    fprintf("\n")

end

if strcmp(params.ErrorRate,'SER')
    [~,ind_alpha_opt] = min(SER_alpha);
elseif strcmp(params.ErrorRate,'BER')
    [~,ind_alpha_opt] = min(BER_alpha);
end

alpha_opt = params.Alpha(ind_alpha_opt);
fprintf("Optimum alpha = %.2f\n\n",alpha_opt);

if params.PlotResults
    % Plots
    lw = 1.5;
    ms = 8;
    
    if strcmp(params.ErrorRate,'SER')
        figure
        semilogy(params.Alpha,SER_alpha,'r-s','LineWidth',lw,'MarkerSize',ms)
        grid on
        xlabel('\alpha')
        ylabel('SER')
        title(['SNR = ',num2str(params.SNR), ' dB, B = ', num2str(B), ', T = ', num2str(T), ', M = ', num2str(M), ', N = ', num2str(N)])
    end

    if strcmp(params.ErrorRate,'BER')
        figure
        semilogy(params.Alpha,BER_alpha,'b-s','LineWidth',lw,'MarkerSize',ms)
        grid on
        xlabel('\alpha')
        ylabel('BER')
        title(['SNR = ',num2str(params.SNR), ' dB, B = ', num2str(B), ', T = ', num2str(T), ', M = ', num2str(M), ', N = ', num2str(N)])
    end
end

end