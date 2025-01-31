%-------------------------------------------------------------------------%
%
% Script:
% SER and BER performance of multi-antenna Grass-Lattice [1] and Exp-Map [2] 
% constellations.
%
% [1] D. Cuevas, C. Beltran, M. Gutierrez, I. Santamaria and V. Tucek, 
% "Structured Multi-Antenna Grassmannian Constellations for Noncoherent
% Communications," IEEE 13rd Sensor Array and Multichannel Signal
% Processing Workshop (SAM), Corvallis, OR, USA, Jul. 2024,
% doi: 10.1109/SAM60225.2024.10636457
%
% [2] I. Kammoun, A. M. Cipriano and J. -C. Belfiore, "Non-Coherent Codes
% over the Grassmannian," in IEEE Transactions on Wireless Communications, 
% vol. 6, no. 10, pp. 3657-3667, Oct. 2007, doi: 10.1109/TWC.2007.06059
%
%-------------------------------------------------------------------------%
%
% Author: Diego Cuevas Fernandez
% Department: Advanced Signal Processing Group (GTAS)
% University: Universidad de Cantabria, Spain
% Year: 2025
%
%-------------------------------------------------------------------------%
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/
%
%-------------------------------------------------------------------------%

clear
close all
clc

addpath('./functions')

%% PARAMETERS
T = 6; % time slots > max(M,N);
M = 2; % transmit antennas
N = 4; % receive antennas

%% GRASS-LATTICE PARAMETERS
B_GL = 2; % bits/lattice symbol
NumCodewords = 2^(2*M*(T-M)*B_GL); % constellation size
Nbits = 2*M*(T-M)*B_GL; % num bits per constellation point
SpecEff = Nbits/T; % spectral efficiency in bits/sec/Hz

params = struct();
params.SNR = 15; % in dB
params.NumSim = 1e3;
params.Alpha = 0.01:0.02:0.25;
params.PlotResults = 0;
params.ErrorRate = 'BER';

alpha = FindAlphaOpt(T,M,N,B_GL,params); % find optimal value of parameter 
                                         % alpha (some tables with the optimal
                                         % values of alpha for different 
                                         % scenarios are provided in [3])

P_GL = 2^B_GL; % P points equispaced btw [alpha, 1-alpha]
lattice_GL = alpha + (0:P_GL-1)*(1-2*alpha)./(P_GL-1); % lattice used 
                                                       % for each real 
                                                       % component ([3,Eq.10])

%% EXP-MAP PARAMETERS
B_expmap = 2 * B_GL; % bits/lattice symbol (same spectral efficiency as Grass-Lattice)
Q = 2^B_expmap; % Q-QAM (Coherent code)

%% SIMULATION PARAMETERS
SNRdB = 10:5:30; % SNR in dBs
snr = 10.^(SNRdB/10); % snr in linear scale

NoiseVar = M./(T*snr); % noise variance

NumSim = 1e3*ones(size(SNRdB)); % Number of MC simulations
for cc = 1:length(SNRdB)
    if SNRdB(cc) >=20
        NumSim(cc) = 1*1e4;
    end
    if SNRdB(cc) >=30
        NumSim(cc) = 5*1e4;
    end
    if SNRdB(cc) >=40
        NumSim(cc) = 1e5;
    end
end

%% SIMULATIONS 
SER_GL = zeros(1,length(SNRdB)); % symbol-error rate Grass-Lattice
SER_ExpMap = zeros(1,length(SNRdB)); % symbol error rate Exp-Map
BER_GL = zeros(1,length(SNRdB)); % bit-error-rate Grass-Lattice
BER_ExpMap = zeros(1,length(SNRdB)); % bit error rate Exp-Map

for cc = 1:length(SNRdB) % SNR loop

    disp(['SNR =  ' int2str(SNRdB(cc))])

    NoiseVarIter = NoiseVar(cc);
    SymbolsError_GL = 0; % initialize symbol error counter Grass-Lattice
    BitsError_GL = 0; % initialize bit error counter Grass-Lattice
    SymbolsError_ExpMap = 0; % initialize symbol error counter Exp-Map
    BitsError_ExpMap = 0; % initialize bit error counter Exp-Map

    for ss = 1:NumSim(cc) % MC loop (could be paralelized using "parfor")

        if ~mod(ss,floor(NumSim(cc)/10)), fprintf('.'); end

        % Bit string to transmit
        tx_bits = randi([0,1],Nbits,1);

        % TX Signal (Grass-Lattice)
        tx_bitsbin = Gray2Bin(tx_bits,B_GL);
        tx_symbols_GL = bit2int(tx_bitsbin,B_GL)+1; % bit-to-symbol mapper
        X_GL = GrassLatticeEncoding(M,lattice_GL,tx_symbols_GL);

        % TX Signal (Exp-Map)
        tx_symbols_expmap = bit2int(tx_bits,B_expmap); % bit-to-symbol mapper
        X_expmap = ExpMapEncoding(M,Q,tx_symbols_expmap); % TX signal

        % Rayleigh MIMO channel (fixed during the coherence block)
        H = randn(M,N)/sqrt(2) + 1i*randn(M,N)/sqrt(2);
     
        % AWGN
        Noise = sqrt(NoiseVarIter/2)*randn(T,N) + 1i*sqrt(NoiseVarIter/2)*randn(T,N);

        % RX Signal
        Y_GL = X_GL * H + Noise;
        Y_expmap = X_expmap * H + Noise;

        % Decoder and Demapper (Grass-Lattice)
        rx_symbols_GL = GrassLatticeDecoding(M,lattice_GL,Y_GL);
        rx_bits_GL_bin = int2bit(rx_symbols_GL-1,B_GL);
        rx_bits_GL = Bin2Gray(rx_bits_GL_bin, B_GL);

        % Decoder and Demapper (Exp-Map)
        rx_symbols_expmap = ExpMapDecoding(M,Q,Y_expmap);
        rx_bits_expmap = int2bit(rx_symbols_expmap,B_expmap);

        % Symbol error counters
        SymbolsError_GL = SymbolsError_GL + (1-all(rx_symbols_GL==tx_symbols_GL));
        SymbolsError_ExpMap = SymbolsError_ExpMap + (1-all(rx_symbols_expmap==tx_symbols_expmap));
        
        % Bit error counters
        BitsError_GL = BitsError_GL + sum(rx_bits_GL~=tx_bits);
        BitsError_ExpMap = BitsError_ExpMap + sum(rx_bits_expmap~=tx_bits);
        
    end
    
    % SER computation
    SER_GL(cc) = SymbolsError_GL / NumSim(cc);
    SER_ExpMap(cc) = SymbolsError_ExpMap / NumSim(cc);

    % BER computation
    BER_GL(cc) = BitsError_GL / (Nbits * NumSim(cc));
    BER_ExpMap(cc) = BitsError_ExpMap / (Nbits * NumSim(cc));

    fprintf('\n')

end


%% PLOT RESULTS
fs = 11;
lw = 1.5;
ms = 8;

% SER curve
figure(1);clf;semilogy(SNRdB,SER_GL,'r-^','MarkerSize',ms,'LineWidth',lw);hold on;
semilogy(SNRdB,SER_ExpMap,'b-o','MarkerSize',ms,'LineWidth',lw)
xlabel('SNR (dB)');
ylabel('SER');
title(['T = ', num2str(T), ', M = ', num2str(M), ', N = ', num2str(N), ', \eta = ', num2str(SpecEff)])
legend('Grass-Lattice','Exp-Map')
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
grid on

% BER curve
figure(2);clf;semilogy(SNRdB,BER_GL,'r-^','MarkerSize',ms,'LineWidth',lw);hold on;
semilogy(SNRdB,BER_ExpMap,'b-o','MarkerSize',ms,'LineWidth',lw)
xlabel('SNR (dB)');
ylabel('BER');
title(['T = ', num2str(T), ', M = ', num2str(M), ', N = ', num2str(N), ', \eta = ', num2str(SpecEff)])
legend('Grass-Lattice','Exp-Map')
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
grid on