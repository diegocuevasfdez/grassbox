%-------------------------------------------------------------------------%
%
% Script:
% SER and BER performance of single-antenna Grass-Lattice [1-2], 
% Exp-Map [3] and Cube-Split [4] constellations. 
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
% [3] I. Kammoun, A. M. Cipriano and J. -C. Belfiore, "Non-Coherent Codes
% over the Grassmannian," in IEEE Transactions on Wireless Communications, 
% vol. 6, no. 10, pp. 3657-3667, Oct. 2007, doi: 10.1109/TWC.2007.06059
%
% [4] K. -H. Ngo, A. Decurninge, M. Guillaud and S. Yang, "Cube-Split: A
% Structured Grassmannian Constellation for Non-Coherent SIMO
% Communications," in IEEE Transactions on Wireless Communications, vol.
% 19, no. 3, pp. 1948-1964, Mar. 2020, doi: 10.1109/TWC.2019.2959781.
%
%-------------------------------------------------------------------------%
%
% Authors: Diego Cuevas Fernandez
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

clc
clear
format compact

addpath('./functions')

%% PARAMETERS
M = 1; % number of transmit antennas
N = 1; % number of receive antennas
T = 2; % coherence time

%% GRASS-LATTICE PARAMETERS
B_GL = 2; % number of bits per real dimension
P_GL = 2^B_GL; % P = 2^B points equispaced btw [alpha, 1-alpha]

NumCodewordsGL = 2^(2*(T-M)*M*B_GL); % constellation size Grass-Lattice
NbitsGL = 2*(T-M)*M*B_GL; % num bits per codeword
SpecEffGL = NbitsGL / T; % spectral efficiency in bits/sec/Hz

%% EXP-MAP PARAMETERS
B_expmap = 2 * B_GL; % bits/lattice symbol (same spectral efficiency as Grass-Lattice)
Q = 2^B_expmap; % Q-QAM (coherent code)

%% CUBE-SPLIT PARAMETERS
B_CS = B_GL; % number of bits per grid

NumCodewordsCS = T * 2^(2*(T-M)*M*B_CS); % constellation size Cube-Split
NbitsCS = log2(T) + 2*(T-M)*M*B_CS; % num bits per codeword
SpecEffCS = NbitsCS / T; % spectral efficiency in bits/sec/Hz

%% SIMULATION PARAMETERS
SNRdB = 0:5:30; % signal-to-noise ratio [dB]
snr = 10.^(SNRdB/10); % snr in linear scale

NoiseVar = M./(T*snr); % noise variance

EbN0_GL = 10*log10(snr./SpecEffGL); % Eb/N0 [dB] Grass-Lattice
EbN0_CS = 10*log10(snr./SpecEffCS); % Eb/N0 [dB] Cube-Split

NumSim = 1e4*ones(size(SNRdB)); % number of MC simulations
for cc = 1:length(SNRdB)
    if SNRdB(cc) >=15
        NumSim(cc) = 5*1e4;
    end
end

%% GRASS-LATTICE CONSTELLATION
params = struct();
params.SNR = 15; % in dB
params.NumSim = 1e3;
params.Alpha = 0.01:0.01:0.25;
params.PlotResults = 0;
params.ErrorRate = 'BER';

alpha = FindAlphaOpt(T,M,N,B_GL,params); % find optimal value of parameter 
                                         % alpha (some tables with the optimal
                                         % values of alpha for different 
                                         % scenarios are provided in [2])

latticeGL = alpha + (0:P_GL-1)*(1-2*alpha)./(P_GL-1); % lattice used 
                                                      % for each real 
                                                      % component ([2,Eq.7])

%% SIMULATION
SER_GL = zeros(1,length(SNRdB)); % symbol-error rate Grass-Lattice
SER_ExpMap = zeros(1,length(SNRdB)); % bit-error-rate Exp-Map
SER_CS = zeros(1,length(SNRdB)); % bit-error-rate Cube-Split
BER_GL = zeros(1,length(SNRdB)); % bit-error rate Grass-Lattice
BER_ExpMap = zeros(1,length(SNRdB)); % bit-error rate Exp-Map
BER_CS = zeros(1,length(SNRdB)); % bit-error-rate CS

for cc = 1:length(SNRdB) % SNR loop

    disp(['SNR =  ' int2str(SNRdB(cc))])
    
    NoiseVarIter = NoiseVar(cc);
    SymbolsError_GL = 0; % symbol error counter Grass-Lattice
    SymbolsError_ExpMap = 0; % symbol error counter Exp-Map
    SymbolsError_CS = 0; % symbol error counter Cube-Split
    BitsError_GL = 0; % bit error counter Grass-Lattice
    BitsError_ExpMap = 0; % bit error counter Exp-Map
    BitsError_CS = 0; % bit error counter Cube-Split

    for ss = 1:NumSim(cc) % MC loop (could be paralelized using "parfor")

        if ~mod(ss,floor(NumSim(cc)/10)), fprintf('.'); end

        % Encoder and Mapper (Grass-Lattice)
        tx_bits_GL_expmap = randi([0,1],[2*M*(T-M)*B_GL,1]); % bit string to transmit 
                                                           % (Gray encoded)
        tx_bitsbin_GL = Gray2Bin(tx_bits_GL_expmap,B_GL);
        tx_symbols_GL = bit2int(tx_bitsbin_GL,B_GL)+1; % bit-to-symbol mapper
        X_GL = GrassLatticeEncoding(M,latticeGL,tx_symbols_GL); % TX signal

        % Encoder and Mapper (Exp-Map)
        tx_symbols_expmap = bit2int(tx_bits_GL_expmap,B_expmap); % bit-to-symbol mapper
        X_expmap = ExpMapEncoding(M,Q,tx_symbols_expmap); % TX signal

        % Encoder and Mapper (Cube-Split)
        tx_bits_cell_index_CS = randi([0,1],[log2(T),1]); % cell bits
        tx_cell_index_CS = bit2int(tx_bits_cell_index_CS,log2(T)) + 1;

        tx_bits_symbols_CS = tx_bits_GL_expmap; % coordinate bits
        tx_bitsbin_CS = Gray2Bin(tx_bits_symbols_CS,B_CS);
        tx_symbols_CS = bit2int(tx_bitsbin_CS,B_CS) + 1; % bit-to-symbol mapper

        X_CS = CubeSplitEncoding(B_CS,tx_cell_index_CS,tx_symbols_CS); % TX signal

        % Rayleigh MIMO Channel (fixed during the coherence block)
        H = randn(M,N)/sqrt(2) + 1i*randn(M,N)/sqrt(2);

        % AWGN
        Noise = sqrt(NoiseVarIter/2)*randn(T,N) + 1i*sqrt(NoiseVarIter/2)*randn(T,N);

        % RX Signal
        Y_GL = X_GL * H + Noise;
        Y_expmap = X_expmap * H + Noise;
        Y_CS = X_CS * H + Noise;

        % Decoder and Demapper (Grass-Lattice)
        rx_symbols_GL = GrassLatticeDecoding(M,latticeGL,Y_GL);
        rx_bitsbin_GL = int2bit(rx_symbols_GL-1,B_GL);
        rx_bits_GL = Bin2Gray(rx_bitsbin_GL,B_GL);

        % Decoder and Demapper (Exp-Map)
        rx_symbols_expmap = ExpMapDecoding(M,Q,Y_expmap);
        rx_bits_expmap = int2bit(rx_symbols_expmap,B_expmap);

        % Decoder and Demapper (Cube-Split)
        [rx_cell_index_CS, rx_symbols_CS] = CubeSplitDecoding(B_CS,Y_CS);
        rx_bits_cell_index_CS = int2bit(rx_cell_index_CS-1,log2(T)); % estimated cell bits
        rx_bitsbin_CS = int2bit(rx_symbols_CS-1,B_CS); 
        rx_bits_symbols_CS = Bin2Gray(rx_bitsbin_CS,B_CS); % estimated coordinate bits

        % Symbol error counters
        SymbolsError_GL = SymbolsError_GL + any(rx_symbols_GL~=tx_symbols_GL);
        SymbolsError_ExpMap = SymbolsError_ExpMap + any(rx_symbols_expmap~=tx_symbols_expmap);
        SymbolsError_CS = SymbolsError_CS + ((rx_cell_index_CS~=tx_cell_index_CS) || any(rx_symbols_CS~=tx_symbols_CS));

        % Bit error counters
        BitsError_GL = BitsError_GL + sum(rx_bits_GL~=tx_bits_GL_expmap);
        BitsError_ExpMap = BitsError_ExpMap + sum(rx_bits_expmap~=tx_bits_GL_expmap);
        BitsError_CS = BitsError_CS + sum(rx_bits_cell_index_CS~=tx_bits_cell_index_CS) + sum(rx_bits_symbols_CS~=tx_bits_symbols_CS);

    end

    % SER computation
    SER_GL(cc) = SymbolsError_GL / NumSim(cc);
    SER_ExpMap(cc) = SymbolsError_ExpMap / NumSim(cc);
    SER_CS(cc) = SymbolsError_CS / NumSim(cc);

    % BER computation
    BER_GL(cc) = BitsError_GL / (NumSim(cc)*NbitsGL);
    BER_ExpMap(cc) = BitsError_ExpMap / (NumSim(cc)*NbitsGL); % same spectral efficiency as Grass-Lattice
    BER_CS(cc) = BitsError_CS / (NumSim(cc)*NbitsCS);

    fprintf('\n')

end


%% PLOT RESULTS
fs = 11;
lw = 1.5;
ms = 8;

% SER curve
figure(1);clf;semilogy(EbN0_GL,SER_GL,'r-^','MarkerSize',ms,'LineWidth',lw);hold on;
semilogy(EbN0_GL,SER_ExpMap,'g-s','MarkerSize',ms,'LineWidth',lw);
semilogy(EbN0_CS,SER_CS,'b-o','MarkerSize',ms,'LineWidth',lw);
xlabel('E_b/N_0 (dB)'); ylabel('SER');
title(['T = ' int2str(T) ', ' 'M = ' int2str(M) ', ' 'N = ' int2str(N) ', ' '\eta_{GL} = ' int2str(SpecEffGL) ', ' '\eta_{CS} = ' int2str(SpecEffCS)])
legend('Grass-Lattice','Exp-Map','Cube-Split')
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
grid on

% BER curve
figure(2);clf;semilogy(EbN0_GL,BER_GL,'r-^','MarkerSize',ms,'LineWidth',lw); hold on;
semilogy(EbN0_GL,BER_ExpMap,'g-s','MarkerSize',ms,'LineWidth',lw);
semilogy(EbN0_CS,BER_CS,'b-o','MarkerSize',ms,'LineWidth',lw);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title(['T = ' num2str(T) ', ' 'M = ' num2str(M) ', ' 'N = ' num2str(N) ', ' '\eta_{GL} = ' int2str(SpecEffGL) ', ' '\eta_{CS} = ' int2str(SpecEffCS)])
legend('Grass-Lattice','Exp-Map','Cube-Split')
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
grid on

