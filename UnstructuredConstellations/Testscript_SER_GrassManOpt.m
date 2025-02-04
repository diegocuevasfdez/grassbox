%-------------------------------------------------------------------------%
%
% Script:
% SER performance of GrassManOpt [1] constellation.
%
% [1] D. Cuevas, C. Beltran, I. Santamaria, V. Tucek and G. Peters,
% "A Fast Algorithm for Designing Grassmannian Constellations," 25th
% International ITG Workshop on Smart Antennas, French Riviera, France,
% Nov. 2021.
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
addpath('./BestPackings')

%% CONSTELLATION PARAMETERS
NumCodewords = 16; % constellation size
M = 2; % transmit antennas
N = 2; % receive antennas
T = 4; % time slots > max(M,N);

%% GENERATE CONSTELLATION
% Optimization parameters
opt_params = struct();
opt_params.Alpha = 1.1;
opt_params.MaxIter = 200;
opt_params.hini = 0.1;
opt_params.hfin = 1e-4;
opt_params.MinImp = 1e-5;
opt_params.Distance = 'chordal';
opt_params.verbose = false;

% Starting point
C_ini = SubspacesGenerationRandom(T,M,NumCodewords); % random subspaces generation

% Optimization algorithm
[C,separation] = ManifoldOptimization(C_ini,opt_params);

% Load BestPackings
% Instead of generating a new constellation, we can load one
% Cbest_struct = load('Cbest4x2x16.mat');
% C = Cbest_struct.Cbest;

%% SIMULATION PARAMETERS
SNR = [0 3 6 9 12 15 18 21]; % signal-to-noise ratio [dB]
NumSim = [2e4 2e4 2e4 2e4 2e4 2e5 2e5 2e5]; % num channel realizations (coherence blocks)

NoiseVar = (M/T)*10.^(-SNR/10); % noise variance

%% SIMULATION    
SER = zeros(1,length(SNR)); % symbol error rate

for cc = 1:length(SNR) % SNR loop

    disp(['SNR =  ' int2str(SNR(cc))])
    NoiseVarIter = NoiseVar(cc);
    SymbolsError = 0; % symbol-error-rate counter

    for ss = 1:NumSim(cc) % MC loop (could be paralelized using "parfor")
        
        if ~mod(ss,floor(NumSim(cc)/10)), fprintf('.'); end
        
        % TX Signal
        TX_codeword = randi(NumCodewords);
        X = C(:,:,TX_codeword);
        
        % Rayleigh MIMO channel (fixed during the coherence block)
        H = randn(M,N)/sqrt(2) + 1i*randn(M,N)/sqrt(2); 
        
        % AWGN 
        Noise = sqrt(NoiseVarIter/2)*randn(T,N) + 1i*sqrt(NoiseVarIter/2)*randn(T,N);
        
        % RX Signal
        Y = X*H + Noise;
        
        % ML Detector
        RX_codeword = MLGrassDetector(C,Y);

        % Error counter
        SymbolsError = SymbolsError + (RX_codeword ~= TX_codeword);

    end
    
    SER(cc) = SymbolsError / NumSim(cc); % SER computation

    fprintf('\n')

end

%% PLOT RESULTS
fs = 11;
lw = 1.5;
ms = 8;

% SER curve
figure(1);clf;semilogy(SNR,SER,'r-^','MarkerSize',ms,'LineWidth',lw);hold on;
xlabel('SNR (dB)');
ylabel('SER');
legend('GrassManOpt');
title(['T = ' num2str(T) ', ' 'M = ' num2str(M) ', ' 'N = ' num2str(N) ', ' '|C| = ' int2str(NumCodewords)])
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
grid on
