%-------------------------------------------------------------------------%
% Script:
% SER performance of GrassManOpt constellation.
%
% Author: Diego Cuevas Fernandez
% Department: Advanced Signal Processing Group (GTAS)
% University: Universidad de Cantabria, Spain
% Year: 2025
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/
%-------------------------------------------------------------------------%

clear
close all
clc

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
C_ini = subspaces_generation_random(T,M,NumCodewords); % random subspaces generation

% Optimization algorithm
[C,separation] = manifold_optimization(C_ini,opt_params);

%% SIMULATION PARAMETERS
SNR = [0 3 6 9 12 15 18 21]; % signal-to-noise ratio [dB]
num_sim = [2e4 2e4 2e4 2e4 2e4 2e5 2e5 2e5];   % num channel realizations (coherence blocks)

NoiseVar = (M/T)*10.^(-SNR/10); % noise variance

SER = zeros(1,length(SNR)); % symbol error rate

%% SIMULATION    
for cc = 1:length(SNR)
    disp(['SNR =  ' int2str(SNR(cc))])
    VarNoise = NoiseVar(cc);
    for ss = 1:num_sim(cc)
        
        if ~mod(ss,floor(num_sim(cc)/10)), fprintf('.'); end
        
        % TX Signal
        TX_Codeword = randi(NumCodewords);
        X = C(:,:,TX_Codeword);
        
        % Rayleigh MIMO channel (fixed during the coherence block)
        H = randn(M,N)/sqrt(2) + 1i*randn(M,N)/sqrt(2); 
        
        % AWGN 
        Noise = sqrt(VarNoise/2)*randn(T,N) + 1i*sqrt(VarNoise/2)*randn(T,N);
        
        % RX Signal
        Y = X*H + Noise;
        
        % ML Detector
        RX_codeword = MLGrassDetector(C,Y);

        % Errors counter
        SER(cc) = SER(cc) + (RX_codeword ~= TX_Codeword);
    end
    fprintf('\n')
end

SER = SER ./ num_sim; % SER computation

%% PLOT RESULTS
lw = 1.5;
ms = 8;
figure(1);clf;semilogy(SNR,SER,'r-^','MarkerSize',ms,'LineWidth',lw);hold on;
xlabel('SNR (dB)');
ylabel('SER');
legend('GrassManOpt');
grid on;
hold off;