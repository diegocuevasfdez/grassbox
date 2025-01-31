function th = thnt(n,t)

%THNT Implementation of function "t * h_n(t)" used in mapping vartheta2
%for multi-antenna Grass-Lattice encoding. This function is proposed in
%[1, Lemma 1].
%
% [1] D. Cuevas, C. Beltran, M. Gutierrez, I. Santamaria and V. Tucek, 
% "Structured Multi-Antenna Grassmannian Constellations for Noncoherent
% Communications," IEEE 13rd Sensor Array and Multichannel Signal
% Processing Workshop (SAM), Corvallis, OR, USA, Jul. 2024,
% doi: 10.1109/SAM60225.2024.10636457
%
%   TH = THNT(N, T) computes the value of function "t * h_n(t)".
%
%   INPUTS:
%   Name: n
%   Type: positive integer scalar
%   Description: order of function to be evaluated ("n" parameter in 
%   [1, Lemma 1]).
%
%   Name: t
%   Type: double precision scalar
%   Description: argument of function to be evaluated ("t" parameter in 
%   [1, Lemma 1].
%
%   OUTPUT
%   Name: th
%   Type: double precision scalar in [0,1]
%   Description: value of "t * h_n(t)". 
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

aux=@(th) (n+1)*(th).^(2*n)-n*(th).^(2*n+2)-1+exp(-t.^2)*sum(t.^(0:2:2*n-2)./factorial(0:n-1));
th=fzero(aux,[0,1]);

end