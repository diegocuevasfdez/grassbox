function rx_symbols = GrassLatticeDecoding(M, lattice, Y)

%GRASSLATTICEDECODING Decoding for the Grass-Lattice constellation. This
%constellation design was proposed in the following papers:
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
%   RX_SYMBOLS = GRASSLATTICEDECODING(M, LATTICE, Y) estimates the 
%   transmitted symbols from the received signal "Y" on a given "lattice"
%   and for the chosen value of transmit antennas "M".
%
%   INPUTS:
%   Name: M
%   Type: positive integer scalar
%   Description: number of transmit antennas. Grass-Lattice mapping has
%   been derived for either 1 or 2 antennas at the transmitter.
%
%   Name: lattice
%   Type: 1x2M(T-M) double precision row vector
%   Description: lattice for the unit cube coordinates. If the lattice
%   has P points, the code uses B=log2(P) bits per coordinate. We decode a
%   total of 2M(T-M) coordinates, where T is the coherence time.
%
%   Name: Y
%   Type: TxN double precision matrix
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

[T,~] = size(Y);  % T is the coherence time
 
% DENOISING STEP
[U,~,~] = svd(Y);
C = U(:,1:M);

if M == 1

    vr = C; 
    
    % INVERT MAPPING VARTHETA3
    
    % [2, Sec. IV B) Step 1]
    w = vr(2:end)*abs(vr(1))/(vr(1)*norm(vr));
    
    % INVERT MAPPING VARTHETA2
    
    % [2, Sec. IV B) Step 2]
    nw = norm(w);
    sr = sqrt(icdf('chi2',(nw.^(2*(T-1))), 2*(T-1))/2);
    zr = sr*w/nw;
    
    % INVERT MAPPING VARTHETA1
    
    % [2, Sec. IV B) Step 3]
    ar = real(zr);
    br = imag(zr);
    ahat = cdf('normal',ar,0,1/sqrt(2));
    bhat = cdf('normal',br,0,1/sqrt(2));
    
    % HARD DECISIONS

    % [2, Sec. IV B) Step 4]
    A = repmat(ahat,1,length(lattice))-repmat(lattice,T-1,1);
    A = A.^2;
    [~,mininda] = min(A,[],2);

    B = repmat(bhat,1,length(lattice))-repmat(lattice,T-1,1);
    B = B.^2;
    [~,minindb] = min(B,[],2);

    rx_symbols = zeros(2*(T-1),1);
    rx_symbols(1:2:end) = mininda;
    rx_symbols(2:2:end) = minindb;

elseif M == 2
    
    % INVERT MAPPING VARTHETA3
    
    % [3, Sec. IV B) Step 2]
    C1 = C(1:M,1:M); % upper square block
    C2 = C(M+1:end,1:M); % lower block
    [u,~,v] = svd(C1);
    Q = u*v';
    W = C2*Q';
    
    % INVERT MAPPING VARTHETA2

    % [3, Sec. IV B) Step 3]
    p = W(:,1);
    t = W(:,2);
    v = p'*t/(norm(p)^2)*p;
    qbar = t-v;
    
    % [3, Sec. IV B) Step 4]
    uu = fzero(@(uu) (1-exp(-uu^2))^(1/2)...
        - norm(v)/sqrt((1-norm(p)^2)*(1-norm(qbar)^2)),[0,100]);
    u = uu*v/norm(v);
    
    % [3, Sec. IV B) Step 5]
    m = fzero(@(m) thnt(T-2,m)-norm(p),[0,100]);
    n = fzero(@(n) thnt(T-3,n)-norm(qbar),[0,100]);
    r = m*p/norm(p);
    sbar = n*qbar/norm(qbar);

    % [3, Sec. IV B) Step 6]
    s = sbar+u;
    
    % INVERT MAPPING VARTHETA3
    
    % [3, Sec. IV B) Step 7]
    ar_hat = cdf('normal',real(r),0,1/sqrt(2)); ar_hat = ar_hat(:);
    ai_hat = cdf('normal',imag(r),0,1/sqrt(2)); ai_hat = ai_hat(:); % first column
    br_hat = cdf('normal',real(s),0,1/sqrt(2)); br_hat = br_hat(:);
    bi_hat = cdf('normal',imag(s),0,1/sqrt(2)); bi_hat = bi_hat(:); % second column

    % HARD DECISIONS
    
    % [3, Sec. IV B) Step 8]
    A = repmat(ar_hat,1,length(lattice))-repmat(lattice,T-2,1);
    A = A.^2;
    [~,mininda] = min(A,[],2);

    B = repmat(ai_hat,1,length(lattice))-repmat(lattice,T-2,1);
    B = B.^2;
    [~,minindb] = min(B,[],2);

    C = repmat(br_hat,1,length(lattice))-repmat(lattice,T-2,1);
    C = C.^2;
    [~,minindc] = min(C,[],2);

    D = repmat(bi_hat,1,length(lattice))-repmat(lattice,T-2,1);
    D = D.^2;
    [~,minindd] = min(D,[],2);

    rx_symbols = zeros(2*M*(T-M),1);
    rx_symbols(1:2:M*(T-M)) = mininda;
    rx_symbols(2:2:M*(T-M)) = minindb;
    rx_symbols(M*(T-M)+1:2:2*M*(T-M)) = minindc;
    rx_symbols(M*(T-M)+2:2:2*M*(T-M)) = minindd;

else

    error('Number of transmit antennas not supported')

end

end