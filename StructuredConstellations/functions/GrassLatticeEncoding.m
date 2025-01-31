function X = GrassLatticeEncoding(M,lattice,tx_symbols)

%GRASSLATTICEENCODING Encoding for the Grass-Lattice constellation. This
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
%   X = GRASSLATTICEENCODING(M, LATTICE, TX_SYMBOLS) generates the 
%   Grass-Lattice encoded codeword on a given "lattice" and for the desired 
%   "tx_symbols".
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
%   has P points, the code uses B=log2(P) bits per coordinate. We encode a
%   total of 2M(T-M) coordinates, where T is the coherence time.
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

T = length(tx_symbols)/(2*M) + M; % coherence time

values = repmat(lattice, 2*M*(T-M),1); % generate all lattices
linidx = sub2ind(size(values), 1:size(values,1), tx_symbols');  
xtilde = values(linidx);xtilde=xtilde(:); % take values corresponding to transmitted symbols on each lattice

if M == 1

    a = xtilde(1:2:end); a = a(:);
    b = xtilde(2:2:end); b = b(:);

    % MAPPING VARTHETA1 
    % Inverse transform sampling to generate N(0,1)
    
    % [2, Sec. IV A) Step 1]
    an = norminv(a,0,1/sqrt(2));
    bn = norminv(b,0,1/sqrt(2));
    z = an + 1i*bn;

    % MAPPING VARTHETA2
    % Compute w= zf_{T-1}(||z||)
    
    % [2, Sec. IV A) Step 2]
    normz = norm(z);
    fz = cdf('chi2',2*normz^2, 2*(T-1))^(1/(2*(T-1)))/normz;  % an alternative way to calculate f(z) proposed in [2,Remark 1]
    w = z*fz;

    % MAPPING VARTHETA3
    
    % [2, Sec. IV A) Step 3]
    X = [sqrt(1-w'*w); w]; % encoded vector in Gr(1,C^T) [2,Lemma 4]

elseif M == 2

    ar = xtilde(1:2:M*(T-M)); ar = ar(:);
    ai = xtilde(2:2:M*(T-M)); ai = ai(:);
    br = xtilde(M*(T-M)+1:2:2*M*(T-M)); br = br(:);
    bi = xtilde(M*(T-M)+2:2:2*M*(T-M)); bi = bi(:);

    % MAPPING VARTHETA1
    % Inverse transform sampling to generate N(0,1)
    
    % [3, Sec. IV A) Step 1]
    arn = norminv(ar,0,1/sqrt(2));
    ain = norminv(ai,0,1/sqrt(2)); % first column
    brn = norminv(br,0,1/sqrt(2));
    bin = norminv(bi,0,1/sqrt(2)); % second column

    r = arn + 1i*ain; 
    s = brn + 1i*bin; % matrix Z = [r s] with CN(0,1) entries

    % MAPPING VARTHETA2
    
    % [3, Sec. IV A) Step 2]
    nr = norm(r);
    u = r'*s/nr^2*r;
    sbar = s-u;
    
    % [3, Sec. IV A) Step 3]
    p = r/nr*thnt(T-2,nr);
    nsbar = norm(sbar);
    if norm(sbar)==0
        qbar = 0*sbar;
    else
        qbar = sbar/nsbar*thnt(T-3,nsbar);
    end
    
    % [3, Sec. IV A) Step 4]
    if norm(u)==0
        v = 0*u;
    else
        f1_normu = 1/norm(u) * (1-exp(-norm(u)^2))^(1/2); % [3, Eq.(8)]
        v = u * f1_normu * sqrt((1-norm(p)^2)*(1-norm(qbar)^2)); % [3, Eq.(7)]
    end

    % [3, Sec. IV A) Step 5]
    W = [p qbar+v]; % this matrix is uniformly distributed in the set of matrices of operator norm at most 1

    % MAPPING VARTHETA3
    
    % [3, Sec. IV A) Step 6]
    X = [sqrtm(eye(2)-W'*W); W]; % encoded vector in Gr(2,C^T) proposed in [2, Corollary 1]

else

    error('Number of transmit antennas not supported');

end

end