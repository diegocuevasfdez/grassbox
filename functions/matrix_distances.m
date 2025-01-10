function [DistanceMatrix, Indices, MinimumDistance, Ind_row, Ind_col] = matrix_distances(Codebook, distance)

%MATRIX_DISTANCES Distance matrix computation between all
%pairs of elements (subspaces) of a given codebook.
%
%   DISTANCEMATRIX = MATRIX_DISTANCES(CODEBOOK,DISTANCE)
%   computes the K x K specified distance matrix between all pairs of
%   elements (subspaces) of a given Grassmannian codebook of dimensions
%   T x M x K, where T is the coherence time, M is the number of antennas
%   and K is the number of codewords. Choices for DISTANCE are:
%
%           'chordal' - Chordal distance (default), d(X,Y) = 
%           ||X*X'-Y*Y'||_F/sqrt(2) = sqrt(M-sum(d_m^2)).
%
%           'procrustes' - Procrustes distance, d(X,Y) = sqrt(M-||X'*Y||_*)
%           *sqrt(2) =  sqrt(M-sum(d_m)), where ||X'*Y||_* is the nuclear 
%           norm (sum of cosines of principal angles).
%
%           'spectral' - Spectral distance, d(X,Y) = min_k (sin theta_k) =
%           sqrt(1-||X'*Y||_2^2).
%
%           'coherence' - Coherence criterium, d(X,Y) = sum(log(1 - d_m^2))
%           . Proposed in J. Alvarez-Vizoso, D. Cuevas, C. Beltran, 
%           I. Santamaria, V. Tucek and G. Peters, "Coherence-based 
%           Subspace Packings for MIMO Noncoherent Communications,"
%           30th European Signal Processing Conference (EUSIPCO), Belgrade,
%           Serbia, 2022, doi: 10.23919/EUSIPCO55093.2022.9909658.
%   
%   [DISTANCEMATRIX,INDICES] = MATRIX_DISTANCES(CODEBOOK, DISTANCE)
%   returns a vector containing the indices of the closest codeword to 
%   each of them.
%
%   [..., MINIMUMDISTANCE] = MATRIX_DISTANCES(CODEBOOK, DISTANCE)
%   returns the minimum pairwise distance among all codewords.
%
%   [..., IND_ROW, IND_COL] = MATRIX_DISTANCES(CODEBOOK, DISTANCE)
%   returns the indeces of the pair of codewords that produce the minimum
%   pairwise distance.
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

[~,M,K] = size(Codebook);

if nargin == 1
    distance = 'chordal';
end

DistanceMatrix = zeros(K); 
for i = 1:K
    for j = i+1:K
        if       strcmp(distance,'chordal')
            d = sqrt(M-norm(Codebook(:,:,i)'*Codebook(:,:,j),'fro')^2);

        elseif   strcmp(distance, 'procrustes')
            d = sqrt(M-sum(svd(Codebook(:,:,i)'*Codebook(:,:,j))));

        elseif   strcmp(distance, 'spectral')
            d = sqrt(1-norm(Codebook(:,:,i)'*Codebook(:,:,j)));

        elseif   strcmp(distance, 'coherence')
            d = svd(Codebook(:,:,i)'*Codebook(:,:,j));
            d = sum(log(1- d.^2));

        else  % TBC with other distances, so far we take chordal distance
            d = sqrt(M-norm(Codebook(:,:,i)'*Codebook(:,:,j),'fro')^2);

        end

        DistanceMatrix(i,j) = d;
        DistanceMatrix(j,i) = d;
    end
end

DistanceMatrix = real(DistanceMatrix) + 1e10*eye(K);

[~,Indices] = min(DistanceMatrix,[],2); % indexes for minimum distances 
                                        % (for each codeword)

[MinimumDistance,I] = min(DistanceMatrix(:)); % total minimum distance 

[Ind_row, Ind_col] = ind2sub(size(DistanceMatrix),I); % (Ind_row,Ind_col) 
                                                      % codeword is the 
                                                      % pair producing the 
                                                      % minimum distance




