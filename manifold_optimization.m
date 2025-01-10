function [C,varargout] = manifold_optimization(C_initial,varargin)

%MANIFOLD_OPTIMIZATION Grassmannian codebook design refinement via
%optimization on the manifold (GrassManOpt). This optimization algorithm is
%proposed in D. Cuevas, C. Beltran, I. Santamaria, V. Tucek and G. Peters,
%"A Fast Algorithm for Designing Grassmannian Constellations," 25th
%International ITG Workshop on Smart Antennas, French Riviera, France,
%2021.
%
%   C = MANIFOLD_OPTIMIZATION(C_INITIAL)
%   Constellation design that maximizes the 
%   minimum pairwise distance (or criterion) between codewords.
%
%   [C,SEPARATION] = MANIFOLD_OPTIMIZATION(C_INITIAL) returns the
%   minimum pairwise distance among all codewords.
%
%   [ ... ] = MANIFOLD_OPTIMIZATION(..., params) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by MANIFOLD_OPTIMIZATION. Parameters are:
%
%   'Alpha' - Parameter that controls the precision of the
%   optimization algorithm. Default is 1.001 (must be >1).
%
%   'MaxIter' - Maximum number of iterations allowed for
%   optimization refinement step. Default is 200.
%
%   'hini' - Initial step-size for the optimization algorithm.
%   Default is 0.1.
%
%   'hfin' - Minimum step-size for the optimization algorithm.
%   Default is 1e-4.
%
%   'MinImp' - Minimum improvement in terms of distance between codewords
%   that should be achieved in every iteration for the optimization 
%   algorithm to continue. Default is 1e-5.
%
%   'Distance' - Distance used in the optimization algorithm. Choices are:
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
%           'coherence' - Coherence criterion, d(X,Y) = sum(log(1 - d_m^2))
%           . Proposed in J. Alvarez-Vizoso, D. Cuevas, C. Beltran, 
%           I. Santamaria, V. Tucek and G. Peters, "Coherence-based 
%           Subspace Packings for MIMO Noncoherent Communications,"
%           30th European Signal Processing Conference (EUSIPCO), Belgrade,
%           Serbia, 2022, doi: 10.23919/EUSIPCO55093.2022.9909658.
%
% This file is part of the Grassmannian Constellations Toolbox for Matlab.
% https://github.com/diegocuevasfdez/grassbox/

% Default values
opt_params = struct();
opt_params.Alpha = 1.001;
opt_params.MaxIter = 200;
opt_params.hini = 0.1;
opt_params.hfin = 1e-4;
opt_params.MinImp = 1e-5;
opt_params.Distance = 'chordal';
opt_params.verbose = false;


if nargin < 1
    error(message('TooFewInputs'));
elseif nargin == 2
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'Alpha'
                opt_params.Alpha = param_value;
            case 'MaxIter'
                opt_params.MaxIter = param_value;
            case 'hini'
                opt_params.hini = param_value;
            case 'hfin'
                opt_params.hfin = param_value;
            case 'MinImp'
                opt_params.MinImp = param_value;
            case 'Distance'
                opt_params.Distance = param_value;
            case 'verbose'
                opt_params.verbose = param_value;
        end
    end
elseif nargin > 2
    error(message('TooManyInputs'));
end
 
% Get codebook dimensions
T = size(C_initial,1);
M = size(C_initial,2);
num_codewords = size(C_initial,3);

% Initial codebook
C = C_initial;

% Compute initial minimum pairwise distance
[~, indices, separation,~,~] = matrix_distances(C, opt_params.Distance);
fprintf('Initial separation=%f\n',separation);

% GrassManOpt algorithm
h = opt_params.hini; % initial step-size

Matrixaux = zeros(M);
Matrixaux(1,1) = 1;
break_flag = 0;
Cnew = C;

for k = 1:opt_params.MaxIter

    improvement = 0;

    while (~improvement) && h > opt_params.hfin

        h = h/opt_params.Alpha; % decrease step-size

        for i = 1:num_codewords

            j = indices(i); % index of closest codeword to i-th codeword
            Qi = C(:,:,i); % i-th codeword
            Qj = C(:,:,j); % j-th codeword

            % Gradient computation
            if      strcmp(opt_params.Distance,'chordal')
                direction_i = -(eye(T)-Qi*Qi')*(Qj*(Qj'*Qi));

            elseif  strcmp(opt_params.Distance,'procrustes')
                direction_i = -(eye(T)-Qi*Qi')*Qj*(sqrtm(Qj'*(Qi*Qi')*Qj)^(-1))*(Qj'*Qi);

            elseif strcmp(opt_params.Distance,'spectral')
                [Uaux,~,Vaux] = svd(Qj'*Qi);
                directionaux = Uaux*Matrixaux*Vaux';
                direction_i = -(eye(T)-Qi*Qi')*(Qj*(directionaux));

            elseif strcmp(opt_params.Distance,'coherence')
                aux = Qj'*Qi;
                directionaux = aux/(eye(M)-aux'*aux);
                direction_i = -(eye(T)-Qi*Qi')*(Qj*(directionaux));

            else % TBC with other distances, so far we take the chordal distance
                direction_i = -(eye(T)-Qi*Qi')*(Qj*(Qj'*Qi));
            end

            direction_i = direction_i/norm(direction_i,'fro');
            Qi = Qi+h*direction_i; % separate i-th codeword in the direction given by the gradient
            [Qi,~] = qr(Qi); % retract i-th codeword to the manifold
            Cnew(:,:,i) = Qi(:,1:M);

        end

        [~, indicesnew, separation2, ~, ~] = matrix_distances(Cnew, opt_params.Distance); % recompute minimum pairwise distance

        if separation2 > separation

            C = Cnew; % if there is improvement --> update codebook

            if (separation2-separation) < opt_params.MinImp
               break_flag = 1; % stop optimization
            end

            improvement = 1;
            separation = separation2;
            h = opt_params.Alpha*h; % increase step-size
            indices = indicesnew;

        else
            improvement = 0;          
        end
    end

    if opt_params.verbose == true
    fprintf('Progress percentage %d, h= %f, Actual separation=%f\n',...
                floor(100*k/opt_params.MaxIter),h,separation);
    end
    
    if break_flag % if improvement is below MinImp --> end optimization
       break; 
    end
end

fprintf('Final separation=%f\n',separation);

opt_params.h = h;
varargout{1} = separation;
varargout{2} = opt_params;  % final parameters (after optimization)

end