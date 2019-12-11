function EkofiskInitialEnsemble_perm3D(makeTrue,ensembleSize, randomNumber)

% This function generate an initial ensemble for the full ekofisk field.
% The ensemble is generated for porosity, log-permeability, net-to-gross,
% z-multipliers for base of layers [1,8,11,12,15,18], fault multipliers,
% oil-water contacts, multipliers for scaling of rel-perm curves, and
% region multipliers. Total number of generated parameters is 148159.
%
% For more information we refer to the paper: 
%
% Lorentzen, R., Luo, X., Bhakta, T., Valestrand, R.: "History matching
% ekofisk reservoir and petroelastic models using seismic impedance with
% correlated noise" Submitted to SPE Journal.
% 
% We also ask that the above paper is cited in publications aided by this
% code.
%
% Input:
% - ensembleSize : Number of generated ensemble members (defalult 100)
% - randomNumber : Number used to initialize random generator 
%                  (default 1.2345e5)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C): IRIS (International Research Institute of Stavanger), 2017. 
% Contact: Rolf.Lorentzen@iris.no


% set random generator
if nargin < 3
    rng(1.2345e5);
else
    if isscalar(randomNumber) && isreal(randomNumber) && ...
       rem(randomNumber,1) == 0 && randomNumber >= 0
        rng(randomNumber);
    else
        error('Random seed must be a nonnegative integer')
    end
end

% number of ensemble members
if nargin < 2
    N = 100;
else
    if isscalar(ensembleSize) && isreal(ensembleSize) && ...
       rem(ensembleSize,1) == 0 && ensembleSize > 0
        N = ensembleSize;
    else
        error('Ensemble size must be a positive integer')
    end
end

if nargin > 0 && makeTrue==1
    N = N + 1;
end

% set geostatistical parameters and initialize
ekofisk = EkofiskGeostat_perm3D();
A = ekofisk.actnum; % mask for active gridcells
N_A = sum(A); % number of active gridcells
D = ekofisk.dim; % reservoir dimention
N_L = D(1)*D(2); % number of gridcells in a layer
N_F = prod(D); % total number of gridcells
N_S = 3; % number of static fields
ensemble = zeros(N_A*N_S,N); % initialize ensemble

% mean values for poro, log-permx and ntg
M(:,1) = ekofisk.permxLogMean;
M(:,2) = ekofisk.permxLogMean;
M(:,3)= ekofisk.permxLogMean;

% std for permx, permy and permz
S(:,1) = ekofisk.permxStd;
S(:,2) = ekofisk.permyStd;
S(:,3) = ekofisk.permzStd;

% mean correlation lengths (ranges)
C(:,1) = ekofisk.permxRange;
C(:,2) = ekofisk.permxRange;
C(:,3) = ekofisk.permxRange;

% std for correlation lengths
if isfield('ekofisk','C_S')
    C_S = ekofisk.C_S;
else
    C_S = 2; 
end

% Poro / Permx correlation
R1 = ekofisk.PermxPermyCorr;
R2 = ekofisk.PermxPermzCorr;


% initial ensemble for Poro, Permx,
for I = 1:N
    
    % Permx
    B = 1;
    E = N_A;
    X1 = GaussianWithVariableParameters(D,zeros(N_F,1),1,C(:,1),C_S);
    ensemble(B:E,I) = M(:,1) + S(:,1).*X1(A==1);
    
    % Permy
    B = B + N_A;
    E = E + N_A;
    X2 = GaussianWithVariableParameters(D,zeros(N_F,1),1,C(:,2),C_S);
    X = R1*X1 + sqrt(1-R1^2)*X2;
    ensemble(B:E,I) = M(:,2) + S(:,2).*X(A==1);
    
    % Permz
    B = B + N_A;
    E = E + N_A;
    X3 = GaussianWithVariableParameters(D,zeros(N_F,1),1,C(:,3),C_S);
    X = R2*X1 + sqrt(1-R2^2)*X3;
    ensemble(B:E,I) = M(:,3) + S(:,3).*X(A==1);
   
end

% Define bounds
G = ones(N_A,1);
S_LB = [ekofisk.permxLB*G; ekofisk.permxLB*G; ekofisk.permxLB*G];
S_UB = [ekofisk.permxUB*G;ekofisk.permxUB*G; ekofisk.permxUB*G];

% Adjust parameters with bounds
ensemble = AdjustVariableWithInBounds(ensemble,S_LB,S_UB); %#ok<*NASGU>

% condition on well logs
%%% conditioned_reducedLayer = 1, otherwise not
%conditionInitialEnsembleEkofisk(ensemble);

% Save the ensemble
 save('ekofiskInitEnsPerm3D.mat', 'ensemble');


%-----------------------------------------------
% Subfunction 
%----------------------------------------------- 
function [x,corrLength] = GaussianWithVariableParameters(fieldDim,meanValue, ...
			  Sdev,meanCorrLength,stdCorrLength)

% Setup a gaussian field with correlation length drawn from a
% normal distribution. The horizontal layers are genereated
% independently
%
% Input:
% - fieldDim       : dimension of the field.
% - meanValue      : the mean value of the field (vector of the size
%                    of the field).
% - Sdev           : standard deviation of the field.
% - meanCorrLength : mean correlation length.
% - stdCorrLength  : standard deviation of the correlation length.
%
% Output:
% - x              : the generated field.
% - corrLenght     : the drawn correlation length.

corrLength=AddGnoise(meanCorrLength,stdCorrLength,1);
if length(fieldDim)<3
    
    x=meanValue+FastGaussian(fieldDim,Sdev,corrLength);
    
else
    
    layerDim=prod(fieldDim(1:2));
    
    % Initialization
    x=meanValue;
    if length(Sdev)==1
        for I=1:fieldDim(3)
            x(1+(I-1)*layerDim:I*layerDim,1)=meanValue(1+(I-1)*layerDim: ...
                I*layerDim)+FastGaussian(fieldDim(1:2),Sdev,corrLength);
            
            % Generate new correlation length for the next layer
            corrLength=AddGnoise(meanCorrLength,stdCorrLength,1);
        end
    else
        for I=1:fieldDim(3)
            x(1+(I-1)*layerDim:I*layerDim,1)=meanValue(1+(I-1)*layerDim: ...
                I*layerDim)+ ...
                FastGaussian(fieldDim(1:2),Sdev(1+(I-1)*layerDim:I* ...
                layerDim),corrLength);
            
            % Generate new correlation length for the next layer
            corrLength=AddGnoise(meanCorrLength,stdCorrLength,1);
        end
    end
    
end


%-----------------------------------------------
% Subfunction 
%----------------------------------------------- 
function [variable,n] = AdjustVariableWithInBounds(variable, ...
						lowerbound,upperbound)

% The function returns variable such that 
%    lowerbound <= variable(i) <= upperbound,
% i.e. if variable(i)<lowerbound before calling the function, then
% variable(i)=lowerbound as output. Similarly for upperbound, both
% with inequality the opposite way. If there is no
% lowerbound/upperbound, set the corresponding matrix empty.
%
% Input:
% - variable   : Vector or matrix. Variables (or an ensemble of variable 
%                samples) to be checked and truncated (if necessary) 
% lowerbound   : Scalar or vector. Lower bound(s) for the variables to be
%                checked
% upperbound   : Scalar or vector. Upper bound(s) for the variables to be
%                checked
% Output:
% - variable   : Variables after check/truncation
% - n          : Number of truncations

if nargin < 3
  error(['adjustVariableWithInBounds is only implemented for three' ...
	 ' arguments'])
end

ne=size(variable,2);
n=0;
if ~isempty(lowerbound)
    if max(size(lowerbound))==1
        n=sum(variable<lowerbound);
        variable(variable<lowerbound)=lowerbound;
    else
        lowerbound=kron(lowerbound,ones(1,size(variable,2)));
        n=sum(variable<lowerbound);
        variable(variable<lowerbound)=lowerbound(variable<lowerbound);
    end
end

if ~isempty(upperbound)
    if max(size(upperbound))==1
        n=n+sum(variable>upperbound);
        variable(variable>upperbound)=upperbound;
    else
        upperbound=kron(upperbound,ones(1,ne));
        n=n+sum(variable>upperbound);
        variable(variable>upperbound)=upperbound(variable>upperbound);
    end
end


%-----------------------------------------------
% Subfunction 
%----------------------------------------------- 
function [Y,RTSIGMA] = AddGnoise( Ytrue,SIGMA,SQ )

% Add noise, normally distributed, with covariance given by SIGMA.
%
% Input:
% - Ytrue    Original signal.
% - SIGMA    Specification of covariance matrix of noise.  May be
%            entered as scalar, vector or full matrix. If it SIGMA
%            is a vector, then it is interpreted as the covariance
%            matrix is diag(SIGMA).
% - SQ       If present, determine whether SIGMA or SIGMA*SIGMA' is used
%            as the covariance matrix.  Thus, if the square root of
%            the covariance matrix has allready been calculated
%            previously, work may be saved by setting SQ to 1.
%
% Output:
% - Y        Signal with noise added.
% - RTSIGMA  The square root of SIGMA; RTSIGMA*RTSIGMA' = SIGMA.
%            (Helpful if it is cumbersome to compute).

% Compute the normally distributed noise, with covariance matrix
% SIGMA or SIGMA*SIGMA'.
try
    
    if nargin > 2 && SQ == 1
        % Use SIGMA*SIGMA' as covariance matrix
        RTSIGMA = SIGMA ;
        if min(size(SIGMA)) == 1
            % SIGMA is a scalar or vector
            error = RTSIGMA.*randn(size(Ytrue)) ;
        else
            error = RTSIGMA*randn(size(RTSIGMA,2),1) ;
        end
    else
        % Use SIGMA as covariance matrix
        if min(size(SIGMA))==1
            % SIGMA is entered as a scalar or a vector
            RTSIGMA = realsqrt(SIGMA);
            error = RTSIGMA.*randn(size(Ytrue));
        else
            [RTSIGMA,p] = chol(SIGMA); % The matrix must be transposed.
            if p>0
                disp('Problem with Cholesky factorization')
                disp(['p = ',num2str(p)]);
                RTSIGMA = real(sqrtm(SIGMA));
                disp('Finnaly - we got a square root!')
            end
            RTSIGMA = RTSIGMA';
            error = RTSIGMA*randn(size(Ytrue));
        end %if
    end
    % Add the noise:
    Y = Ytrue+error;
  
catch err
    
  disp('Error in AddGnoise');
  disp('Size Ytrue:')
  size(Ytrue)
  disp('Size SIGMA:')
  size(SIGMA)
  rethrow(err);
  
end

