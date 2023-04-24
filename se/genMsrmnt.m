function [noisyMeas, CompleteObsIdx, sigma_mod] = genMsrmnt(casename, Strtgy, P, Q, V,...
    MpcInd, baseMVA, sigma, varargin)
% Generate measurements from the power data extracted from the test
% case.
%
% Inputs:
% ------------------------------------------------------------------------
% casename: (string) test case file name.
% Strtgy  : (string) strategy of choosing observable states.
% P       : (struct) containing real power values. 
%   .PG   : (array) real generated power.
%   .PD   : (array) real power load demand.
%   .PS   : (array) real power shunt consumption.
%   .PF   : (array) real power from-to bus.
%   .PT   : (array) real power to-from bus.
%
% Q       : (struct) containing reactive power values.
%   .QG   : (array) reactive generated power.
%   .QD   : (array) reactive power load demand.
%   .QS   : (array) reactive power shunt consumption.
%   .QF   : (array) reactive power from-to bus.
%   .QT   : (array) reactive power to-from bus.
%
% V       : (struct) bus voltages.
%   .Vm   : (array) voltage magnitudes.
%   .Va   : (array) voltage angles.
%
% MpcInd  : (struct) containing indicies of the buses for above values.
%   .PG   : (array) real generated power index.
%   .PD   : (array) real power load demand index.
%   .PS   : (array) real power shunt consumption index.
%   .PF   : (array) real power from-to bus index.
%   .PT   : (array) real power to-from bus index.
%   .QG   : (array) reactive generated power index.
%   .QD   : (array) reactive power load demand index.
%   .QS   : (array) reactive power shunt consumption index.
%   .QF   : (array) reactive power from-to bus index.
%   .QT   : (array) reactive power to-from bus index.
%
% NOTE: The all the values are ordered as per the mpc structure.
%
% baseMVA : (scalar) power base value.
%
% sigma   : (struct) standard deviations of the specified measurements.
% Unused measurements are left as empty arrays.
%
% Outputs:
% ------------------------------------------------------------------------
% noisyMeas: (struct) noisy measurements of the specified measurements. 
% Unused measurements are left as empty arrays.
%
% Idx      :  (struct) indices of the specified measurements. Unused 
% measurements are left as empty arrays.
%
% sigma_mod: (struct) standard deviation formated to suit the doSE.m
% function format. 


% Return the observable states indices
CompleteObsIdx = selectObsrvblStates(casename, Strtgy, MpcInd);

% Extract specific data
meas = extrctMeas(P, Q, V, baseMVA, MpcInd, CompleteObsIdx);

% Add noise 
if isempty(varargin)
    [noisyMeas, sigma_mod] = addNoise(meas, sigma);
else
    [noisyMeas, sigma_mod] = addNoise(meas, sigma, varargin);
end

%% Supporting functions
function meas = extrctMeas(P, Q, V, baseMVA, CmpltInd, ObsvrIdx)
% Extract measured data referenced by the given indicies.
%
% Inputs:
% ------------------------------------------------------------------------
% P, Q, V, baseMVA: as defined in genMsrmnt.
% CmpltInd: (struct) complete indices of values in the test case.
% ObsvrIdx: (struct) indices of the observable states in the test case.
% Other indices are left empty.


% Available field names in the given index structure
fn_idx = fieldnames(ObsvrIdx);

for i = 1:numel(fn_idx)
    
    % Measurement type
    meas_type = fn_idx{i}(end-1:end);
    
    % Convert indicies to appropriate reference
    if strcmp(meas_type(2),'G')
        meas_ind = ObsvrIdx.(fn_idx{i});
    else
        meas_ind = ismember(CmpltInd.(meas_type), ObsvrIdx.(fn_idx{i}));
    end
    
   
    switch meas_type(1)
        
        % Real power
        case 'P'
            meas.(meas_type) = P.(meas_type)(meas_ind)/baseMVA;
        % Reactive power
        case 'Q'
            meas.(meas_type) = Q.(meas_type)(meas_ind)/baseMVA;
        % Voltage
        case 'V'
            meas.(meas_type) = V.(meas_type)(meas_ind);
                  
    end
    
end

function [noisyMeas, sigma_mod] = addNoise(meas, sigma, varargin)
% Add noise to data specified by the provided standard deviation.
%
% Inputs:
% ------------------------------------------------------------------------
% meas : (struct) containing the measurements from the power system defined
% in the test case.
% sigma: (struct) containing the standard devaition of the corresponding
% measurements.
% Outputs:
% ------------------------------------------------------------------------
% noisyMeas: (struct) containing noisy measurements from the power system 
% defined in the test case.    
% sigma_mod: (struct) standard deviation formated to suit the doSE.m
% function format. 

% Account for variable number of input arguments
if nargin == 4
    nMeas = cell2mat(varargin{1});
    VarType = varargin{2};
elseif nargin == 3
    nMeas = cell2mat(varargin{1});
    VarType = 'equal';
elseif nargin == 2
    nMeas = 1;
    VarType = 'equal';
else
    error('Unsupported number of arguments!')    
end

% Get field names
fn_sigma = fieldnames(sigma);

% Create field names for "sigma" structure
ErasFun = @(x) erase(x, 'sigma_');
fn = cellfun(ErasFun, fn_sigma ,'UniformOutput',false);

% Extension to matrix 
for i = 1:numel(fn_sigma)
   if isempty(sigma.(fn_sigma{i}))
       noisyMeas.(fn{i}) = meas.(fn{i});
   else
       % Extend sigma
       if strcmpi( (VarType), 'equal' )
           sigma_mod.(fn_sigma{i}) = sigma.(fn_sigma{i}) * ones(numel(meas.(fn{i})), 1);
       elseif strcmpi( VarType, 'distinct' )
           sigma_mod.(fn_sigma{i}) = sigma.(fn_sigma{i});
       end

       % Add noise
       noisyMeas.(fn{i}) = meas.(fn{i}) + sigma_mod.(fn_sigma{i}).*randn(numel(meas.(fn{i})),nMeas);
   end
end


