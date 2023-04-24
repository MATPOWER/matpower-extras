function  [CmpltObsvrblIdx, ObsvrblIdx] = selectObsrvblStates(casename, Strtgy)
% Select states in the power system such that the system is observable.
% Returns the indices of such states.
% % Input(s):
% * casename: (string) name of the test case.
% * Strtgy  : (string) name of the strategy to select the observable states.
% Output(s):
% * CmpltObsvrblIdx: (struct) complete set of indices including the
% observable states indices in the test case.
% * ObsvrblIdx    : (struct) observable state indices of the selected states 
% per the test case file.

switch lower(Strtgy)
     
    case 'critical'
        switch casename
            case 'case14'
                % Measurements from Davoudi Thesis for 14 bus case
                ObsvrblIdx.idx_zVm = [1];
                ObsvrblIdx.idx_zPG = [1,2,4,5];
                ObsvrblIdx.idx_zQG = [1,2,4,5];
                ObsvrblIdx.idx_zPD = [9,10,11,12,13,14];
                ObsvrblIdx.idx_zQD = [9,10,11,12,13,14];
                ObsvrblIdx.idx_zPF = [6,10,14];
                ObsvrblIdx.idx_zQF = [6,10,14];
                
            case 'case_ieee30'
                % Buses where selected as follows:
                % 1- All generator buses except ref.
                % 2- All demad/load buses.
                % 3- Power flows to make the above buses set. i.e., to
                % cover all the buses in the system.
                % 4- Power flow 9-10 (index 14) was added since the
                % flow and connection 9-11 looks redundant.
                % 5- Vm at 30. Bus 30 had a zero pivot in the QR
                % decomposition. As such all zero pivots where included
                % one-by-one and it turned out that it is sufficient to
                % include bus 30 alone!
                ObsvrblIdx.idx_zVm = [30];
                ObsvrblIdx.idx_zPG = [2,3,4,5,6];
                ObsvrblIdx.idx_zQG = [2,3,4,5,6];
                ObsvrblIdx.idx_zPD = [3,4,7,10,12,14,15,16,17,18,19,20,21,23,24,26,29,30];
                ObsvrblIdx.idx_zQD = [3,4,7,10,12,14,15,16,17,18,19,20,21,23,24,26,29,30];
                ObsvrblIdx.idx_zPF = [9,13,14,31,34,37,36];
                ObsvrblIdx.idx_zQF = [9,13,14,31,34,37,36];
        end

    case 'davoudi'
        switch casename
            case 'case14'
                % Measurements from Davoudi Thesis for 14 bus case
                ObsvrblIdx.idx_zVm = [1,2,3,4,8];
                ObsvrblIdx.idx_zPG = [1,2,4,5];
                ObsvrblIdx.idx_zQG = [1,2,4,5];
                ObsvrblIdx.idx_zPD = [9,10,11,12,14];
                ObsvrblIdx.idx_zQD = [9,10,11,12,14];
                ObsvrblIdx.idx_zPF = [1,2,3,4,6,7,8,10,13];
                ObsvrblIdx.idx_zQF = [1,2,3,4,6,7,8,10,13];
                
            case 'case_ieee30'
                error('Not implemented!')

        end
    case 'bo'
        % Measurements from original test_se_14bus.m file by Bo.
        ObsvrblIdx.idx_zPD = [];
        ObsvrblIdx.idx_zQD = [];
        ObsvrblIdx.idx_zPF = [1;3;8;9;10;13;15;16;17;19];
        ObsvrblIdx.idx_zPT = [4;5;7;11];
        ObsvrblIdx.idx_zPG = [1;2;3;4;5];
        ObsvrblIdx.idx_zVa = [];
        ObsvrblIdx.idx_zQF = [1;3;8;9;10;13;15;19];
        ObsvrblIdx.idx_zQT = [4;5;7;11];
        ObsvrblIdx.idx_zQG = [1;2];
        ObsvrblIdx.idx_zVm = [2;3;6;8;10;14];
        
    case 'kerdchuen'
    switch casename
        case 'case14'
            % Measurements from Kerdchuen2009 paper with modifications. Bus 7
            % was removed and Vm at Bus 1 was added.
            ObsvrblIdx.idx_zPG = [4];
            ObsvrblIdx.idx_zQG = [4];
            ObsvrblIdx.idx_zPD = [4, 5, 13];
            ObsvrblIdx.idx_zQD = [4, 5, 13];
            ObsvrblIdx.idx_zPF = [1,6,11,14,15,16,18,19,20];
            ObsvrblIdx.idx_zQF = [1,6,11,14,15,16,18,19,20];
            ObsvrblIdx.idx_zPT = [];
            ObsvrblIdx.idx_zVa = [];
            ObsvrblIdx.idx_zQT = [];
            ObsvrblIdx.idx_zVm = [1];
        
        case 'case30'
            % Measurements from Shahraeini2011 paper
            % Not Working
            ObsvrblIdx.idx_zPG = [2, 11, 27];
            ObsvrblIdx.idx_zQG = [2, 11, 27];
            ObsvrblIdx.idx_zPD = [6, 10, 12, 25, 3, 4, 9, 15, 20, 28, 29];
            ObsvrblIdx.idx_zQD = [6, 10, 12, 25, 3, 4, 9, 15, 20, 28, 29];
            ObsvrblIdx.idx_zPF = [2,3,8,13,14,16,19,20,23,26,27,29,30,31,34,38,40];
            ObsvrblIdx.idx_zQF = [2,3,8,13,14,16,19,20,23,26,27,29,30,31,34,38,40];
            ObsvrblIdx.idx_zPT = [];
            ObsvrblIdx.idx_zVa = [];
            ObsvrblIdx.idx_zQT = [];
            ObsvrblIdx.idx_zVm = [];
        
    end
             
    otherwise
        error('Not implemented!')
end

% Create a complete index structure 
CmpltObsvrblIdx = cmpltInd(ObsvrblIdx);


%% Supporting function(s)
function CmpltIdx = cmpltInd(idx)
% Fill the missing entries of the index structure with empty array ([]).

% All field names in the index structures 
idx_fields = {'idx_zPF', 'idx_zPT','idx_zQF', 'idx_zQT', ...
    'idx_zPG', 'idx_zQG', 'idx_zPD', 'idx_zQD', 'idx_zVm','idx_zVa'};

% Available field names in the given index structure
fn_idx = fieldnames(idx);
CmpltIdx = idx;

for i = 1:numel(idx_fields)
    % Check if required field name is given
   if ~any(strcmp(idx_fields{i}, fn_idx))
      % Add empty field name
      CmpltIdx.(idx_fields{i}) = [];
   end
end

