function [mpm, success] = OLTC(mpm,idx_branches_in_oltc)
%% OLTC model for changing the tap-ratios

%%
% Roger Germann, Kerry Jansen, Denis Mikhaylov 02-05.2020
% Updated by Lukas Ortmann on 23-11.2020

%% Copyright (c) 2023 ETH Zurich, Automatic Control Laboratory, Lukas Ortmann
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its
% contributors may be used to endorse or promote products derived from
% this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

define_constants;
%% Check that OLTCs are enabled
if mpm.oltc.enable
    
    number_branches_oltc = size(idx_branches_in_oltc,1);
    
    %% Set the OLTC parameters
    V_deadband = mpm.oltc.V_deadband;
    TAP_stepsize = mpm.oltc.TAP_stepsize;
    TAP_lbound = 0.9 + TAP_stepsize;
    TAP_ubound = 1.1 - TAP_stepsize;
    %% Set the voltage setpoints
    V_setpoint = ones(number_branches_oltc,1);
    
    
    %% Get the indices of the buses the OLTCs are connected to on
    %% 'from and 'to' side
    
    
    % Initialise an empty array to save the from and to bus ids where
    % the OLTCs are positioned
    idx_bus_oltc_f = zeros(number_branches_oltc,1);
    idx_bus_oltc_t = zeros(number_branches_oltc,1);
    % Go through all branches
    for i = 1 : number_branches_oltc
        idx_bus_oltc_f(i) = find(mpm.mpc.bus(:,1) == mpm.mpc.branch(idx_branches_in_oltc(i),F_BUS));
        idx_bus_oltc_t(i) = find(mpm.mpc.bus(:,1) == mpm.mpc.branch(idx_branches_in_oltc(i),T_BUS));
    end
    
    %% Determine bus indices of lower voltage side of the OLTCs
    BASEf = mpm.mpc.bus(idx_bus_oltc_f,BASE_KV);
    BASEt = mpm.mpc.bus(idx_bus_oltc_t,BASE_KV);
    
    % BASE_KV on from side lower
    idx_bus_f_lower = idx_bus_oltc_f(BASEf<BASEt);
    idx_branches_f_lower = idx_branches_in_oltc(BASEf<BASEt);
    % BASE_KV on to side lower
    idx_bus_t_lower = idx_bus_oltc_t(BASEf>BASEt);
    idx_branches_t_lower = idx_branches_in_oltc(BASEf>BASEt);
    
    if mpm.oltc.verbose == true
        TAP_position = zeros(number_branches_oltc,1);
        VMf0 = mpm.mpc.bus(idx_bus_oltc_f,VM);
        VMt0 = mpm.mpc.bus(idx_bus_oltc_t,VM);
    end
    
    %% Change the tap ratios as long as voltage violations occur
    %% or until the maximum number of iterations is reached
    i = 0;
    TAPs_changed = 1;
    while TAPs_changed && i < mpm.oltc.max_iter
        i = i + 1;
        TAPs_changed = 0;
        % Solve the power flow
        [mpm.mpc, success] = runpf(mpm.mpc,mpm.mpopt);
        % If runpf does not converge, return directly
        if success == false
            fprintf("runpf did not converge on %d iteration. \n", i);
            return
        end
        if mpm.oltc.verbose == true
            TAP_before = mpm.mpc.branch(idx_branches_in_oltc, TAP);
        end
        
        % Get voltage magnitude on from and to side in p.u.
        VMf = mpm.mpc.bus(idx_bus_f_lower,VM);
        VMt = mpm.mpc.bus(idx_bus_t_lower,VM);
        
        % Voltage on FROM side lower, OLTC on FROM side
        % Check voltage violations and change tap ratios
        idx_decr = idx_branches_f_lower(VMf > V_setpoint + V_deadband);
        idx_incr = idx_branches_f_lower(VMf < V_setpoint - V_deadband);
        % Decrease TAP ratios, check lower bound
        TAP_decr = mpm.mpc.branch(idx_decr, TAP) > TAP_lbound;
        mpm.mpc.branch(idx_decr, TAP) = mpm.mpc.branch(idx_decr, TAP) - TAP_decr * TAP_stepsize;
        % Increase TAP ratios, check upper bound
        TAP_incr = mpm.mpc.branch(idx_incr, TAP) < TAP_ubound;
        mpm.mpc.branch(idx_incr, TAP) = mpm.mpc.branch(idx_incr, TAP) + TAP_incr * TAP_stepsize;
        
        % Check if any TAP ratios were changed
        TAPs_changed = TAPs_changed | any(TAP_incr) | any(TAP_decr);
        
        % Voltage on TO side lower, OLTC on TO side
        % Check voltage violations and change tap ratios
        idx_decr = idx_branches_t_lower(VMt < 1 - V_deadband);
        idx_incr = idx_branches_t_lower(VMt > 1 + V_deadband);
        % Decrease TAP ratios, check lower bound
        TAP_decr = mpm.mpc.branch(idx_decr, TAP) > TAP_lbound;
        mpm.mpc.branch(idx_decr, TAP) = mpm.mpc.branch(idx_decr, TAP) - TAP_decr * TAP_stepsize;
        % Increase TAP ratios, check upper bound
        TAP_incr = mpm.mpc.branch(idx_incr, TAP) < TAP_ubound;
        mpm.mpc.branch(idx_incr, TAP) = mpm.mpc.branch(idx_incr, TAP) + TAP_incr * TAP_stepsize;
        
        % Check if any TAP ratios were changed
        TAPs_changed = TAPs_changed | any(TAP_incr) | any(TAP_decr);
        
        % Calculate meaningful data if needed
        if mpm.oltc.verbose == true
            TAP_after = mpm.mpc.branch(idx_branches_in_oltc, TAP);
            TAP_position = TAP_position + double(TAP_before~=TAP_after);
        end
    end
end

% Print out warning, if max_iter is reached (max number of tap changes
% were made)
if i >= mpm.oltc.max_iter
    fprintf("Maximum number of OLTC-iteratios (tap-changes) reached: %d \n", mpm.oltc.max_iter);
end

% Print out meaningful data if needed
if mpm.oltc.verbose == true
    TAPs_changed
    %         [TAP_before'; TAP_after'];
    %         [mpm.mpc.bus(idx_bus_oltc_f,1) mpm.mpc.bus(idx_bus_oltc_t,1) BASEf BASEt TAP_before TAP_after];
    %         VMf1 = mpm.mpc.bus(idx_bus_oltc_f,VM);
    %         VMt1 = mpm.mpc.bus(idx_bus_oltc_t,VM);
    %         [VMf0 VMf1 VMt0 VMt1];
    TAP_position
end


