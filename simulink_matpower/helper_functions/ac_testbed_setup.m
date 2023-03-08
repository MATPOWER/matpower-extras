%% Initialisation of the Simulink simulation

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


%% From MATPOWER: Defines all the field names within the MATPOWER case struct
define_constants;

%% Load MATPOWER case from file
mpm.mpc = loadcase(case39); % used for the example

%% Configure MATPOWER options to reduce console clutter
mpm.mpopt = mpoption(mpoption, 'verbose', 0); % avoid printing during iteration
mpm.mpopt.out.all=0; % avoid printing reports after simulations
mpm.mpopt.pf.tool=10^-12;

%% Enforce reactive power constraints
mpm.mpopt.pf.enforce_q_lims = 1;

%% presolve the grid once
mpm.mpc = runpf(mpm.mpc,mpm.mpopt);

%% Set up OLTC configuration
mpm.oltc.enable = 0;
mpm.oltc.V_deadband = 0.015;
mpm.oltc.TAP_stepsize = 0.00625;
mpm.oltc.verbose = false;
mpm.oltc.max_iter = 100;

%% Define Inputs for Simulation
std_load = mpm.mpc.bus(:,PD) + 1j*mpm.mpc.bus(:,QD); %In MVA
std_gen = mpm.mpc.gen(:,PG) + 1j*mpm.mpc.gen(:,QG); %In MVA


%% Get the indexes of the branches that are transformers with OLTCs
        
% Get the indexes of all OLTCs and not just the Blocaux area
idx_branches_in_oltc = find_OLTC_indexes(mpm.mpc);
nmb_idx_branches_in_oltc = length(idx_branches_in_oltc);
zeros_num_tap_ratios = zeros(nmb_idx_branches_in_oltc,1);

%% Get the number of generators
number_generators = length(mpm.mpc.gen(:,1));
zeros_num_gens = zeros(size(mpm.mpc.gen(:,1))); % needed to initialize the output of get_gens otherwise Simulink doesn't play along.

%% Get the number of branches
number_branches = length(mpm.mpc.branch(:,1));
zeros_num_branches = zeros(size(mpm.mpc.branch(:,1))); % needed to initialize the output of ACPF and get_buses otherwise Simulink doesn't play along.

%% Get the number of buses
number_buses = length(mpm.mpc.bus(:,1));
zeros_num_buses = zeros(size(mpm.mpc.bus(:,1))); % needed to initialize the output of get_buses otherwise Simulink doesn't play along.

