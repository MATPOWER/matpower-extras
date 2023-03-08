function [vm_out, va_out, pg_out, qg_out, pd_out,qd_out, pf_out,qf_out, pt_out,qt_out, success] ...
    = extrinsic_ACPF(vm_in, va_in, pg_in, qg_in, pd_in, qd_in)
%% This function is a wrapper function for MATPOWER's runpf so that it can
%% be used in Simulink. It also implements the analysis of OLTCs in steady
%% state. 

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
%%
    define_constants;
    %% Get mpm struct from base workspace
    mpm = evalin('base','mpm');
    idx_branches_in_oltc = evalin('base','idx_branches_in_oltc');
    %% Read values from Simulink
    mpm.mpc.bus(:,PD) = pd_in; %In MW
    mpm.mpc.bus(:,QD) = qd_in; %In MVAr
    mpm.mpc.bus(:,VM) = vm_in; %In p.u.
    mpm.mpc.bus(:,VA) = va_in; %In degrees
    mpm.mpc.gen(:,PG) = pg_in; %In MW
    mpm.mpc.gen(:,QG) = qg_in; %In MVAr
    
    
    %% Solve power flow
    % standard matpower solver
    if ~mpm.oltc.enable
        [mpm.mpc, success] = runpf(mpm.mpc,mpm.mpopt);
    % standard matpower solver + OLTC
    else
        [mpm, success] = OLTC(mpm,idx_branches_in_oltc);
    end

    %% Write values back to Simulink
    pg_out = mpm.mpc.gen(:,PG); %In MW
    qg_out = mpm.mpc.gen(:,QG); %In MVAr
    vm_out = mpm.mpc.bus(:,VM); %In p.u
    va_out = mpm.mpc.bus(:,VA); %In degrees
    pd_out = mpm.mpc.bus(:,PD); %In MW
    qd_out = mpm.mpc.bus(:,QD); %In MVAr
    pf_out = mpm.mpc.branch(:,PF); %In MW
    qf_out = mpm.mpc.branch(:,QF); %In MVAr
    pt_out = mpm.mpc.branch(:,PT); %In MW
    qt_out = mpm.mpc.branch(:,QT); %In MVAr
    %% Store mpm struct in base workspace
    assignin('base','mpm',mpm);
end