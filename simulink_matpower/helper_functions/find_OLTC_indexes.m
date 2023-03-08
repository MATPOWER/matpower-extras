function [idx_branches_in_oltc] = find_OLTC_indexes(mpc)

%% This function takes the Matpower case file and finds the indexes of the branches that are transformers with on load tap changers (OLTCs)
% The indexes are used by the OLTC function to model the behavior of OLTCs

% CHANGE THIS CODE TO FIND THE INDEXES OF YOUR TAP CHANGERS


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
% A branch is connecting two buses with different voltage level that branch is a transformer
% In my case transformers between the following voltage levels have OLTCs

% 225 kV & 90 kV
% 225 kV & 63 kV
% 380 kV & 63 kV 
% 380 kV & 90 kV

%%

nmb_branches = size(mpc.branch,1); % number of branches in the case file

idx_branches_in_oltc = [];

for i=1:nmb_branches % go through all the branches
    from_bus = mpc.branch(i,1); % What is the bus the branch is starting from
    to_bus = mpc.branch(i,2);   % What is the bus the branch is ending at
    from_bus_index = find(mpc.bus(:,1)==from_bus); % find the index of the from_bus
    to_bus_index = find(mpc.bus(:,1)==to_bus); % find the index of the to_bus
    % with the indexes we can check the voltage levels of the from and to
    % bus and find out which branches are transformers with OLTCs
%     if isempty(to_bus_index) || isempty(from_bus_index)
%         continue
%     end
if ~isempty(from_bus_index) & ~isempty(to_bus_index)
    if mpc.bus(from_bus_index,10) == mpc.bus(to_bus_index,10)
        % the voltages level are equal => no transformer and therefore no
        % OLTC
    elseif ((mpc.bus(from_bus_index,10) == 225 && mpc.bus(to_bus_index,10) == 90) || (mpc.bus(from_bus_index,10)== 90 && mpc.bus(to_bus_index,10)==225))
        % One bus is 225 and one is 90. Or one bus is 90 and the other one is 225.
        % Therefore we have a transformer with an OLTC
        idx_branches_in_oltc = [idx_branches_in_oltc; i];
    elseif ((mpc.bus(from_bus_index,10) == 225 && mpc.bus(to_bus_index,10) == 63) || (mpc.bus(from_bus_index,10)== 63 && mpc.bus(to_bus_index,10)==225))
        idx_branches_in_oltc = [idx_branches_in_oltc; i];
    elseif ((mpc.bus(from_bus_index,10) == 380 && mpc.bus(to_bus_index,10) == 63) || (mpc.bus(from_bus_index,10)== 63 && mpc.bus(to_bus_index,10)==380))
        idx_branches_in_oltc = [idx_branches_in_oltc; i];
    elseif ((mpc.bus(from_bus_index,10) == 380 && mpc.bus(to_bus_index,10) == 90) || (mpc.bus(from_bus_index,10)== 90 && mpc.bus(to_bus_index,10)==380))
        idx_branches_in_oltc = [idx_branches_in_oltc; i];
    else
        % nothing to do. This branch is a line or there is no tap changer on this transformer
    end
end
    


end


end

