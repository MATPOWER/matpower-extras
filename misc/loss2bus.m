function bus_loss = loss2bus(mpc)
%LOSS2BUS   Accumulates branch real power losses at downstream buses.
%
%   BUS_LOSS = LOSS2BUS(MPC)
%
%   Takes a solved AC power flow case as input and returns an
%   NB x 1 vector containing all branch active power losses
%   accumulated to the bus at the downstream end of the branch.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman
%   Copyright (c) 2014 by Ray Zimmerman
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

define_constants;
    
%% create external to internal bus map
nb  = size(mpc.bus, 1);     %% number of buses
nl  = size(mpc.branch, 1);  %% number of branches
mb  = max(abs([mpc.bus(:, BUS_I); mpc.gen(:, GEN_BUS); ...
            mpc.branch(:, F_BUS); mpc.branch(:, T_BUS)]));
e2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, mb, 1);

%% assign losses to downstream buses
loss = mpc.branch(:, PF) + mpc.branch(:, PT);
bus_loss = zeros(nb, 1);
for j = 1:nl    %% need to use loop to accumulate for multiple lines per bus
    if mpc.branch(j, PF) >= mpc.branch(j, PT)
        b = T_BUS;
    else
        b = F_BUS;
    end
    ib = e2i(mpc.branch(j, b));
    bus_loss(ib) = bus_loss(ib) + loss(j);
end
