function [r, chgs] = make_opf_feasible(mpc, mpopt)
%MAKE_OPF_FEASIBLE Attempts to relax constraints to make an OPF feasible.
%   [RESULTS, CHGS] = MAKE_OPF_FEASIBLE(MPC, MPOPT) 
%
%   Attempts to automate the process of finding a feasible OPF solution when
%   starting with an infeasible case.
%
%   MPC - initial (possibly infeasible) MATPOWER case struct
%   MPOPT - (optional) MATPOWER options struct.
%
%   Work-in-progress:   CHGS are currently not returned.
%                       At this point, the code *is* the documentation.
%
%   See also RUNOPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
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

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%%-----  process inputs  -----
if nargin < 2
    mpopt = mpoption();
end

%% initialize
msg = '';
done = 0;
rc = 0;     %% return code

%%-----  on-line generation capacity <= fixed load  -----
ig = find(~isload(mpc.gen) & mpc.gen(:, GEN_STATUS) > 0);
total_Pmax = sum(mpc.gen(ig, PMAX));
total_Pd   = total_load(mpc.bus, [], 'all');
if total_Pmax <= total_Pd
    msg = sprintf('%stotal fixed load (%g MW) > total on-line generation capacity (%g MW)\n', ...
        msg, total_Pd, total_Pmax);
    r = mpc;
    r.success = 0;
    r.f = NaN;
    rc = -2;        %% inadequate on-line generation capacity
    done = 1;
end

%%-----  check for infeasible limits  -----
if ~done
    %% generator P limits
    g = find(mpc.gen(:, PMIN) > mpc.gen(:, PMAX) & mpc.gen(:, GEN_STATUS) > 0);
    if ~isempty(g)
        for k = 1:length(g)
            msg = sprintf('%sPmax (%g) < Pmin (%g) for generator %d at bus %d, ', ...
                msg, mpc.gen(g(k), PMAX), mpc.gen(g(k), PMIN), g(k), mpc.gen(g(k), GEN_BUS));
            %% decide which limit to move based on where Pg is
            if mpc.gen(g(k), PG) >= mpc.gen(g(k), PMIN)
                mpc.gen(g(k), PMAX) = mpc.gen(g(k), PMIN);      %% move Pmax
                msg = sprintf('%smove Pmax\n', msg);
            elseif mpc.gen(g(k), PG) <= mpc.gen(g(k), PMAX)
                mpc.gen(g(k), PMAX) = mpc.gen(g(k), PMIN);      %% move Pmin
                msg = sprintf('%smove Pmax\n', msg);
            else
                tmp = mpc.gen(g(k), PMAX);                      %% swap 'em
                mpc.gen(g(k), PMAX) = mpc.gen(g(k), PMIN);
                mpc.gen(g(k), PMIN) = tmp;
                msg = sprintf('%sswap Pmin and Pmax\n', msg);
            end
        end
        rc = 2;         %% fixed generation capacity limits
    end
    
    %% generator Q limits
%     rc = 3;         %% fixed reactive generation capacity limits

    %% voltage angle limits
%     rc = 4;         %% fixed voltage angle limits

    %% voltage magnitude limits
%     rc = 5;         %% fixed voltage magnitude limits

    %% negative branch flow limits
%     rc = 6;         %% fixed negative branch flow limits

    %% run initial (infeasible?) case
    r = runopf(mpc, mpopt);
    if r.success
        done = 1;
    end
end

%%-----  attempt with short term branch ratings  -----
if ~done
    %% get branch limits
    rate_a = mpc.branch(:, RATE_A);
    rate_b = mpc.branch(:, RATE_B);
    rate_c = mpc.branch(:, RATE_C);
    rate_a(rate_a == 0) = Inf;
    rate_b(rate_b == 0) = Inf;
    rate_c(rate_c == 0) = Inf;

    %% set short term limits
    rating = max([rate_a rate_b], [], 2);
    rating(isinf(rating)) = 0;
    mpc.branch(:, RATE_A) = rating;
    r = runopf(mpc, mpopt);
    if r.success
        done = 1;
        rc = 7;         %% using short-term branch ratings
        msg = sprintf('%susing short-term branch ratings\n', msg);
    end
end

%%-----  attempt with emergency branch ratings  -----
if ~done
    %% set emergency limits
    rating = max([rate_a rate_b rate_c], [], 2);
    rating(isinf(rating)) = 0;
    mpc.branch(:, RATE_A) = rating;
    r = runopf(mpc, mpopt);

    if r.success
        done = 1;
        rc = 8;         %% using emergency branch ratings
        msg = sprintf('%susing emergency branch ratings\n', msg);
    end
end

%%-----  attempt without branch flow limits  -----
if ~done
    %% try with no line limits
    mpc.branch(:, RATE_A) = 0;
    r = runopf(mpc, mpopt);
                
    if ~r.success
        done = 1;
        msg = sprintf('%sstill infeasible without line limits\n', msg);
    else
        %% the following lines should be eliminated when I do further
        %% relaxations
        msg = sprintf('%susing no branch limits\n', msg);
        rc = 9;         %% using no branch ratings
    end
end

%%-----  try again, relaxing only violated limits  -----
if ~done

end

%% include msg in output
r.msg = msg;
r.rc  = rc;
if mpopt.verbose
    fprintf('%s', r.msg);
end
