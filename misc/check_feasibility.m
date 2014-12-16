function [v, f, hn, gn, Al, Au, xl, xu] = check_feasibility(mpc, mpopt)
%CHECK_FEASIBILITY  Returns the maximum constraint violation in p.u.
%   V = CHECK_FEASIBILITY(MPC)
%   [V, F, HN, GN, AL, AU, XL, XU] = CHECK_FEASIBILITY(MPC, MPOPT)
%
%   Not thoroughly tested.

%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
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

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

if nargin < 2
    mpopt = mpoption;
end

[mpc, mpopt] = opf_args(mpc, mpopt);
mpc = ext2int(mpc);
om = opf_setup(mpc, mpopt);
om = build_cost_params(om);

%% unpack data
mpc = get_mpc(om);
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(mpc.bus, 1);      %% number of buses
nl = size(mpc.branch, 1);   %% number of branches
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs

%% linear constraints
[A, l, u] = linear_constraints(om);

%% bounds on optimization vars
[x, xmin, xmax] = getv(om);

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);

%% set y variables (ONLY IMPLEMENTED FOR PG)
if ny > 0
    ipwl = find(mpc.gencost(:, MODEL) == PW_LINEAR);
    ig = vv.i1.Pg:vv.iN.Pg;
    x(vv.i1.y:vv.iN.y) = totcost(mpc.gencost(ipwl, :), x(ig(ipwl)) * mpc.baseMVA);
end

%% find branches with flow limits
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

%%-----  run opf  -----
f_fcn = @(x)opf_costfcn(x, om);
gh_fcn = @(x)opf_consfcn(x, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
hess_fcn = @(x, lambda, cost_mult)opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);

[f, df, d2f] = f_fcn(x);
[hn, gn, dhn, dgn] = gh_fcn(x);

Ax = A * x;
Au = Ax - u;
Al = l - Ax;

xu = x - xmax;
xl = xmin - x;

v = max([hn; abs(gn); Au; Al; xu; xl]);
