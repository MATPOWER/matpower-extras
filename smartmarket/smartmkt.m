function [cq, cp, bus, gen, branch, f, dispatch, success] = ...
            smartmkt(baseMVA, bus, gen, gencost, branch, areas, q, p, mkt, max_p, u0, t, mpopt)
%SMARTMKT  Runs the PowerWeb smart market.
%   [cq, cp, bus, gen, branch, f, dispatch, success] = smartmkt(baseMVA, bus, gen, ...
%   branch, areas, gencost, q, p, max_p, u0, t, mpopt) runs the ISO smart market.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialization  -----
%% default arguments
if nargin < 13
    mpopt = mpoption;       %% use default options
end

%% options
verbose = mpopt(31);
margin = 1.05;          %% must have 5% capacity margin if considering losses

%% parse market code
code        = mkt - 1000;
adjust4loc  = fix(code/100);    code = rem(code, 100);
auction_type= fix(code/10);

%% set power flow formulation based on market
mpopt = mpoption(mpopt, 'PF_DC', adjust4loc == 2);

if adjust4loc ~= 1
    margin = 1;         %% no margin needed without a network or with loss-less power flow
    if adjust4loc == 0
        error('The non-network version of the smart market has not yet been implemented.');
    end
end
success = 0;

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
    GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[AREA_I, PRICE_REF_BUS] = idx_area;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;
[QUANTITY, PRICE, FCOST, VCOST, SCOST, PENALTY] = idx_disp;

%% initialize dispatch
[ng, np] = size(q);
dispatch = zeros(ng, PENALTY);

%% find import generators
imports = find(gen(:, GEN_STATUS) == -1);

if ~isempty(imports)
	%% modify p so that imports come in at reservation price + gap
	%% or, alternatively, max offered price + gap
	%% and all actual offers above reservation price are kept out
	
	%% find max offer price less than max_p
	gap = 5;        %% set the gap between reservation price and max import price
	on = find(gen(:, GEN_STATUS) > 0);  %% which generators are on?
	on_p = p(on,:);                     %% get submatrix to do vector mods
	in  = find(on_p <= max_p);          %% which blocks are in?
	out = find(on_p >  max_p);          %% which blocks are out?
	max_offered_p = max(on_p(in));      %% maximum offered price
	on_p(out) = on_p(out) + gap;        %% make sure blocks above max_p are still
										%% eliminated by reservation price + gap
	p(on, :) = on_p;                    %% put back submatrix
	
	%% set price at which imports come in
	p(imports, :) = ones(length(imports), np) * (max_offered_p + gap);
	% p(imports, :) = ones(length(imports), np) * (max_p + gap);
	
	%% bump up max price
	max_p = max_p + gap;
end

%% set up cost info & generator limits
%% eliminates offers (but not bids) above max_p
[gen, genoffer] = off2case(gen, gencost, q, p, max_p);

%% set import cost equal to import offers (replace zero placeholders)
if ~isempty(imports)
    tmp = min([ size(gencost, 2) size(genoffer, 2) ]);
    gencost(imports, 1:tmp) = genoffer(imports, 1:tmp);
end

%% make copy of initial bus voltages and generator outputs
bus0 = bus(:, [VM, VA]);
gen0 = gen(:, [PG, QG, VG]);

%% compute total load capacity
load_capacity   = sum(bus(:, PD));

%% move Pmin and Pmax limits out slightly to avoid problems
%% with lambdas caused by rounding errors when corner point
%% of cost function lies at exactly Pmin or Pmax
if any(find(genoffer(:, MODEL) == PW_LINEAR))
    gg = find( gen(:, PMAX) > 0 );      %% skip variable loads
    gen(gg, PMIN) = gen(gg, PMIN) - 100 * mpopt(16) * ones(size(gg));
    gen(gg, PMAX) = gen(gg, PMAX) + 100 * mpopt(16) * ones(size(gg));
end

%%-----  solve the optimization problem  -----
withreserves    = 0;
while 1
    %% check for sufficient generation
    on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
    gen_capacity    = sum(gen(on, PMAX));
    if gen_capacity < load_capacity * margin & ~withreserves
        %%-----  insufficient generation, try again with imports  -----
        success = 0;
        if verbose
            fprintf('\nSMARTMARKET: insufficient generation');
        end
    else
        %%-----  sufficient generation, attempt OPF  -----
        %% restore voltages and gen outputs from original case
        bus(:, [VM, VA]) = bus0;
        gen(:, [PG, QG, VG]) = gen0;
        
        %% attempt OPF
        [bus, gen, branch, f, success, et] =  uopf(baseMVA, bus, gen, genoffer, ...
                    branch, areas, mpopt);
        if verbose & ~success
            fprintf('\nSMARTMARKET: non-convergent UOPF');
        end
    end
    
    %% should we retry?
    if ~success & ~withreserves & ~isempty(imports)
        if verbose
            fprintf('\nSMARTMARKET: turning on import generators\n\n');
        end
        gen(imports, GEN_STATUS) = 2 * ones(size(imports));
        withreserves = 1;
    else
        break;
    end
end

%%-----  compute quantities, prices & costs  -----
%% compute quantities & prices
if success      %% OPF solved case fine
    quantity    = gen(:, PG);
    %% turn down verbosity one level for call to auction
    if verbose
        mpopt = mpoption(mpopt, 'VERBOSE', verbose-1);
    end
    [cq, cp] = auction(bus, gen, gencost, q, p, max_p, auction_type, mpopt);
    if size(cq, 2) == 1
        k = find( cq );
    else
        k = find( sum( cq' )' );
    end
    price = cp(:, 1);           %% need this for prices for gens that are shut down
    price(k) = sum( cq(k, :)' .* cp(k, :)' )' ./ sum( cq(k, :)' )';
else        %% did not converge even with imports
    quantity    = zeros(ng, 1);
    price       = max_p * ones(ng, 1);
end

%% compute costs in $ (note, NOT $/hr)
fcost   = t * totcost(gencost, zeros(ng, 1) );                          %% fixed costs
vcost   = t * totcost(gencost, quantity     ) - fcost;                  %% variable costs
scost   = (~u0 & gen(:, GEN_STATUS) >  0) .* gencost(:, STARTUP) + ...  %% startup costs
          ( u0 & gen(:, GEN_STATUS) <= 0) .* gencost(:, SHUTDOWN);      %% shutdown costs

%% store in dispatch
dispatch(:, [QUANTITY PRICE FCOST VCOST SCOST]) = [quantity price fcost vcost scost];

return;
