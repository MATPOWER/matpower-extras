function [co, cb, bus, gen, branch, f, dispatch, success] = ...
            smartmkt(mpc, offers, bids, mkt, mpopt)
%SMARTMKT  Runs the PowerWeb smart market.
%   [co, cb, bus, gen, branch, f, dispatch, success] = smartmkt(mpc, ...
%		offers, bids, mkt, mpopt) runs the ISO smart market.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2006 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialization  -----
%% default arguments
if nargin < 5
    mpopt = mpoption;       %% use default options
end

%% options
verbose = mpopt(31);

%% initialize some stuff
G = find( ~isload(mpc.gen) );       %% real generators
L = find(  isload(mpc.gen) );       %% dispatchable loads
if isfield(offers, 'Q') | isfield(bids, 'Q')
    haveQ = 1;
else
    haveQ = 0;
end

if haveQ & ~isempty(L) & mkt.auction_type ~= 0 & ...
        mkt.auction_type ~= 1 & mkt.auction_type ~= 5
    error(['Combined active/reactive power markets with constant power factor ', ...
            'dispatchable loads are only implemented for auction types 0, 1 & 5']);
end

%% set power flow formulation based on market
mpopt = mpoption(mpopt, 'PF_DC', strcmp(mkt.OPF, 'DC'));

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[AREA_I, PRICE_REF_BUS] = idx_area;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[QUANTITY, PRICE, FCOST, VCOST, SCOST, PENALTY] = idx_disp;

%% set up cost info & generator limits
mkt.lim = pricelimits(mkt.lim, isfield(offers, 'Q') | isfield(bids, 'Q'));
[gen, genoffer] = off2case(mpc.gen, mpc.gencost, offers, bids, mkt.lim);

%% move Pmin and Pmax limits out slightly to avoid problems
%% with lambdas caused by rounding errors when corner point
%% of cost function lies at exactly Pmin or Pmax
if any(find(genoffer(:, MODEL) == PW_LINEAR))
    gg = find( ~isload(gen) );      %% skip dispatchable loads
    gen(gg, PMIN) = gen(gg, PMIN) - 100 * mpopt(16) * ones(size(gg));
    gen(gg, PMAX) = gen(gg, PMAX) + 100 * mpopt(16) * ones(size(gg));
end

%%-----  solve the optimization problem  -----
%% attempt OPF
[bus, gen, branch, f, success, et] =  uopf(mpc.baseMVA, mpc.bus, gen, ...
            genoffer, mpc.branch, mpc.areas, mpopt);
if verbose & ~success
    fprintf('\nSMARTMARKET: non-convergent UOPF');
end

%%-----  compute quantities, prices & costs  -----
%% compute quantities & prices
ng = size(gen, 1);
if success      %% OPF solved case fine
    %% get nodal marginal prices from OPF
    gbus    = gen(:, GEN_BUS);                     %% indices of buses w/gens
    npP 	= max([ size(offers.P.qty, 2) size(bids.P.qty, 2) ]);
    lamP    = diag(bus(gbus, LAM_P)) * ones(ng, npP);    %% real power prices
    lamQ    = diag(bus(gbus, LAM_Q)) * ones(ng, npP);    %% reactive power prices
    
    %% compute fudge factor for lamP to include price of bundled reactive power
    pf   = zeros(ng, npP);                        %% for loads Q = pf * P
    Qlim =  (gen(L, QMIN) == 0) .* gen(L, QMAX) + ...
            (gen(L, QMAX) == 0) .* gen(L, QMIN);
    pf(L) = Qlim ./ gen(L, PMIN);

    gtee_prc.offer = 1;         %% guarantee that cleared offers are >= offers
    Poffer = offers.P;
    Poffer.lam = lamP(G,:);
    Poffer.total_qty = gen(G, PG);
    
    Pbid = bids.P;
    Pbid.total_qty = -gen(L, PG);
    if haveQ
        Pbid.lam = lamP(L,:);   %% use unbundled lambdas
        gtee_prc.bid = 0;       %% allow cleared bids to be above bid price
    else
        Pbid.lam = lamP(L,:) + diag(pf(L)) * lamQ(L,:); %% used bundled lambdas
        gtee_prc.bid = 1;       %% guarantee that cleared bids are <= bids
    end

    [co.P, cb.P] = auction(Poffer, Pbid, mkt.auction_type, mkt.lim.P, gtee_prc);

    if haveQ
        npQ = max([ size(offers.Q.qty, 2) size(bids.Q.qty, 2) ]);
        
        %% get nodal marginal prices from OPF
        lamP    = diag(bus(gbus, LAM_P)) * ones(ng, npQ);    %% real power prices
        lamQ    = diag(bus(gbus, LAM_Q)) * ones(ng, npQ);    %% reactive power prices

%         %% compute fudge factor for lamP to include price of bundled reactive power
%         pf   = zeros(ng, npQ);                        %% for loads P = pf * Q
%         Qlim =  (gen(L, QMIN) == 0) .* gen(L, QMAX) + ...
%                 (gen(L, QMAX) == 0) .* gen(L, QMIN);
%         kk = find(Qlim);
%         pf(L(kk)) = gen(L(kk), PMIN) ./ Qlim(kk);
    
        Qoffer = offers.Q;
        Qoffer.lam = lamQ;      %% use unbundled lambdas
%         Qoffer.lam(L,:) = lamQ(L,:) + diag(pf(L)) * lamP(L,:);
        Qoffer.total_qty = (gen(:, QG) > 0) .* gen(:, QG);
        
        Qbid = bids.Q;
        Qbid.lam = lamQ;        %% use unbundled lambdas
%         Qbid.lam(L,:) = lamQ(L,:) + diag(pf(L)) * lamP(L,:);
        Qbid.total_qty = (gen(:, QG) < 0) .* -gen(:, QG);

        [co.Q, cb.Q] = auction(Qoffer, Qbid, mkt.auction_type, mkt.lim.Q, gtee_prc);
    end

    quantity    = gen(:, PG);
    price       = zeros(ng, 1);
    price(G)    = offers.P.prc(:, 1);   %% need these for prices for
    price(L)    = bids.P.prc(:, 1);     %% gens that are shut down
    if npP == 1
        k = find( co.P.qty );
        price(G(k)) = co.P.prc(k, :);
        k = find( cb.P.qty );
        price(L(k)) = cb.P.prc(k, :);
    else
        k = find( sum( co.P.qty' )' );
        price(G(k)) = sum( co.P.qty(k, :)' .* co.P.prc(k, :)' )' ./ sum( co.P.qty(k, :)' )';
        k = find( sum( cb.P.qty' )' );
        price(L(k)) = sum( cb.P.qty(k, :)' .* cb.P.prc(k, :)' )' ./ sum( cb.P.qty(k, :)' )';
    end
else        %% did not converge even with imports
    quantity    = zeros(ng, 1);
    price       = mkt.lim.P.max_offer * ones(ng, 1);
    co.P.qty = zeros(size(offers.P.qty));
    co.P.prc = zeros(size(offers.P.prc));
    cb.P.qty = zeros(size(bids.P.qty));
    cb.P.prc = zeros(size(bids.P.prc));
    if haveQ
        co.Q.qty = zeros(size(offers.Q.qty));
        co.Q.prc = zeros(size(offers.Q.prc));
        cb.Q.qty = zeros(size(bids.Q.qty));
        cb.Q.prc = zeros(size(bids.Q.prc));
    end
end


%% compute costs in $ (note, NOT $/hr)
fcost   = mkt.t * totcost(mpc.gencost, zeros(ng, 1) );      %% fixed costs
vcost   = mkt.t * totcost(mpc.gencost, quantity     ) - fcost;  %% variable costs
scost   =   (~mkt.u0 & gen(:, GEN_STATUS) >  0) .* ...
                mpc.gencost(:, STARTUP) + ...               %% startup costs
            ( mkt.u0 & gen(:, GEN_STATUS) <= 0) .* ...
                mpc.gencost(:, SHUTDOWN);                   %% shutdown costs

%% store in dispatch
dispatch = zeros(ng, PENALTY);
dispatch(:, [QUANTITY PRICE FCOST VCOST SCOST]) = [quantity price fcost vcost scost];

return;
