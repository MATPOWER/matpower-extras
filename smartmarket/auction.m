function [cq, cp] = auction(bus, gen, gencost, q, p, max_p, auction_type, mpopt)
%AUCTION  Clear auction based on OPF results.
%   [cq, cp] = auction(bus, gen, gencost, q, p, max_p, auction_type, mpopt)
%   Clears a set of bids and offers based on the results of an OPF, where the
%   pricing is adjusted for network losses and binding constraints. There are
%   8 types of auctions implemented here, specified by auction_type.
%
%      0 - discriminative pricing (price equal to offer or bid)
%      1 - last accepted offer auction
%      2 - first rejected offer auction
%      3 - last accepted bid auction
%      4 - first rejected bid auction
%      5 - first price auction (marginal unit, offer or bid, sets the price)
%      6 - second price auction (if offer is marginal, price set by
%             min(FRO,LAB), else max(FRB,LAO)
%      7 - split the difference pricing (set by last accepted offer & bid)
%      8 - LAO sets seller price, LAB sets buyer price
%
%   Cleared offer prices (but not bid prices) are clipped to max_p.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialization  -----
%% default parameters
if nargin < 8
    mpopt = mpoption;       %% use default options
    if nargin < 7
        auction_type = 1;       %% use 1st price by default
    end
end

%% options
verbose = mpopt(31);
if have_fcn('minopf')
    zero_tol = 1e-5;
else
    zero_tol = 0.1;     %% fmincon is SO bad with prices that it is
                        %% NOT recommended for use with auction.m
end
big_num = 1e6;

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
    GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[AREA_I, PRICE_REF_BUS] = idx_area;

%% initialize some stuff
[ng, np] = size(q);
cq      = zeros(size(q));                   %% cleared quantity offers
cp      = zeros(size(p));                   %% cleared price offers
in      = find(gen(:, GEN_STATUS) >= 0);    %% which gens are in the game (non-reserve)
                                            %% (NOTE: can include some that are off-line)
G       = find( gen(in, PMIN) >= 0 | gen(in, PMAX) > 0 );   %% real generators
L       = find( gen(in, PMIN) < 0 & gen(in, PMAX) <= 0 );   %% variable loads
qin     = q(in, :);             %% reduce quantity offers ...
pin     = p(in, :);             %% ... and price offers to those gens that are "in"
cqin    = cq(in, :);
cpin    = cp(in, :);

%% determine cleared quantities
Pg = gen(in, PG);           %% active generation from generators that are "in"
Pg(L) = -Pg(L);             %% make quantity positive for loads
accept = zeros(size(qin));
for i = 1:length(in)        %% generator in(i)
    for j = 1:np                %% block j
        if qin(i, j)                    %% ignore zero quantity offers
            %% compute fraction of the block accepted ...
            accept(i, j) = (Pg(i) - sum(qin(i, 1:j-1))) / qin(i, j);
            %% ... clipped to the range [0, 1]  (i.e. 0-100%)
            if accept(i, j) > 1
                accept(i, j) = 1;
            elseif accept(i, j) < 1.0e-5
                accept(i, j) = 0;
            end
            cqin(i, j) = qin(i, j) * accept(i, j);
        end
    end
end
on  = (accept  > 0);
off = (accept == 0);
all_ones  = ones(size(pin));
all_zeros = zeros(size(pin));

%% get nodal marginal prices from OPF
gbus    = gen(in, GEN_BUS);                     %% indices of buses w/gens
lamP    = diag(bus(gbus, LAM_P)) * all_ones;    %% real power prices
lamQ    = diag(bus(gbus, LAM_Q)) * all_ones;    %% reactive power prices

%% compute fudge factor for lamP to include price of bundled reactive power
nz          = find( gen(L, PG) );
pf          = all_zeros;                        %% for loads Q = pf * P
pf(L(nz))   = gen(L(nz), QG) ./ gen(L(nz), PG);
Qfudge      = all_zeros;
Qfudge(L,:) = diag(pf(L)) * lamQ(L,:);

%%-----  compute shift matrices to add to lamP to get desired pricing  -----
%% If refP is the lambda P at an arbitrary reference bus, and normP is
%% the set of offers/bids normalized to the refbus location with bids
%% normalized to include bundled reactive power as well, then
%%     gap = normP - refP
%% In other words, (gap + refP) can be used to form the normal
%% non-locational bid/offer stack. And since refP is just a constant
%% it doesn't actually change the shape of the stacks.
gap = pin - lamP - Qfudge;

DISC_P  = gap;
DISC_PQ = gap + Qfudge;

LAO = on(G,:) .* gap(G,:) - off(G,:) * big_num;
LAO( find(LAO(:) > zero_tol) ) = -big_num;  %% don't let gens @ Pmin set price
LAO = max( LAO(:) ) * all_ones;

FRO = off(G,:) .* gap(G,:) + on(G,:) * big_num;
FRO = min( FRO(:) ) * all_ones;

if ~isempty(L)
    LAB_P = on(L,:) .* gap(L,:) + off(L,:) * big_num;
    LAB_P = min( LAB_P(:) ) * all_ones;
    
    FRB_P = off(L,:) .* gap(L,:) - on(L,:) * big_num;
    FRB_P = max( FRB_P(:) ) * all_ones;
else
    LAB_P = big_num * all_ones;
    FRB_P = LAB_P;
end
LAB_PQ = LAB_P + Qfudge;
FRB_PQ = FRB_P + Qfudge;

%% generator and load price shifts for different auction types
if auction_type == 0
    G_shift = DISC_P;
    L_shift = DISC_PQ;
elseif auction_type == 1
    G_shift = LAO;
    L_shift = LAO + Qfudge;
elseif auction_type == 2
    G_shift = FRO;
    L_shift = FRO + Qfudge;
elseif auction_type == 3
    G_shift = LAB_P;
    L_shift = LAB_PQ;
elseif auction_type == 4
    G_shift = FRB_P;
    L_shift = FRB_PQ;
elseif auction_type == 5
    G_shift = all_zeros;
    L_shift = Qfudge;
elseif auction_type == 6
    if abs(LAO(1,1)) < zero_tol
        G_shift = min(FRO(1,1),LAB_P(1,1)) * all_ones;
        L_shift = min(FRO(1,1),LAB_P(1,1)) + Qfudge;
    else
        G_shift = max(LAO(1,1),FRB_P(1,1)) * all_ones;
        L_shift = max(LAO(1,1),FRB_P(1,1)) + Qfudge;
    end
elseif auction_type == 7
    G_shift = (LAO + LAB_P) / 2;
    L_shift = (LAO + LAB_P) / 2 + Qfudge;
elseif auction_type == 8
    G_shift = LAO;
    L_shift = LAB_PQ;
end

%% compute the prices
cpin(G,:) = lamP(G,:) + G_shift(G,:);
cpin(L,:) = lamP(L,:) + L_shift(L,:);

%% clip cleared prices by offers and bids
clip = on .* (pin - cpin);
cpin(G,:) = cpin(G,:) + (clip(G,:) > zero_tol) .* clip(G,:);
cpin(L,:) = cpin(L,:) + (clip(L,:) < -zero_tol) .* clip(L,:);

%% clip cleared offer prices by max_p
cpin(G,:) = cpin(G,:) + (cpin(G,:) > max_p) .* (max_p - cpin(G,:));

%% make prices uniform after clipping (except for discrim auction)
%% since clipping may only affect a single block of a multi-block generator
if auction_type ~= 0 & np > 1
    cpin(G,:) = diag(max(cpin(G,:)')) * all_ones(G,:);  %% equal to largest price in row
    cpin(L,:) = diag(min(cpin(L,:)')) * all_ones(L,:);  %% equal to smallest price in row
end

cq(in, :) = cqin;
cp(in, :) = cpin;

%%-----  plot offers, adjusted offers, prices  -----
if verbose
    %% compute locational adjustment for generators
    ref_p = max(lamP(:,1));
    % ref_p = max(bus(:, LAM_P));
    adjustment = ref_p - lamP(:,1);                     %% locational price adjustment
    
    %% form 1-d vectors of all valid offers & corresponding adjustments, acceptance
    i = find(qin & pin <= max_p);
    qq = qin(i);                            %% quantity vector
    pp = pin(i);                            %% price vector
    adj = diag(adjustment) * all_ones;      %% adjustment matrix
    adj = adj(i);                           %% adjustment vector
    adj_pp = pp + adj;                      %% adjusted offer price vector
    acc = accept(i);                        %% acceptance vector
    final_p = cpin;
    final_p = final_p(i);
    
    %% sort by adjusted price
    [junk, k] = sort(adj_pp);
    qq = qq(k);
    pp = pp(k);
    adj = adj(k);
    adj_pp = adj_pp(k);
    acc = acc(k);
    final_p = final_p(k);

    %% plot it
    subplot(1,1,1);
    xd = [0; cumsum(qq)];
    [xs, ys1] = stairs(xd, [pp; pp(length(pp))]);
    [xs, ys2] = stairs(xd, [adj_pp; adj_pp(length(adj_pp))]);
    [xs, ys3] = stairs(xd, [acc; acc(length(acc))]);
    [xs, ys4] = stairs(xd, [final_p; final_p(length(final_p))]);
    [xs1, ys5] = stairs([0 sum(Pg) sum(qq)], [0 100 100]);
    plot(xs, ys2, 'r--');       %% adjusted offers
    axis([0 1.1*max(xd) 0 max([110, 1.1*max(adj_pp)])]);
    hold on;
    plot(xs, ys1, 'y-.');       %% unadjusted offers
    plot(xs, 100*ys3, 'b:');    %% percent of block accepted
    plot(xs, ys4, 'g');         %% price paid
%   plot(xs1, ys5, 'b');
    hold off;
    drawnow;
end

return;
