function [co, cb] = auction(offers, bids, auction_type, limit_prc, gtee_prc)
%AUCTION  Clear auction based on OPF results (qty's and lambdas).
%   [co, cb] = auction(offers, bids, auction_type, limit_prc, gtee_prc)
%   Clears a set of bids and offers based on the results of an OPF, where the
%   pricing is adjusted for network losses and binding constraints.
%   The arguments offers and bids are structs with the following fields:
%       qty       - m x n, offer/bid quantities, m offers/bids, n blocks
%       prc       - m x n, offer/bid prices
%       lam       - m x n, corresponding lambdas
%       total_qty - m x 1, total quantity cleared for each offer/bid
%
%   There are 8 types of auctions implemented, specified by auction_type.
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
%   Whether or not cleared offer (bid) prices are guaranteed to be greater
%   (less) than or equal to the corresponding offer (bid) price is specified by
%   a flag gtee_prc.offer (gtee_prc.bid). The default is value true.
%   
%   Offer/bid and cleared offer/bid min and max prices are specified in the
%   limit_prc struct with the following fields:
%       max_offer
%       min_bid
%       max_cleared_offer
%       min_cleared_bid
%   Offers (bids) above (below) max_offer (min_bid) are treated as withheld
%   and cleared offer (bid) prices above (below) max_cleared_offer
%   (min_cleared_bid) are clipped to max_cleared offer (min_cleared_bid) if
%   given. All of these limit prices are ignored if the field is missing
%   or is empty.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialization  -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% initialize some stuff
zero_tol = 1e-5;
% zero_tol = 0.1;   %% fmincon is SO bad with prices that it is
                    %% NOT recommended for use with auction.m
big_num = 1e6;
if isempty(bids)
    bids = struct(  'qty', [], ...
                    'prc', [], ...
                    'lam', [], ...
                    'total_qty', [] );
end
if nargin < 4 || isempty(limit_prc)
    limit_prc = struct( 'max_offer', [], 'min_bid', [], ...
                        'max_cleared_offer', [], 'min_cleared_bid', [] );
else
    if ~isfield(limit_prc, 'max_offer'),         limit_prc.max_offer = [];         end
    if ~isfield(limit_prc, 'min_bid'),           limit_prc.min_bid = [];           end
    if ~isfield(limit_prc, 'max_cleared_offer'), limit_prc.max_cleared_offer = []; end
    if ~isfield(limit_prc, 'min_cleared_bid'),   limit_prc.min_cleared_bid = [];   end
end

if nargin < 5 || isempty(gtee_prc)
    gtee_prc = struct(  'offer', 1, 'bid', 1    );
else
    if ~isfield(gtee_prc, 'offer'), gtee_prc.offer = 1; end
    if ~isfield(gtee_prc, 'bid'),   gtee_prc.bid = 1;   end
end

[nro, nco] = size(offers.qty);
[nrb, ncb] = size(bids.qty);

%% determine cleared quantities
if isempty(limit_prc.max_offer)
    [co.qty, o.on, o.off] = clear_qty(offers.qty, offers.total_qty);
else
    mask = offers.prc <= limit_prc.max_offer;
    [co.qty, o.on, o.off] = clear_qty(offers.qty, offers.total_qty, mask);
end
if isempty(limit_prc.min_bid)
    [cb.qty, b.on, b.off] = clear_qty(bids.qty, bids.total_qty);
else
    mask = bids.prc <= limit_prc.min_bid;
    [cb.qty, b.on, b.off] = clear_qty(bids.qty, bids.total_qty, mask);
end

%% initialize cleared prices
co.prc  = zeros(nro, nco);              %% cleared offer prices
cb.prc  = zeros(nrb, ncb);              %% cleared bid prices

%%-----  compute shift values to add to lam to get desired pricing  -----
%% The locationally adjusted offer/bid price, when normalized to an arbitrary
%% reference location where lambda is equal to ref_lam, is:
%%      norm_prc = prc + (ref_lam - lam)
%%      norm_prc = prc * (ref_lam / lam)
%% Then we can define the difference between the normalized offer/bid prices
%% and the ref_lam as:
%%      diff = norm_prc - ref_lam = prc - lam
%% This diff represents the gap between the marginal unit (setting lambda)
%% and the offer/bid price in question.
offer_diff = offers.prc - offers.lam;
bid_diff   = bids.prc   - bids.lam;

%% shift.LAO + lambda is equal to the last accepted offer
shift.LAO = o.on .* offer_diff - o.off * big_num;
shift.LAO( shift.LAO(:) > zero_tol ) = -big_num;    %% don't let gens @ Pmin set price
shift.LAO = max( shift.LAO(:) );

%% shift.FRO + lambda is equal to the first rejected offer
shift.FRO = o.off .* offer_diff + o.on * big_num;
shift.FRO = min( shift.FRO(:) );

if nrb
    %% shift.LAB + lambda is equal to the last accepted bid
    shift.LAB = b.on .* bid_diff + b.off * big_num;
    shift.LAB = min( shift.LAB(:) );
    
    %% shift.FRB + lambda is equal to the first rejected bid
    shift.FRB = b.off .* bid_diff - b.on * big_num;
    shift.FRB = max( shift.FRB(:) );
else
    shift.LAB = big_num;
    shift.FRB = shift.LAB;
end

%% cleared offer/bid prices for different auction types
if auction_type == 0        %% discriminative
    co.prc = offers.prc;
    cb.prc = bids.prc;
elseif auction_type == 1    %% LAO
    co.prc = offers.lam + shift.LAO;
    cb.prc = bids.lam   + shift.LAO;
elseif auction_type == 2    %% FRO
    co.prc = offers.lam + shift.FRO;
    cb.prc = bids.lam   + shift.FRO;
elseif auction_type == 3    %% LAB
    co.prc = offers.lam + shift.LAB;
    cb.prc = bids.lam   + shift.LAB;
elseif auction_type == 4    %% FRB
    co.prc = offers.lam + shift.FRB;
    cb.prc = bids.lam   + shift.FRB;
elseif auction_type == 5    %% 1st price
    co.prc = offers.lam;
    cb.prc = bids.lam;
elseif auction_type == 6    %% 2nd price
    if abs(shift.LAO) < zero_tol
        co.prc = offers.lam + min(shift.FRO,shift.LAB);
        cb.prc = bids.lam   + min(shift.FRO,shift.LAB);
    else
        co.prc = offers.lam + max(shift.LAO,shift.FRB);
        cb.prc = bids.lam   + max(shift.LAO,shift.FRB);
    end
elseif auction_type == 7    %% split the difference
    co.prc = offers.lam + (shift.LAO + shift.LAB) / 2;
    cb.prc = bids.lam   + (shift.LAO + shift.LAB) / 2;
elseif auction_type == 8    %% LAO seller, LAB buyer
    co.prc = offers.lam + shift.LAO;
    cb.prc = bids.lam   + shift.LAB;
end

%% guarantee that cleared offer prices are >= offers
if gtee_prc.offer
    clip = o.on .* (offers.prc - co.prc);
    co.prc = co.prc + (clip > zero_tol) .* clip;
end

%% guarantee that cleared bid prices are >= bids
if gtee_prc.bid
    clip = b.on .* (bids.prc - cb.prc);
    cb.prc = cb.prc + (clip < -zero_tol) .* clip;
end

%% clip cleared offer prices by limit_prc.max_cleared_offer
if ~isempty(limit_prc.max_cleared_offer)
    co.prc = co.prc + (co.prc > limit_prc.max_cleared_offer) .* ...
        (limit_prc.max_cleared_offer - co.prc);
end

%% clip cleared bid prices by limit_prc.min_cleared_bid
if ~isempty(limit_prc.min_cleared_bid)
    cb.prc = cb.prc + (cb.prc < limit_prc.min_cleared_bid) .* ...
        (limit_prc.min_cleared_bid - co.prc);
end

%% make prices uniform after clipping (except for discrim auction)
%% since clipping may only affect a single block of a multi-block generator
if auction_type ~= 0
    %% equal to largest price in row
    if nco > 1
        co.prc = diag(max(co.prc')) * ones(nro,nco);
    end
    if ncb > 1
        cb.prc = diag(min(cb.prc')) * ones(nrb,ncb);
    end
end


function [cqty, on, off] = clear_qty(qty, total_cqty, mask)
%CLEAR_QTY  Computed cleared offer/bid quantities from totals.
%  Inputs:
%   qty        - m x n, offer/bid quantities, m offers/bids, n blocks
%   total_cqty - m x 1, total cleared quantity for each offer/bid
%   mask       - m x n, boolean indicating which offers/bids are valid (not withheld)
%  Outputs:
%   cqty       - m x n, cleared offer/bid quantities, m offers/bids, n blocks
%   on         - m x n, 1 if partially or fully accepted, 0 if rejected
%   off        - m x n, 1 if rejected, 0 if partially or fully accepted

[nr, nc] = size(qty);
accept  = zeros(nr,nc);
cqty    = zeros(nr,nc);
for i = 1:nr            %% offer/bid i
    for j = 1:nc        %% block j
        if qty(i, j)    %% ignore zero quantity offers/bids
            %% compute fraction of the block accepted ...
            accept(i, j) = (total_cqty(i) - sum(qty(i, 1:j-1))) / qty(i, j);
            %% ... clipped to the range [0, 1]  (i.e. 0-100%)
            if accept(i, j) > 1
                accept(i, j) = 1;
            elseif accept(i, j) < 1.0e-5
                accept(i, j) = 0;
            end
            cqty(i, j) = qty(i, j) * accept(i, j);
        end
    end
end

if nargin == 3
    accept = mask .* accept;
end

on  = (accept  > 0);
off = (accept == 0);
