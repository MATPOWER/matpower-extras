function [mpc_out, co, cb, f, dispatch, success, et] = ...
                runmarket(mpc, offers, bids, mkt, mpopt, fname, solvedcase)
%RUNMARKET  Runs PowerWeb-style smart market.
%
%   [mpc_out, co, cb, f, dispatch, success, et] = ...
%           runmarket(mpc, offers, bids, mkt, mpopt, fname, solvedcase)
%
%   Computes the new generation and price schedules (cleared offers and bids)
%   based on the offers and bids submitted. See 'help off2case' for a
%   description of the offers and bids arguments. mkt is a struct with the
%   following fields:
%       auction_type - market used for dispatch and pricing
%       t            - time duration of the dispatch period in hours
%       u0           - vector of gen commitment status from prev period
%       lim          - offer/bid/price limits (see 'help pricelimits')
%       OPF          - 'AC' or 'DC', default is 'AC'
%
%   mpopt is an optional MATPOWER options vector (see 'help mpoption' for
%   details). The values for the auction_type field are defined as follows:
%
%      0 - discriminative pricing (price equal to offer or bid)
%      1 - last accepted offer auction
%      2 - first rejected offer auction
%      3 - last accepted bid auction
%      4 - first rejected bid auction
%      5 - first price auction (marginal unit, offer or bid, sets the price)
%      6 - second price auction (if offer is marginal, then price is set
%                                   by min(FRO,LAB), if bid, then max(FRB,LAO)
%      7 - split the difference pricing (price set by last accepted offer & bid)
%      8 - LAO sets seller price, LAB sets buyer price
%
%   The default auction_type is 5, where the marginal block (offer or bid)
%   sets the price. The default lim sets no offer/bid or price limits. The
%   default previous commitment status u0 is all ones (assume everything was
%   running) and the default duration t is 1 hour. The results may
%   optionally be printed to a file (appended if the file exists) whose name
%   is given in fname (in addition to printing to STDOUT). Optionally
%   returns the final values of the solved case in mpc_out, the cleared
%   offers and bids in co and cb, the objective function value f, the old
%   style dispatch matrix, the convergence status of the OPF in success, and
%   the elapsed time et. If a name is given in solvedcase, the solved case
%   will be written to a case file in MATPOWER format with the specified
%   name.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialize  -----
%% default arguments
if nargin < 7
    solvedcase = '';                        %% don't save solved case
    if nargin < 6
        fname = '';                         %% don't print results to a file
        if nargin < 5
            mpopt = mpoption;               %% use default options
            if nargin < 4
                mkt = [];      %% use default market
                if nargin < 3
                    bids = struct([]);
                    if nargin < 2
                        offers = struct([]);
                        if nargin < 1
                            mpc = 'case9';  %% default data file is 'case9.m'
                        end
                    end
                end
            end
        end
    end
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% read data & convert to internal bus numbering
mpc = loadcase(mpc);
[i2e, mpc.bus, mpc.gen, mpc.branch] = ...
    ext2int(mpc.bus, mpc.gen, mpc.branch);

%% assign default arguments
if isempty(mkt)
    mkt = struct( 'OPF', [], 'auction_type', [], 'lim', [], 'u0', [], 't', []);
end
if ~isfield(mkt, 'OPF') || isempty(mkt.OPF)
    mkt.OPF = 'AC';         %% default OPF is AC
end
if ~isfield(mkt, 'auction_type') || isempty(mkt.auction_type)
    mkt.auction_type = 5;   %% default auction type is first price
end
if ~isfield(mkt, 'lim') || isempty(mkt.lim)
    mkt.lim = pricelimits([], isfield(offers, 'Q') || isfield(bids, 'Q'));
end
if ~isfield(mkt, 'u0') || isempty(mkt.u0)
    mkt.u0 = ones(size(mpc.gen, 1), 1); %% default for previous gen commitment, all on
end
if ~isfield(mkt, 't') || isempty(mkt.t)
    mkt.t = 1;              %% default dispatch duration in hours
end

%% if offers not defined, use gencost
if isempty(offers) || isempty(offers.P.qty)
    [q, p] = case2off(mpc.gen, mpc.gencost);

    %% find indices for gens and variable loads
    G = find( ~isload(mpc.gen) );   %% real generators
    L = find(  isload(mpc.gen) );   %% variable loads
    offers = struct( 'P', struct( 'qty', q(G, :), 'prc', p(G, :) ) );
    bids   = struct( 'P', struct( 'qty', q(L, :), 'prc', p(L, :) ) );
end

%% start the clock
t0 = clock;

%% run the market
[co, cb, bus, gen, branch, f, dispatch, success] = ...
        smartmkt(mpc, offers, bids, mkt, mpopt);

%% compute elapsed time
et = etime(clock, t0);

%% convert back to original bus numbering & print results
gencost = mpc.gencost;
baseMVA =  mpc.baseMVA;
[bus, gen, branch] = int2ext(i2e, bus, gen, branch);
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        printmkt(baseMVA, bus, gen, branch, f, mkt.t, dispatch, success, et, fd, mpopt);
        fclose(fd);
    end
end
printmkt(baseMVA, bus, gen, branch, f, mkt.t, dispatch, success, et, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, baseMVA, bus, gen, branch, gencost);
end

if nargout
    mpc_out = struct(   'baseMVA', baseMVA, ...
                        'bus',     bus, ...
                        'gen',     gen, ...
                        'gencost', gencost, ...
                        'branch',  branch   );
end
