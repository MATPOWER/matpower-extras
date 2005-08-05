function [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
                runmkt(casename, q, p, mkt, max_p, u0, t, mpopt, fname, solvedcase)
%RUNMKT  Runs smart market for PowerWeb, computing a new generation
%        schedule from a set of offers and bids.
%
%   [baseMVA, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
%           runmkt(casename, q, p, mkt, max_p, u0, t, mpopt, fname, solvedcase)
%
%   Computes the new generation and price schedules based on the offers
%   submitted, where offers are specified by q and p, mkt tells it what
%   type of market to use, max_p is the price cap, u0 is a vector
%   containing the commitment status of each generator from the previous
%   period (for computing startup/shutdown costs), t is the time duration
%   of the dispatch period in hours, and mpopt is a MATPOWER options vector
%   (see 'help mpoption' for details). Uses default options if mpopt is not
%   given. The rows in q and p correspond to the rows in gen and gencost,
%   and each column corresponds to another block in the marginal offer or
%   bid. The market codes are defined as the sum of the
%   following numbers:
%        1000                - all markets
%         100 * adjust4loc   - adjust4loc = 0 to ignore network,
%                              1 to compute locational adjustments via AC OPF,
%                              2 to compute them via DC OPF
%          10 * auction_type - where the values for auction_type are as follows:
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
%   If p or q are empty or not given, they are created from the generator
%   cost function. The default market code is 1150, where the marginal
%   block (offer or bid) sets the price. The default max_p is 500, the
%   default u0 is all ones (assume everything was running) and the default
%   duration t is 1 hour. The results may optionally be printed to a file
%   (appended if the file exists) whose name is given in fname (in addition
%   to printing to STDOUT). Optionally returns the final values of baseMVA,
%   cq, cp, bus, gen, gencost, branch, f, dispatch, success, and et. If a
%   name is given in solvedcase, the solved case will be written to a case file
%   in MATPOWER format with the specified name with a '.m' extension added.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialize  -----
%% default arguments
if nargin < 10
    solvedcase = '';                        %% don't save solved case
    if nargin < 9
        fname = '';                         %% don't print results to a file
        if nargin < 8
            mpopt = mpoption;               %% use default options
            if nargin < 7
                t = [];                     %% use default dispatch period duration (hours)
                if nargin < 6
                    u0 = [];                %% use default for previous gen commitment
                    if nargin < 5
                        max_p = 500;        %% use default price cap
                        if nargin < 4
                            mkt = [];       %% use default market
                            if nargin < 3
                                q = []; p = []; %% p & q not defined (use gencost)
                                if nargin < 1
                                    casename = 'case9'; %% default data file is 'case9.m'
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if isempty(mkt)
    mkt = 1150;             %% default market is LAO EMPIRE market
end
if isempty(max_p)
    max_p = 500;            %% default price cap is 500
end
if isempty(t)
    t = 1;                  %% default dispatch duration in hours
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
    RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch, areas, gencost] = loadcase(casename);
[i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);

%% finish assigning default arguments
if isempty(u0)
    u0 = ones(size(gen, 1), 1); %% default for previous gen commitment, all on
end

%% if q and p not defined, use gencost
if isempty(q) | isempty(p)
    [q, p] = case2off(gen, gencost);
end

%% start the clock
t0 = clock;

%% run the market
[cq, cp, bus, gen, branch, f, dispatch, success] = ...
        smartmkt(baseMVA, bus, gen, gencost, branch, areas, q, p, mkt, max_p, u0, t, mpopt);

%% compute elapsed time
et = etime(clock, t0);

%% convert back to original bus numbering & print results
[bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        printmkt(baseMVA, bus, gen, branch, f, t, dispatch, success, et, fd, mpopt);
        fclose(fd);
    end
end
printmkt(baseMVA, bus, gen, branch, f, t, dispatch, success, et, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, baseMVA, bus, gen, branch, areas, gencost);
end

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
