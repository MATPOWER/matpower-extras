function [MVAbase, bus, gen, gencost, branch, f, dispatch, success, et] = ...
				runmkt(casename, q, p, mkt, max_p, u0, t, mpopt, fname, solvedcase)
%RUNMKT  Runs smart market for PowerWeb, computing a new generation
%        schedule from a set of offers.
%
%   [baseMVA, bus, gen, gencost, branch, f, dispatch, success, et] = ...
%           runmkt(casename, q, p, mkt, max_p, u0, t, mpopt, fname, solvedcase)
%
%   Computes the new generation and price schedules based on the offers
%   submitted, where offers are specified by q and p, mkt tells it what
%   type of market to use, max_p is the reservation price, u0 is a vector
%   containing the commitment status of each generator from the previous
%   period (for computing startup/shutdown costs), t is the time duration
%   of the dispatch period in hours, and mpopt is a MATPOWER options vector
%   (see 'help mpoption' for details). Uses default options if mpopt is not
%   given. The rows in q and p correspond to the rows in gen and gencost,
%   and each column corresponds to another data point or block in the
%   marginal offer. The market codes are defined as the sum of the
%   following numbers:
%       10000               - all markets
%        1000 * discrete    - discrete = 1 for block offers, 0 for continuous
%                             linear offers
%         100 * adjust4loc  - adjust4loc = 0 to ignore network,
%							  1 to compute locational adjustments via AC OPF,
%                             2 to compute them via DC OPF
%          10 * which_price - which_price = 1 for last accepted offer, 2 for
%                             first rejected offer
%   If p or q are empty or not given, they are created from the generator
%   cost function. The default is market code is 11110, for the standard
%   LAO EMPIRE auction. The default max_p is 500, the default u0 is all ones
%   (assume everything was running) and the default duration t is 1 hour.
%   The results may optionally be printed to a file (appended
%   if the file exists) whose name is given in fname (in addition to
%   printing to STDOUT). Optionally returns the final values of baseMVA,
%   bus, gen, gencost, branch, f, dispatch, success, and et. If a name is
%   given in solvedcase, the solved case will be written to a case file
%   in MATPOWER format with the specified name with a '.m' extension added.

%   MATPOWER Version 2.5b3
%   by Ray Zimmerman, PSERC Cornell    9/21/99
%   Copyright (c) 1996-1999 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%-----  initialize  -----
%% default arguments
if nargin < 10
	solvedcase = '';							%% don't save solved case
	if nargin < 9
		fname = '';								%% don't print results to a file
		if nargin < 8
			mpopt = mpoption;					%% use default options
			if nargin < 7
				t = [];							%% use default dispatch period duration (hours)
				if nargin < 6
					u0 = [];					%% use default for previous gen commitment
					if nargin < 5
						max_p = 500;			%% use default reservation price
						if nargin < 4
							mkt = [];			%% use default market
							if nargin < 3
								q = []; p = [];	%% p & q not defined (will use gencost)
								if nargin < 1
									casename = 'case';	%% default data file is 'case.m'
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
	mkt = 11110;				%% default market is LAO EMPIRE market
end
if isempty(max_p)
	max_p = 500;				%% default reservation price is 500
end
if isempty(t)
	t = 1;						%% default dispatch duration in hours
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch, area, gencost] = loadcase(casename);
[i2e, bus, gen, branch, area] = ext2int(bus, gen, branch, area);

%% finish assigning default arguments
if isempty(u0)
	u0 = ones(size(gen, 1), 1);	%% default for previous gen commitment, all on
end

%% if q and p not defined, use gencost
if isempty(q) | isempty(p)
	discrete = fix((mkt-10000)/1000);
	[q, p] = case2off(gen, gencost, discrete);
end

%% start the clock
t0 = clock;

%% run the market
[bus, gen, branch, f, dispatch, success] = ...
		smartmkt(baseMVA, bus, gen, gencost, branch, area, q, p, mkt, max_p, u0, t, mpopt);

%% compute elapsed time
et = etime(clock, t0);

%% convert back to original bus numbering & print results
[bus, gen, branch, area] = int2ext(i2e, bus, gen, branch, area);
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
	savecase(solvedcase, baseMVA, bus, gen, branch, area, gencost);
end

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
