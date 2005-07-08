function [gen, gencost] = off2case(gen, gencost, q, p, max_p)
%OFF2CASE  Updates case variables gen & gencost from quantity & price offers.
%   [gen, gencost] = off2case(gen, gencost, q, p, max_p) updates gen &
%   gencost variables based on quantity and price offers and the market type
%   specified. Updates PMIN, PMAX, GEN_STATUS and all cost info except STARTUP
%   and SHUTDOWN. Any quantity offered above max_p will be ignored.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
    RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%% save old gencost
oldgencost = gencost;
[ng, np]    = size(q);

%% do conversion
pmin = zeros(ng, 1);
pmax = zeros(ng, 1);
qmin = zeros(ng, 1);
qmax = zeros(ng, 1);

gencost             = zeros(ng, COST + 2*np - 1);
gencost(:, MODEL)   = PW_LINEAR * ones(ng, 1);
for i = 1:ng
    qq = q(i, :)';          %% column vector of quantity offers for gen i
    pp = p(i, :)';          %% column vector of price offers for gen i

    if isload(gen(i, :))
        %% strip zero quantities, and flip bids to turn them into fake offers
        valid = find(qq);
        n = length(valid);
        qq = qq(valid(n:-1:1)); %% column vector of quantities of valid offers for gen i
        pp = pp(valid(n:-1:1)); %% column vector of prices of valid offers for gen i
    else
        %% strip zero quantities and prices over max_p
        valid = find(qq & pp <= max_p);
        qq = qq(valid);         %% column vector of quantities of valid offers for gen i
        pp = pp(valid);         %% column vector of prices of valid offers for gen i
    end
    n = length(qq) + 1;

    %% form piece-wise linear total cost function, set Pmin & Pmax
    if n > 1        %% otherwise, leave all cost info zero (specifically N)
        %% set Pmin and Pmax
        if isload(gen(i, :))    %% it's a load
            pmin(i) = -sum(qq);
            %% scale Q limits to keep original power factor specification
            qmin(i) = gen(i, QMIN) * pmin(i) / gen(i, PMIN);
            qmax(i) = gen(i, QMAX) * pmin(i) / gen(i, PMIN);
%           pmax(i) = gen(i, PMAX);     %% no change
            xx = [0; cumsum(qq)]' + pmin(i);
            yy = [0; cumsum(pp .* qq)]';
            yy = yy - yy(end);
        else                    %% it's a generator
%           pmin(i) = qq(1);
%           pmin(i) = gen(i, PMIN);     %% no change
            pmax(i) = sum(qq);
            xx = [0; cumsum(qq)]';
            yy = [0; cumsum(pp .* qq)]';
        end
    
        %% form piece-wise linear total cost function
        gencost(i, N) = n;
        gencost(i,      COST:2:( COST + 2*n - 1 )) = xx;
        gencost(i,  (COST+1):2:( COST + 2*n     )) = yy;
    end
end

%% copy back startup and shutdown costs
gencost(:, [STARTUP SHUTDOWN]) = oldgencost(:, [STARTUP SHUTDOWN]);

%% set PMIN, PMAX, GEN_STATUS
off = find(gencost(:, N) == 0);                 %% find gens with no valid offers
gen(off, GEN_STATUS) = zeros(length(off), 1);   %% turn them off
G = find( gen(:, GEN_STATUS) > 0 & ~isload(gen) );  %% in-service gens
L = find( gen(:, GEN_STATUS) > 0 &  isload(gen) );  %% in-service disp loads

%% update limits
gen(G, PMAX) = pmax(G);
gen(L, PMIN) = pmin(L);
gen(L, QMIN) = qmin(L);
gen(L, QMAX) = qmax(L);

return;
