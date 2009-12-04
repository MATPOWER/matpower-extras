function printmkt(baseMVA, bus, gen, branch, f, t, dispatch, success, et, fd, mpopt)
%PRINTMKT   Prints results of ISO computation.
%   printmkt(baseMVA, bus, gen, branch, f, t, dispatch, success, et, fd, mpopt)
%   prints results of ISO computation to fd (a file descriptor which defaults
%   to STDOUT). mpopt is a MATPOWER options vector (see 'help mpoption' for
%   details). Uses default options if this parameter is not given. The objective
%   function value is given in f, the duration of the dispatch period (in hours) in
%   t, and the elapsed time in et.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% default arguments
if nargin < 11
    mpopt = mpoption;   %% use default options
    if nargin < 10
        fd = 1;         %% print to stdio by default
    end
end

%% options
OUT_ALL         = mpopt(32);
OUT_RAW         = mpopt(43);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[QUANTITY, PRICE, FCOST, VCOST, SCOST, PENALTY] = idx_disp;

%% parameters
ng = size(gen, 1);

%%----- print the stuff -----
pay = dispatch(:, PRICE) .* dispatch(:, QUANTITY) * t;
cost = dispatch(:, FCOST) + dispatch(:, VCOST) + dispatch(:, SCOST) + dispatch(:, PENALTY);
if OUT_ALL
    %% dispatch data
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Market Summary                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nDispatch period duration: %.2f hours', t);
    fprintf(fd, '\nGen  Bus     Pg      Price    Revenue   Fix+Var   Strt/Stp   Total    Earnings');
    fprintf(fd, '\n #    #     (MW)    ($/MWh)     ($)     Cost ($)  Cost ($)  Cost ($)     ($)  ');
    fprintf(fd, '\n---  ---  --------  --------  --------  --------  --------  --------  --------');
    for i = 1:size(gen, 1)
        if gen(i, PG)
            fprintf(fd, '\n%3d%5d%9.2f%10.3f%10.2f%10.2f%10.2f%10.2f%10.2f', ...
                i, gen(i, GEN_BUS), dispatch(i, QUANTITY), dispatch(i, PRICE), pay(i), ...
                dispatch(i, FCOST) + dispatch(i, VCOST), ...
                dispatch(i, SCOST), cost(i), pay(i) - cost(i));
        else
            if dispatch(i, SCOST) || dispatch(i, PENALTY)
                fprintf(fd, '\n%3d%5d      -  %10.3f       -         -  %10.2f%10.2f%10.2f', ...
                    i, gen(i, GEN_BUS), dispatch(i, PRICE), dispatch(i, SCOST), ...
                    cost(i), pay(i) - cost(i));
            else
                fprintf(fd, '\n%3d%5d      -  %10.3f       -         -         -         -         -', ...
                    i, gen(i, GEN_BUS), dispatch(i, PRICE));
            end
        end
        if dispatch(i, PENALTY)
            fprintf(fd, '%10.2f penalty (included in total cost)', dispatch(i, PENALTY));
        end
    end
    fprintf(fd, '\n          --------            --------  --------  --------  --------  --------');
    fprintf(fd, '\nTotal:  %9.2f          %10.2f%10.2f%10.2f%10.2f%10.2f', ...
        sum(dispatch(:, QUANTITY)), sum(pay), sum(dispatch(:, FCOST)) + sum(dispatch(:, VCOST)), ...
        sum(dispatch(:, SCOST)), sum(cost), sum(pay-cost));
    if sum(dispatch(:, PENALTY))
        fprintf(fd, '%10.2f penalty (included in total cost)', sum(dispatch(:, PENALTY)));
    end
    fprintf(fd, '\n');
end

%% print raw data for Perl database interface
if OUT_RAW
    fprintf(fd, '----------  raw PW::Dispatch data below  ----------\n');
    fprintf(fd, 'dispatch\n');
    fprintf(fd, '%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n', ...
                [(1:ng)' dispatch(:, [QUANTITY, PRICE, FCOST, VCOST, SCOST, PENALTY]) pay-cost]');
    fprintf(fd, '----------  raw PW::Dispatch data above  ----------\n');
end

%% print remaining opf output
printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt);
