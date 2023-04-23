function test_se_case3bus_Ex
%TEST_SE_CASE3BUS_EX  Test state estimation for 3 bus system based on the
%test case file "case2bus_Ex.m".
%   created by Rui Bo on 2007/11/12
%   modified by Sami Aldalahmeh 2023/15/2
%
%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Rui Bo
%
%   This file is part of MATPOWER/mx-se.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mx-se/ for more info.
%
%%------------------------------------------------------
%% which measurements are available
idx.idx_zPD = [3]; % P3  -> Power demand at bus 3
idx.idx_zQD = [];  
idx.idx_zPF = [2]; % P13 -> 2nd entry in the branch matrix
idx.idx_zPT = [];
idx.idx_zPG = [];
idx.idx_zVa = [];
idx.idx_zQF = [];
idx.idx_zQT = [1]; % Q21 -> 1st entry in the branch matrix
idx.idx_zQG = [2]; % Q2  -> Injected reactive power at bus 2
idx.idx_zVm = [3]; % V3  -> Voltage magnitude at bus 3

%% specify measurements
measure.PD = [-1.181]; % P3  -> Power demand at bus 3
measure.QD = []; 
measure.PF = [0.668];  % P13 -> 2nd entry in the branch matrix
measure.PT = [];
measure.PG = [];
measure.Va = [];
measure.QF = [];
measure.QT = [-0.082]; % Q21 -> 1st entry in the branch matrix
measure.QG = [-0.086]; % Q2  -> Injected reactive power at bus 2
measure.Vm = [0.975];  % V3  -> Voltage magnitude at bus 3

%% specify measurement variances
sigma.sigma_PD = [0.050] ;
sigma.sigma_QD = [] ;     
sigma.sigma_PF = [0.050] ;
sigma.sigma_PT = [];
sigma.sigma_PG = [];
sigma.sigma_Va = [];
sigma.sigma_QF = [];
sigma.sigma_QT = [0.075];
sigma.sigma_QG = [0.075];
sigma.sigma_Vm = 0.01;

%% check input data integrity
nbus = 3;
[success, measure, idx, sigma] = checkDataIntegrity(measure, idx, sigma, nbus);

%% run state estimation
casename = 'case3bus_Ex.m';
type_initialguess = 2; % flat start
[baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] ...
    = run_se(casename, measure, idx, sigma, type_initialguess);
