function test_se_Ex_6_17
%TEST_SE  Test state estimation.
%   created by Rui Bo on 2007/11/12

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Rui Bo
%
%   This file is part of MATPOWER/mx-se.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mx-se/ for more info.

%%------------------------------------------------------
% using data in Problem 6.17 in book 'Computational
% Methods for Electric Power Systems' by Mariesa Crow
%%------------------------------------------------------
%% which measurements are available
idx.idx_zPD = [3]; % P3  -> Power demand at bus 3**addition to code**
idx.idx_zQD = [];  % Q3  -> Reactive Power demand **addition to code**
idx.idx_zPF = [2]; % P13 -> 2nd entry in the branch matrix **addition to code**
idx.idx_zPT = [];
idx.idx_zPG = [];
idx.idx_zVa = [];
idx.idx_zQF = [];
idx.idx_zQT = [1]; % Q21 -> 1st entry in the branch matrix
idx.idx_zQG = [2]; % Q2  -> Injected reactive power at bus 2
idx.idx_zVm = [3]; % V3  -> Voltage magnitude at bus 3

%% specify measurements
measure.PD = [-1.181]; % **addition to code**
measure.QD = []; % **addition to code**
measure.PF = [0.668];
measure.PT = [];
measure.PG = [];
measure.Va = [];
measure.QF = [];
measure.QT = [-0.082];
measure.QG = [-0.086];
measure.Vm = [0.975];

%% specify measurement variances
sigma.sigma_PD = [0.050] ; % **addition to code**
sigma.sigma_QD = [] ;      % **addition to code**
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
% if ~success
%     error('State Estimation input data are not complete or sufficient!');
% end

%% run state estimation
casename = 'case3bus_Ex6_17.m';
type_initialguess = 2; % flat start

fprintf('State estimation using the original doSE implemetination\n')
[baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] ...
    = run_se_err(casename, measure, idx, sigma, type_initialguess);

fprintf('\nState estimation using the modified doSE implemetination\n')
[baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] ...
    = run_se_mod(casename, measure, idx, sigma, type_initialguess);
%     = run_se_ext(casename, measure, idx, sigma, type_initialguess);
