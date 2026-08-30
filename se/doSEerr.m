function [V, converged, iterNum, z, z_est, error_sqrsum] = doSEext(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V0, ref, pv, pq, measure, idx, sigma)
%DOSE  Do state estimation (extended version).
%   created by Rui Bo on 2007/11/12

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Rui Bo
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER/mx-se.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mx-se/ for more info.

%% Debugging
% load se_org
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
    GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% options
tol     = 1e-5; % mpopt.pf.tol;
max_it  = 100;  % mpopt.pf.nr.max_it;
verbose = 0;

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

nb = size(Ybus, 1);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses

%% get non reference buses
nonref = [pv;pq];
% For PV buses the required state is the voltage angle.
% Whereas for PQ buses both the voltage magnitude and angle are required.
nonref_a = [pv;pq]; % angle buses indices
nonref_m = [pq]; % magnitude buses indices

% !!!!!!
% nonref_a = nonref;
% nonref_m = nonref;

% Generators buses
gbus = gen(:, GEN_BUS);
gBusInd = 1:length(gbus);

% Create indencies for state vector.
% The first elements are the voltages' angles then come the magnitudes. 
% Note that for the PQ buses, the states are the angles and the magnitudes,
% hence double the number of buses.
Va_ind = 1:(length(pv)+length(pq));
Vm_ind = length(Va_ind)+1:length(pv)+2*length(pq);

%% form measurement vector 'z'. NOTE: all are p.u. values
z = [
    measure.PD  % **addition to code**
    measure.PF
    measure.PT
    measure.PG
    measure.Va
    measure.QD % **addition to code**
    measure.QF
    measure.QT
    measure.QG
    measure.Vm
    ];

%% form measurement index vectors
idx_zPD = idx.idx_zPD; % **addition to code**
idx_zQD = idx.idx_zQD; % **addition to code**
idx_zPF = idx.idx_zPF;
idx_zPT = idx.idx_zPT;
idx_zPG = idx.idx_zPG;
idx_zVa = idx.idx_zVa;
idx_zQF = idx.idx_zQF;
idx_zQT = idx.idx_zQT;
idx_zQG = idx.idx_zQG;
idx_zVm = idx.idx_zVm;

% idx_zPG = gBusInd(ismember(gbus,idx_zPG));
% idx_zQG = gBusInd(ismember(gbus,idx_zQG));

%% get R inverse matrix
% **addition to code**

if length(sigma.sigma_PF) > 1
    
    sigma_vector = [
        sigma.sigma_PD % **addition to code**
        sigma.sigma_PF
        sigma.sigma_PT
        sigma.sigma_PG
        sigma.sigma_Va
        sigma.sigma_QD % **addition to code**
        sigma.sigma_QF
        sigma.sigma_QT
        sigma.sigma_QG
        sigma.sigma_Vm
        ];
else
    % % **remove from code**
    sigma_vector = [
        sigma.sigma_PD*ones(size(idx_zPD, 1), 1) % **addition to code**
        sigma.sigma_PF*ones(size(idx_zPF, 1), 1)
        sigma.sigma_PT*ones(size(idx_zPT, 1), 1)
        sigma.sigma_PG*ones(size(idx_zPG, 1), 1)
        sigma.sigma_Va*ones(size(idx_zVa, 1), 1)
        sigma.sigma_QD*ones(size(idx_zQD, 1), 1) % **addition to code**
        sigma.sigma_QF*ones(size(idx_zQF, 1), 1)
        sigma.sigma_QT*ones(size(idx_zQT, 1), 1)
        sigma.sigma_QG*ones(size(idx_zQG, 1), 1)
        sigma.sigma_Vm*ones(size(idx_zVm, 1), 1)
        ]; % NOTE: zero-valued elements of simga are skipped
end
sigma_square = sigma_vector.^2;
% R_inv = diag(1./sigma_square);
R_inv = inv(diag(1./sigma_square));

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    
    %% --- compute estimated measurement ---
    % **remove from code**
    Sfe = V(f) .* conj(Yf * V);
    Ste = V(t) .* conj(Yt * V); % 
    
    %% Power injection at buses **addition to code**
%     [Sfe, Ste] = cmptSmat(V, Ybus, branch);
 
    
    %% Power demand at buses 
    Sbus = V .* conj(Ybus * V);
    
    %% compute net injection at generator buses
%     gbus = gen(:, GEN_BUS);
    Sgbus = V(gbus) .* conj(Ybus(gbus, :) * V);
    Sgen = Sgbus * baseMVA + (bus(gbus, PD) + 1j*bus(gbus, QD));    %% inj S + local Sd
    Sgen = Sgen/baseMVA;
    z_est = [ % NOTE: all are p.u. values
        real(Sbus(idx_zPD)); % Demand real power **addition to code**
        real(Sfe(idx_zPF));
        real(Ste(idx_zPT));
        real(Sgen(idx_zPG));
        angle(V(idx_zVa));
        imag(Sbus(idx_zQD)); % Demand reactive power **addition to code**
        imag(Sfe(idx_zQF));
        imag(Ste(idx_zQT));
        imag(Sgen(idx_zQG));
        abs(V(idx_zVm));
    ];

    %% --- get H matrix ---
    [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
    [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V);
%     genbus_row = findBusRowByIdx(bus, gbus);
    genbus_row = gbus;  %% rdz, this should be fine if using internal bus numbering

    %% get sub-matrix of H relating to power demand **addition to code**
    dPD_dVa = real(dSbus_dVa);
    dPD_dVm = real(dSbus_dVm);
    
    dQD_dVa = imag(dSbus_dVa);
    dQD_dVm = imag(dSbus_dVm);
    
    %% get sub-matrix of H relating to line flow    
    dPF_dVa = real(dSf_dVa); % from end
    dQF_dVa = imag(dSf_dVa);   
    dPF_dVm = real(dSf_dVm);
    dQF_dVm = imag(dSf_dVm);
    dPT_dVa = real(dSt_dVa);% to end
    dQT_dVa = imag(dSt_dVa);   
    dPT_dVm = real(dSt_dVm);
    dQT_dVm = imag(dSt_dVm);   
    
    %% get sub-matrix of H relating to generator output
    dPG_dVa = real(dSbus_dVa(genbus_row, :));
    dQG_dVa = imag(dSbus_dVa(genbus_row, :));
    dPG_dVm = real(dSbus_dVm(genbus_row, :));
    dQG_dVm = imag(dSbus_dVm(genbus_row, :));
    
    %% get sub-matrix of H relating to voltage angle
    dVa_dVa = eye(nb);
    dVa_dVm = zeros(nb, nb);
    
    %% get sub-matrix of H relating to voltage magnitude
    dVm_dVa = zeros(nb, nb);
    dVm_dVm = eye(nb);
    % **addition to code**
      H = [
        dPD_dVa(idx_zPD, nonref_a)   dPD_dVm(idx_zPD, nonref_m);% PD deriv. c        
        dPF_dVa(idx_zPF, nonref_a)   dPF_dVm(idx_zPF, nonref_m);
        dPT_dVa(idx_zPT, nonref_a)   dPT_dVm(idx_zPT, nonref_m);
        dPG_dVa(idx_zPG, nonref_a)   dPG_dVm(idx_zPG, nonref_m);
        dVa_dVa(idx_zVa, nonref_a)   dVa_dVm(idx_zVa, nonref_m);
        dQD_dVa(idx_zQD, nonref_a)   dQD_dVm(idx_zQD, nonref_m);% QD deriv. c
        dQF_dVa(idx_zQF, nonref_a)   dQF_dVm(idx_zQF, nonref_m);
        dQT_dVa(idx_zQT, nonref_a)   dQT_dVm(idx_zQT, nonref_m);
        dQG_dVa(idx_zQG, nonref_a)   dQG_dVm(idx_zQG, nonref_m);
        dVm_dVa(idx_zVm, nonref_a)   dVm_dVm(idx_zVm, nonref_m);
        ];
%     H = [
%         dPD_dVa(idx_zPD, nonref)   dPD_dVm(idx_zPD, nonref);% PD deriv. **addition to code**
%         dPF_dVa(idx_zPF, nonref)   dPF_dVm(idx_zPF, nonref);
%         dPT_dVa(idx_zPT, nonref)   dPT_dVm(idx_zPT, nonref);
%         dPG_dVa(idx_zPG, nonref)   dPG_dVm(idx_zPG, nonref);
%         dVa_dVa(idx_zVa, nonref)   dVa_dVm(idx_zVa, nonref);
%         dQD_dVa(idx_zQD, nonref)   dQD_dVm(idx_zQD, nonref);
%         dQF_dVa(idx_zQF, nonref)   dQF_dVm(idx_zQF, nonref);
%         dQT_dVa(idx_zQT, nonref)   dQT_dVm(idx_zQT, nonref);
%         dQG_dVa(idx_zQG, nonref)   dQG_dVm(idx_zQG, nonref);
%         dVm_dVa(idx_zVm, nonref)   dVm_dVm(idx_zVm, nonref);
%         ];
    
    %% compute update step
    J = H'*R_inv*H; % gain matrix
    F = H'*R_inv*(z-z_est); % evalute F(x)
    if ~isobservable(H, pv, pq)
        error('doSE: system is not observable');
    end
    dx = (J \ F);
    %% Debugging
    DEBUG_FLAG = true;
    
    if DEBUG_FLAG
       % Rearange H matrix as ex. 6.17 in book
       Hd = H;
       Hd([1 end],:)=Hd([end 1],:);
       Hd([end end-1],:)=Hd([end-1 end],:);
       fprintf('H matrix in interation #%d\n',i)
       fprintf('H = \n')
       disp(Hd)
       
    end
    
    
    %% check for convergence
    normF = norm(F, inf);
    if verbose > 1
        fprintf('\niteration [%3d]\t\tnorm of mismatch: %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
    end
    
    %% update voltage 
%     Va(nonref) = Va(nonref) + dx(1:size(nonref, 1));
    Va(nonref_a) = Va(nonref_a) + dx(Va_ind); % **addition code**
%     Vm(nonref) = Vm(nonref) + dx(size(nonref, 1)+1:2*size(nonref, 1));
    Vm(nonref_m) = Vm(nonref_m) + dx(Vm_ind);% **addition code**
    
%     y = [Va(nonref_a);Vm(nonref_m)]; % **addition code**
    
    V = Vm .* exp(1j * Va); % NOTE: angle is in radians in pf solver, but in degree in case data
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    
%     %% debugging
%     re = sum([abs(abs(V_cell{i}) - Vm),abs(angle(V_cell{i}) - Va)]);
%     disp(['iter: ', num2str(i), ' Error: ' num2str(re) ])
%     disp('')
end

iterNum = i;

%% get weighted sum of squared errors
error_sqrsum = sum((z - z_est).^2./sigma_square);

%% New function added to code
function BusInd = getInd(ind, pv, pq)
% Return the bus index from the H matrix column index.

Lpv = length(pv);
Lpq = length(pq);

for i = 1:length(ind)
    if ind(i) <= Lpv
        BusInd(i) = pv(ind(i));
    elseif (ind(i) > Lpv) && (ind(i) <= Lpv+Lpq)
        BusInd(i) = pq(ind(i) - Lpv);
    else
        BusInd(i) = pq(ind(i) - Lpv - Lpq);
    end
end

