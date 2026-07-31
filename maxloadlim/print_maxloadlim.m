function print_maxloadlim(mpc,results)
% PRINT_MAXLOADLIM(MPC,RESULTS) prints the results in RESULTS of the 
% maximum loadability problems defined from the base case MPC.

%   MATPOWER
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Camille Hamon
%
%   This file is part of MATPOWER/mx-maxloadlim.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mx-maxloadlim/ for more info.

define_constants;

% Print some global information about the parameters and results of the
% maximum loadability limit problem
mp_printf('\n');
mp_printf('=======================================\n');
mp_printf('      Maximum loadability problem\n');
mp_printf('=======================================\n');
mp_printf('The stability margin is %.2f MW\n',results.stab_marg*results.baseMVA);
mp_printf('The type of bifurcation is %s (%s).\n',...
    results.bif.full_name,results.bif.short_name);
if strcmp(results.bif.short_name,'LIB')
    mp_printf('Generator responsible for LIB: Gen %d connected at bus %d\n',...
        results.bif.gen_sll,results.gen(results.bif.gen_sll,GEN_BUS));
end
mp_printf('\n--------------------------------------------------\n');
mp_printf('   Bus nb    Direction     Load at MLL    Voltages \n');
mp_printf('   ------    ---------     -----------    --------\n');
for i = 1:size(results.bus,1)
    mp_printf('   %4d      %.2f         %8.2f         %4.3f\n',...
        results.bus(i,BUS_I),results.dir_mll(i),results.bus(i,PD),results.bus(i,VM));
end

% Print some global information about the generators and their limits
% First we extract the original voltage set points of the generators
gen_v_set = mpc.gen(:,VG);

mp_printf('\n');
mp_printf('=============================================================\n');
mp_printf('Reactive power production and voltages at the generators\n');
mp_printf('=============================================================\n');
mp_printf(' Bus nb       Qgen        Qmin        Qmax       Vm      Vref\n');
mp_printf('--------    --------    --------    --------    -----    -----\n');
for i = 1:size(results.gen,1)
    mp_printf(' %4d       %8.2f    %8.2f    %8.2f    %4.3f    %4.3f',...
        results.gen(i,GEN_BUS),results.gen(i,QG),results.gen(i,QMIN),results.gen(i,QMAX),...
        results.bus(results.gen(i,GEN_BUS),VM),gen_v_set(i));
    if results.bus(results.gen(i,GEN_BUS),BUS_TYPE) == REF
        mp_printf('    REF');
    end
    if strcmp(results.bif.short_name,'LIB') && any(results.bif.gen_sll == i)
        mp_printf('    LIB');
    end
    mp_printf('\n');
end
