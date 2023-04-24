%% State Estimation for 30/50 Bus Test case

%% Initialization
clear
clc
close all
dbstop if error

%% Query user

usrInput = [];
usrStrtgy = [];
EXIT_FLAG = 0;

while isempty(usrInput)
    % Choose test case
    usrInput = input(['Choose the power system test case.\n',...
        'Enter 14 for the 14-bus test case, \n', ...        
        'or enter ieee30 for 30-bus ieee test case.\n', ...
        'To quit enter q: '],'s');
    
    switch lower(usrInput)
        case 'q'
            EXIT_FLAG = 1;
            break
        case '14'
            casename = 'case14';       
        case 'ieee30'
            casename = 'case_ieee30';
        otherwise
            disp('Unsupported test case!')
            usrInput = [];
    end
end

if (EXIT_FLAG ~= 1)
    % Choose state selection strategy
    fprintf(['It is important to provide the appropriate measurements\n', ...
        'such that the system is observable. We have defined several\n', ...
        'strategies for selecting observable states.\n', ...
        'For more information refer to\n', ...
        'https://github.com/aldalahmeh/matpower-extras/wiki\n'])
end

while isempty(usrStrtgy) && (EXIT_FLAG ~= 1)
    % 14 bus test case
    if strcmp(usrInput, '14')
        usrStrtgy = input(['Available observable states selection strategies are:\n', ...
            '1- Critical.\n',...
            '2- Davoudi.\n', ...
            '3- Bo.\n', ...
            '4- Kerdchuen.\n', ...
            'Enter the number of the selection: ']);
        
        switch usrStrtgy
            case 1
                Strtgy = 'critical';
            case 2
                Strtgy = 'Davoudi';
            case 3
                Strtgy = 'Bo';
            case 4
                Strtgy = 'Kerdchuen';
            otherwise
                disp('Unsupported strategy!')
                usrStrtgy = [];
        end
    
    % IEEE 30 bus test case
    elseif strcmp(usrInput, 'ieee30')
        usrStrtgy = input(['The ony available observable states selection strategies is:\n', ...
            '1- Critical.\n',...
            'Enter the number of the selection: ']);
        
        switch usrStrtgy
            case 1
                Strtgy = 'critical';
            otherwise
                disp('Unsupported strategy!')
                usrStrtgy = [];
        end
    else
        
    end
end

%% Simulation
if EXIT_FLAG ~= 1
    
    sigma.sigma_PD = 0.01;
    sigma.sigma_QD = 0.01;
    sigma.sigma_PG = 0.01;
    sigma.sigma_QG = 0.01;
    sigma.sigma_PF = 0.01;
    sigma.sigma_PT = 0.01;
    sigma.sigma_Va = 0.01;
    sigma.sigma_QF = 0.01;
    sigma.sigma_QT = 0.01;
    sigma.sigma_Vm = 0.01;
    
    %% Extract true data from test case
    [P, Q, V, Ind, baseMVA] = extrctMpcData(casename);
    
    %% Generate noisy measurements
    [meas, idx, sigma_mod] = genMsrmnt(casename, Strtgy, P, Q, V,...
        Ind, baseMVA, sigma);
    
    %% State estimation
    type_initialguess = 2; % flat start
    [baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] ...
        = run_se(casename, meas, idx, sigma_mod, type_initialguess);
    
end