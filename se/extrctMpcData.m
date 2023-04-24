function [P, Q, V, Ind, baseMVA] = extrctMpcData(PowerTestCase, varargin)
% Extract measurements from the IEEE test case "PowerTestCase".
% Input(s):
% ------------------------------------------------------------------------
% * PowerTestCase: string; name of the test case
%               OR struct; MATPOWER load flow output structure.
%
% Output(s):
% ------------------------------------------------------------------------
%  P   : struct; containing real power values.
%   .PG: array; real generated power.
%   .PD: array; real power load demand.
%   .PS: array; real power shunt consumption.
%   .PF: array; real power from-to bus.
%   .PT: array; real power to-from bus.
%
%  Q   : struct; containing reactive power values.
%   .QG: array; reactive generated power.
%   .QD: array; reactive power load demand.
%   .QS: array; reactive power shunt consumption.
%   .QF: array; reactive power from-to bus.
%   .QT: array; reactive power to-from bus.
%
% Ind  : struct; containing indicies of the buses for above values.
%   .PG: array; real generated power index.
%   .PD: array; real power load demand index.
%   .PS: array; real power shunt consumption index.
%   .PF: array; real power from-to bus index.
%   .PT: array; real power to-from bus index.
%   .QG: array; reactive generated power index.
%   .QD: array; reactive power load demand index.
%   .QS: array; reactive power shunt consumption index.
%   .QF: array; reactive power from-to bus index.
%   .QT: array; reactive power to-from bus index.
%
% * V: struct;
%   V.Vm: array of voltage magnitudes.
%   V.Va: array of voltage angles.
%
% NOTE: The all the values are ordered as per the mpc structure.
%
define_constants;
if ischar(PowerTestCase)
    %% Load & run test case
    % Load test case
    mpc = loadcase(PowerTestCase);
    baseMVA = mpc.baseMVA;
    
    
    % Suppress output
    mpopt = mpoption('verbose', 0, 'out.all', 0);
    
    % Run power flow
    results_pf = runpf(mpc, mpopt);
    
elseif isstruct(PowerTestCase)
    if nargin == 2
        mpc = varargin{1};
        results_pf = PowerTestCase;
    else
        error('Incorrect input argumets!')
    end
    
end

%% Extract Power 
%% Generated and load power
[P, Q, Ind] = extrctPwr(mpc);

%% Power flows
PFval = results_pf.branch(:,PF);
PTval = results_pf.branch(:,PT);

QFval = results_pf.branch(:,QF);
QTval = results_pf.branch(:,QT);

%% Extract voltage
Vm = results_pf.bus(:,VM); % magnitudes
Va = results_pf.bus(:,VA); % angles

%% Add to structure
% Real power
P.PF = PFval; % From-to
P.PT = PTval; % To-from

% Reactive power
Q.QF = QFval;
Q.QT = QTval;

% Voltage
V.Vm = Vm;
V.Va = Va;

%% Indicies
Ind.PF = 1:length(mpc.branch(:,1));
Ind.QF = Ind.PF;

Ind.PT = 1:length(mpc.branch(:,2));
Ind.QT = Ind.PT;

Ind.Vm = 1:length(Vm);
Ind.Va = 1:length(Va);

function [P, Q, Ind] = extrctPwr(mpc)
% Extract generated and load power from the test case structure.
% Input(s):
% * mpc: struct; test case structure.
%
% Output(s):
% * P: struct; containing real power values.
%   .PG: array; real generated power.
%   .PD: array; real power load demand.
%   .PS: array; real power shunt consumption. (removed)
%
% * Q: struct; containing reactive power values.
%   .QG: array; reactive generated power.
%   .QD: array; reactive power load demand.
%   .QS: array; reactive power shunt consumption.  (removed)
%
% * Ind: struct; containing indicies of the buses for above values.
%   .PG: array; real generated power indicies.
%   .PD: array; real power load demand indicies.
%   .PS: array; real power shunt consumption indicies.  (removed)
%   .QG: array; reactive generated power indicies.
%   .QD: array; reactive power load demand indicies.
%   .QS: array; reactive power shunt consumption indicies. (removed)
%   .UL: array; unloaded PQ buses indicies.
%
% NOTE: The bus indicies are per the mpc test case.

%% Define column indicies
PG_IND = 2; QG_IND = 3; % Generated power
PD_IND = 3; QD_IND = 4; % Power demand
GS_IND = 5; BS_IND = 6; % Power demand

%% Extract PV & PQ bus indicies
[pvBusInd, pqBusInd, refBusInd] = extrctBus(mpc);
nStates = numel(pvBusInd) + 2*numel(pqBusInd);

%% Generated power
nGen = size(mpc.gen,1);

% Index in the bus matrix
% gInd = ismember(mpc.gen(:,1), pvBusInd );
genBusInd = [refBusInd;pvBusInd];

% Index in the generator matrix
% pGenInd = 1:nGen; %pGenInd = pGenInd(gInd);
pGenInd = genBusInd;
qGenInd = pGenInd;

% Subtracted loads at PV buses
PG = mpc.gen(:,PG_IND) - mpc.bus(genBusInd,PD_IND);
QG = mpc.gen(:,QG_IND) - mpc.bus(genBusInd,QD_IND);

% PG = mpc.gen(pGenInd,PG_IND) - mpc.bus(pvBusInd,PD_IND);
% QG = mpc.gen(qGenInd,QG_IND) - mpc.bus(pvBusInd,QD_IND);

%% Demand power excluding shunt elements consumptions 

% Load bus indicies and power
[PD, pLdInd] = extrctColBus(PD_IND, mpc.bus); % real
PD(ismember(pLdInd, pvBusInd) ) = [];
pLdInd( ismember(pLdInd, pvBusInd) ) = []; % remove PV buses from demand

[QD, qLdInd] = extrctColBus(QD_IND, mpc.bus); % reactive
QD(ismember(qLdInd, pvBusInd) ) = [];
qLdInd( ismember(qLdInd, pvBusInd) ) = []; % remove PV buses from demand

% Shunt bus indicies and power
[PS, gShntInd] = extrctColBus(GS_IND, mpc.bus); % real 
[QS, bShntInd] = extrctColBus(BS_IND, mpc.bus); % reactive

% Collect power values in a structure
P.PG = PG; Q.QG = QG;

% Load power include shunt elements
P.PD = PD; 
P.PD = addShntPwr(P.PD, PS, pLdInd, gShntInd);

Q.QD = QD; 
Q.QD = addShntPwr(Q.QD, QS, qLdInd, bShntInd);
          
% Unloaded PQ buses indicies
unLdPQInd = pqBusInd;
unLdPQInd( ismember(pqBusInd, pLdInd) ) = [];

% Collect power indicies in a structure
Ind.PD = pLdInd;
Ind.QD = qLdInd;
Ind.PG = pGenInd;  Ind.QG = qGenInd;
Ind.UL = unLdPQInd;

function [pvBusInd, pqBusInd, refBusInd] = extrctBus(mpc)
% Extract PV and PQ bus indicies from the test case structure "mpc" as
% provided by MATPOWER.
% Input(s):
% * mpc: struct; IEEE test case structure.
%
% Output(s):
% * pvBusInd: array; PV bus indicies.
% * pqBusInd: array; PQ bus indicies.

allBusInd = mpc.bus(:,1);  % Bus index
allBusType = mpc.bus(:,2); % Bus type 
% PV & PQ busses indicies 
pvBusInd = allBusInd(allBusType == 2); % 2 is PV bus type as in IEEE test case
pqBusInd = allBusInd(allBusType == 1); % 1 is PQ bus type as in IEEE test case
refBusInd = allBusInd(allBusType == 3); % 3 is ref. bus type as in IEEE test case


function [Val, ind] = extrctColBus(colInd, MAT)
% Extract value from the bus matrix
% Input(s):
% * colInd: scalar; column index.
% * MAT: matrix; 
%
% Output(s):
% * ind: array; indicies.
% * Val: array; column values.

ind = find(MAT(:,colInd) ~= 0);
Val = MAT(ind,colInd); 

function LdPwrMod = addShntPwr(LdPwr, ShntPwr, LdInd, ShntInd)
% Add shunt power to the load power.
% Input(s):
% * LdPwr: array; load power.
% * ShntPwr: array; shunt power.
% * LdInd: array; load power indicies.
% * ShntInd: array; shunt power indicies.
%
% Output(s):
% * LdPwrMod: array; modified load power.

LdPwrMod = LdPwr;
indx = ismember( LdInd, ShntInd);
LdPwrMod(indx) = LdPwr(indx) + ShntPwr;
