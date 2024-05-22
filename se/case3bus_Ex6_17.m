function [baseMVA, bus, gen, branch, areas, gencost] = case3bus_Ex6_17
%CASE3BUS_Ex3_9 Case of 3 bus system.
%   From Excecise 3.9 in book 'Computational
%   Methods for Electric Power Systems' 2nd Edition by Mariesa Crow page
%   75.
%   created by Sami Aldalahmeh on 2/15/2023
%
%   MATPOWER

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 1000;

%% bus data
%	bus_i	type	Pd	  Qd	Gs	 Bs	 area  Vm	Va	baseKV	zone Vmax	Vmin
bus = [
	   1     3      0     0     0	 0	  1	   1.02	0	230  	1	 1.02	0.99;
	   2     2      0     0     0	 0	  1	   1.00	0	230  	1	 1.02	0.99;
	   3     1      1200  500   0    0	  1	   1.00	0	230  	1	 1.02	0.99;
];

%% generator data
% Note: 
% 1)It's better of gen to be in number order, otherwise gen and genbid
% should be sorted to make the lp solution output clearly(in number order as well)
% 2)set Pmax to nonzero. set to 999 if no limit 
% 3)If change the order of gen, then must change the order in genbid
% accordingly
%	     bus	Pg	  Qg	Qmax  Qmin	Vg	  mBase	   status	Pmax	Pmin
gen = [
	       1	0      0	999	  -999 	1.02  100   	 1	     999	-999;
	       2	500    0	999   -999  1.00  100    	 1	 999	-999;
% 	       3	0      0	999	  -999	1.00  baseMVA	 1	     999	-999;
];
%gen(:, 9) = 999; % inactive the Pmax constraints

%% branch data
%	        fbus	tbus	r	  x	    b	   rateA	rateB	rateC	ratio	angle	status
branch = [
	           1	   2	0.02  0.3   0.150	0	     0  	0	       0	 0	    1;
	           1	   3	0.01  0.1   0.1  	0	     0  	0	       0	 0	    1;
	           2	   3	0.01  0.1   0.1  	0	     0  	0	       0	 0	    1;
];

%%-----  OPF Data  -----%%
%% area data
areas = [
	1	1;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	1	1.5 	1	0;
	2	0	0	2	1   	2	0;
];
return;
