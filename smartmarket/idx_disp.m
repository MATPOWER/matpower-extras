function [QUANTITY, PRICE, FCOST, VCOST, SCOST, PENALTY] = idx_disp
%IDX_DISP   Defines variables for column indices to dispatch.
%   [QUANTITY, PRICE, FCOST, VCOST, SCOST, PENALTY] = idx_disp

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define the indices
QUANTITY		= 1;	%% quantity produced by generator in MW
PRICE			= 2;	%% market price for power produced by generator in $/MWh
FCOST			= 3;	%% fixed cost in $/MWh
VCOST			= 4;	%% variable cost in $/MWh
SCOST			= 5;	%% startup cost in $
PENALTY			= 6;	%% penalty cost in $

return;
