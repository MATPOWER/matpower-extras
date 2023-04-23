function [Sfe, Ste] = cmptSmat(V, Ybus, branch)
% Compute the complex power matrices (from-to) "Sfe" and the (to-from)
% "Ste.
% Inputs:
% ------------------------------------------------------------------------
% V     : Complex array of bus voltages in the bus system.
% Ybus  : Admittance (Y) matrix of lines in the power system.
% branch: Matrix of branch vlaues produced by the LOADCASE function
% provided by the MATPOWER package.
%
% Outputs:
% ------------------------------------------------------------------------
% Sfe : Complex power matrix of the power (from-to) buses, where the row
% index is the source (from) and column index is the destination (to).
% Ste : Complex power matrix of the power (to-from) buses, where the row
% index is the destination (to) and column index is the source (from).
% 
% by Sami Aldalahmeh 2023/2/15.
%

%%
   % Number of buses
   nb = size(Ybus, 1);
   % Overall complex power matrix
    Sij = V * V'.*conj(Ybus) - (V .* conj(V) * ones(1,nb)).*conj(Ybus);
    
    % Find the branch from-to branch number the branch matric
    brncInd = branch(:,[1 2]);
    
    % Convert to linear indicies
    sijFromInd = sub2ind( size(Sij), brncInd(:,1), brncInd(:,2) );

    % Extract from-to elements
    Sfe = Sij(sijFromInd);
    
    % Find linear indicies for the to-from matrix
    sijToInd = sub2ind( size(Sij), brncInd(:,2), brncInd(:,1) );

     % Extract from-to elements
    Ste = Sij(sijToInd);
    
return
    