function [Sfe, Ste] = cmptSmat(V, Ybus, branch)
% Compute the complex power matrices (from-to) "Sfe" and the (to-from)
% "Ste".

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