function [gen, gencost] = off2case(gen, gencost, q, p, max_p, discrete)
%OFF2CASE  Updates case variables gen & gencost from quantity & price offers.
%   [gen, gencost] = off2case(gen, gencost, q, p, max_p, discrete) updates gen &
%   gencost variables based on quantity and price offers and the market type
%   specified. Updates PMIN, PMAX, GEN_STATUS and all cost info except STARTUP
%   and SHUTDOWN. Any quantity offered above max_p will be ignored.

%   by Ray Zimmerman, PSERC Cornell    3/14/00
%   Copyright (c) 1996-2000 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%% save old gencost
oldgencost = gencost;
[ng, np]	= size(q);

%% do conversion
pmin = zeros(ng, 1);
pmax = zeros(ng, 1);
if discrete			%% block offers => piece-wise linear total cost function
	gencost				= zeros(ng, COST + 2*np - 1);
	gencost(:, MODEL)	= PW_LINEAR * ones(ng, 1);
	for i = 1:ng
		qq = q(i, :)';			%% column vector of quantity offers for gen i
		pp = p(i, :)';			%% column vector of price offers for gen i
	
		%% strip zero quantities and prices over max_p
		valid = find(qq & pp <= max_p);
		qq = qq(valid);			%% column vector of quantities of valid offers for gen i
		pp = pp(valid);			%% column vector of prices of valid offers for gen i
		n = length(qq) + 1;
	
		%% form piece-wise linear total cost function, set Pmin & Pmax
		if n > 1		%% otherwise, leave all cost info zero (specifically N)
			%% set Pmin and Pmax
			if gen(i, PMIN) < 0 & gen(i, PMAX) <=0
				pmin(i) = gen(i, PMAX) - sum(qq);
				pmax(i) = gen(i, PMAX);
				xx = [0; cumsum(qq)]' + pmin(i);
			else
				pmin(i) = qq(1);
				pmax(i) = sum(qq);
				xx = [0; cumsum(qq)]';
			end
		
			%% form piece-wise linear total cost function
			yy = [0; cumsum(pp .* qq)]';
			gencost(i, N) = n;
			gencost(i, 		COST:2:( COST + 2*n - 1	)) = xx;
			gencost(i, 	(COST+1):2:( COST + 2*n		)) = yy;
			
		end
	end
else				%% continuous linear offers => quadratic total cost function
	%% valid offers must be increasing & less than max_p
	valid = find( (q(:, 2) - q(:, 1)) & p(:, 2) >= p(:, 1) );
	
	%% discard portion above max_p
	fix = find( p(valid, 2) > max_p );
	fix = valid(fix);
	q(fix, 2) = ( max_p*ones(size(fix, 1), 1) - p(fix, 1) ) ./ ( p(fix, 2) - p(fix, 1) ) ...
				.* ( q(fix, 2) - q(fix, 1) ) + q(fix, 1);
	p(fix, 2) = max_p * ones(size(fix, 1), 1);

	%% compute quadratic coefficients
	c2 = ( p(valid, 2) - p(valid, 1) ) ./ ( q(valid, 2) - q(valid, 1) ) / 2;
	c1 = p(valid, 1) - 2 * c2 .* q(valid, 1);
	gencost						= zeros(ng, COST + 2);
	gencost(valid, MODEL)		= POLYNOMIAL * ones(ng, 1);
	gencost(valid, N:(COST+2))	= [3 * ones(length(valid), 1) c2 c1 zeros(length(valid), 1)];
	pmin = q(:, 1);
	pmax = q(:, 2);
end

% if mkt == 1					%% offer = linear, cost = quadratic
% elseif mkt == 2				%% offer = block, cost = quadratic
% 	gencost				= zeros(ng, COST + 2);
% 	gencost(:, MODEL)	= POLYNOMIAL * ones(ng, 1);
% 	gencost(:, N)		= 3 * ones(ng, 1);
% 	npts				= 50;	%% fit to 50 points
% % 	plotnum = 0;
% 	for i = 1:ng
% 		qq = q(i, :)';			%% column vector of quantity offers for gen i
% 		pp = p(i, :)';			%% column vector of price offers for gen i
% 	
% 		%% strip zero quantities and prices over max_p
% 		valid = find(qq & pp <= max_p);
% 		qq = qq(valid);			%% column vector of quantities of valid offers for gen i
% 		pp = pp(valid);			%% column vector of prices of valid offers for gen i
% 		n = length(qq);
% 	
% 		if n > 1		%% otherwise, leave all cost info zero (specifically N)
% 			%% form piece-wise linear total cost function
% 			xx = [0; cumsum(qq)];
% 			yy = [0; cumsum(pp .* qq)];
% 			pwlin = mkpp(xx, [pp yy(1:n)]);
% 			
% 			%% fit a quadratic to it & put in gencost
% 			xmax = max(xx);		fitstep = xmax / (npts - 1);
% 			fitx = [0:fitstep:xmax];
% 			fity = ppval(pwlin, fitx);
% 			fitpoly = polyfit(fitx, fity, 2);
% 
% % 			%% plot the fit in total and marginal costs
% % 			contstep = 1.01*xmax / 200;
% % 			contx = [0:contstep:(1.01*xmax)]';
% % 			[xs, ys] = stairs(xx, [pp; pp(n)]);		%% marginal cost blocks
% % 			fitpolyder = polyder(fitpoly);			%% find derivative of cost fit
% % 	% 		p1 = polyfit(xx, yy, 2);				%% fit polynomial to corner points of piecewise linear
% % 	% 		pd1 = polyder(p1);						%% take derivative
% % 			plotnum = plotnum + 1;
% % 			subplot(ng, 2, plotnum);
% % 			plot(xs, ys, 'g');
% % 			axis([0 1.1*max(xx) 0 1.1*max(pp)]);
% % 			hold on;
% % 			plot(xs, ys, 'o');
% % 	% 		plot(contx, polyval(pd1, contx), 'y');
% % 			plot(contx, polyval(fitpolyder, contx), 'r');
% % 			grid on;
% % 			hold off;
% % 			plotnum = plotnum + 1;
% % 			subplot(ng, 2, plotnum);
% % 			plot(contx, ppval(pwlin, contx), 'g');
% % 			grid on;
% % 			hold on;
% % 			plot(xx, yy, 'yo');
% % 	% 		plot(contx, polyval(p1, contx), 'y');
% % 			plot(contx, polyval(fitpoly, contx), 'r');
% % 			hold off;
% % 			pause;
% 
% 			gencost(i, COST:(COST+2)) = fitpoly;
% 		end
% 	end
% 	
% 	pmin = q(:, 1);
% 	pmax = sum(q')';
% elseif mkt == 3				%% offer = block, cost = piece-wise linear
% end

%% copy back startup and shutdown costs
gencost(:, [STARTUP SHUTDOWN]) = oldgencost(:, [STARTUP SHUTDOWN]);

%% set PMIN, PMAX, GEN_STATUS
off = find(gencost(:, N) == 0);					%% find gens with no valid offers
gen(off, GEN_STATUS) = zeros(length(off), 1);	%% turn them off
on = find(gen(:, GEN_STATUS) > 0);				%% find generators which are on
% gen(on, PMIN) = pmin(on);
gen(on, PMAX) = pmax(on);

return;
