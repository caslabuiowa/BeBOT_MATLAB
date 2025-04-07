function [Cpout, Pos, tknotsout] = PiecewisedeCasteljau(Cp, tknots, lambda)
% PIECEWISEDECASTELJAU Evaluates a piecewise Bernstein polynomial at a given lambda
% 
% Inputs:
%   Cp - Control points of the composite Bernstein polynomial
%   tknots - Knot sequence defining piecewise polynomial intervals
%   lambda - The evaluation point in the range [0,1]
%
% Outputs:
%   Cpout - Updated control points incorporating the evaluated result
%   Pos - The position of lambda in the Bernstein polynomial

    % Normalize tknots to the range [0,1] for consistent evaluation
    tknots_norm = tknots / tknots(end);

    % Check if lambda is already in tknots
    if any(abs(tknots_norm - lambda) < 1e-5)
        Cpout = Cp;
        Pos = find(abs(tknots_norm - lambda) < 1e-10, 1);
        tknotsout = tknots;
        return;
    end
    
    % Find the interval where lambda belongs
    k = find(tknots_norm <= lambda, 1, 'last');
    
    % Determine the polynomial degree
    n = length(Cp)/(length(tknots_norm)-1) - 1;
    
    % Extract Bernstein coefficients for the identified interval
    B_coeffs = Cp((k-1)*(n+1) + 1: k*(n+1));
    local_lambda = (lambda - tknots_norm(k)) / (tknots_norm(k+1) - tknots_norm(k)); 
    
    % Apply De Casteljau's algorithm to compute the value at lambda
    [result, Pos] = deCasteljau(B_coeffs, local_lambda);
    mid_idx = ceil(length(result) / 2);
    result = [result(1:mid_idx), result(mid_idx), result(mid_idx+1:end)];
    
    % Insert lambda into the knot vector and rescale to original range
    tknots_norm = sort([tknots_norm, lambda]);
    tknotsout = tknots_norm * tknots(end);
    
    % Update control points by inserting the computed result
  
    Cpout = [Cp(1:(k-1)*(n+1)), result, Cp(k*(n+1) + 1:end)];
end





%% Example code
% clear all
% close all
% 
% Cp = [0 1 -1 2 2 1 0 1 1 3 3 3 3 1 0 -1];
% N = 3;
% tknots = [0 1 3 3.5 6];
% t = linspace(0,tknots(end),1000);
% lambda = 0.6;
% 
% [tnodes,~,~] = PiecewiseBeBOT(N,tknots);
% 
% figure 
% plot(tnodes,Cp,'o','color','r') 
% hold on
% plot(t,PiecewiseBernsteinPoly(Cp,tknots,t),'color','r') 
% 
% [Cpout, Pos, tknotsout] = PiecewisedeCasteljau(Cp, tknots, lambda);
% 
% [tnodes,~,~] = PiecewiseBeBOT(N,tknotsout);
% 
% 
% plot(tnodes,Cpout,'o','color','b') 
% plot(t,PiecewiseBernsteinPoly(Cpout,tknotsout,t),'color','b') 