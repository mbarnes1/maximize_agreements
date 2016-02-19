%MAXIMIZEAGREEMENTS
% Computational proof for an improved correlation clustering approximation
% Matt Barnes
% mbarnes1@cs.cmu.edu
% Carnegie Mellon University
clear, clc, close all
dampening = 0;
tol = 1E-6;
maxiter = 100;
theta = 0.5;
n = 2;
alph = maxAlpha(theta, n);
disp('Iter  Approx   Theta   Alpha')
disp('----------------------------------------------')

iter = 1;
old_approx = 0;
approx = 2*tol;
while iter<maxiter  % approx - old_approx > tol && 
    old_approx = approx;
    alphaUndamped = maxAlpha(theta, n);
    %alph = alphaUndamped;
    alph = (dampening*alph+alphaUndamped)/(1+dampening);
    [ thetaUndamped, approx, duality_gap ] = minTheta( alph, theta );
    %theta = thetaUndamped;
    theta = (dampening*theta + thetaUndamped)/(1+dampening);
    fprintf('% 3u % 8.4f % 7.3f |% 5.3f % 5.3f ... % 3.3f \n', iter, approx, theta, alph(1), alph(2), alph(end))
    iter = iter+1;
end
