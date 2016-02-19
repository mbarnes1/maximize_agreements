function [ theta, approx, duality_gap ] = minTheta( alpha, theta0 )
%MINTHETA Minimize over theta, player 2's strategy
%   Inputs:
%       alpha - [n x 1] vector, sum(alpha) = 1, alpha>0
%       theta0 - scalar in (0, pi/2). Optimization initialization value
%
%   Outputs:
%       theta - scalar in [0, pi/2]. Player 2's optimal strategy
%       approx - Approximation ratio, i.e. the dual solution
%       duality_gap - constant >= 0. Gap b/t dual and primal solutions
validateattributes(alpha, {'numeric'}, {'nonnegative', 'vector'})
if abs(sum(alpha)-1) > 1E-12
    error('Alpha must sum to 1')
end
addpath('convexprog')

verbose = false;
maxiter = 1000;
descentdir = 'newton';
tol = 1E-6;

%% Minimize term 1
obj = @(x) term1(x, alpha);
grad = @(x) gradTerm1(x, alpha);
[theta1, primal1, duality_gap1] = ipsolver(theta0, obj, grad, @constraints, @jacobian, descentdir, tol, maxiter, verbose);
approx1 = primal1 - duality_gap1;

%% Minimze term 2
obj = @(x) term2(x, alpha);
grad = @(x) gradTerm2(x, alpha);
[theta2, primal2, duality_gap2] = ipsolver(theta0, obj, grad, @constraints, @jacobian, descentdir, tol, maxiter, verbose);
approx2 = primal2 - duality_gap2;

%% Return min(term1, term2)
if approx1 < approx2
    theta = theta1;
    approx = approx1;
    duality_gap = duality_gap1;
else
    theta = theta2;
    approx = approx2;
    duality_gap = duality_gap2;
end

if duality_gap > tol*10
    error('Duality gap is too large')
end

end


function [G, H] = gradTerm1(theta, alpha)
%GRADTERM1 Gradient with respect to theta
% Inputs
%       theta - scalar in [0, pi/2)
%       alpha - [n x 1] vector, alpha > 0, sum(alpha) = 1
% Outputs
%       G - scalar gradient w.r.t. theta
%       H - scalar Hessian w.r.t. theta

validateattributes(alpha, {'numeric'}, {'nonnegative', 'real', 'vector'})
validateattributes(theta, {'numeric'}, {'scalar', 'real', '<', pi/2})
if abs(sum(alpha)-1) > 1E-12
    error('Alpha must sum to 1')
end
n = length(alpha)-1;
n0 = 0;

G = sum(alpha.*((1-theta/pi).^(n0:n)*tan(theta)*sec(theta)-(n0:n).*(1-theta/pi).^(n0-1:n-1)*sec(theta)/pi));
H = sum(alpha.*((n0-1:n-1).*(n0:n).*(1-theta/pi).^(n0-2:n-2)*sec(theta)/pi^2+...
    (1-theta/pi).^(n0:n)*(sec(theta)^3 + tan(theta)^2*sec(theta)) - ...
    2*(n0:n).*(1-theta/pi).^(n0-1:n-1)*tan(theta)*sec(theta)/pi));
validateattributes(G, {'numeric'}, {'scalar', 'real'})
validateattributes(H, {'numeric'}, {'scalar', 'real', 'nonnegative'})
end


function [G, H] = gradTerm2(theta, alpha)
%GRADTERM2 Gradient with respect to theta
% Inputs
%       theta - scalar in [0, pi/2)
%       alpha - [n x 1] vector, alpha > 0, sum(alpha) = 1
% Outputs
%       G - scalar gradient w.r.t. theta
%       H - scalar Hessian w.r.t. theta

validateattributes(alpha, {'numeric'}, {'nonnegative', 'real', 'vector'})
validateattributes(theta, {'numeric'}, {'scalar', 'real', '<', pi/2})
if abs(sum(alpha)-1) > 1E-12
    error('Alpha must sum to 1')
end
n = length(alpha)-1;
n0 = 0;

G = -sin(theta)/(1-cos(theta))^2 + ...
    sum(alpha.*((n0:n).*(1-theta/pi).^(n0-1:n-1)/(pi*(1-cos(theta))) + ...
    (1-theta/pi).^(n0:n)*sin(theta)/(1-cos(theta))^2));
H = 1/4*(cos(theta)+2)*csc(theta/2)^4 + ...
    sum(alpha.*(-(n0-1:n-1).*(n0:n).*(1-theta/pi).^(n0-2:n-2)/(pi^2*(1-cos(theta))) - ...
    (1-theta/pi).^(n0:n).*(2*sin(theta)^2/(1-cos(theta))^3 - cos(theta)/(1-cos(theta))^2) - ...
    2*(n0:n).*(1-theta/pi).^(n0-1:n-1)*sin(theta)/(pi*(1-cos(theta))^2)));
validateattributes(G, {'numeric'}, {'scalar', 'real'})
validateattributes(H, {'numeric'}, {'scalar', 'real', 'nonnegative'})
end


function C = constraints(theta)
%CONSTRAINTS
C = [-theta; theta - pi/2];  % theta in (0, pi/2)
end


function [J, W] = jacobian(theta, Z)
%JACOBIAN
%   Inputs
%       theta - scalar in [0, pi/2)
%       Z - [2 x 1] vector of dual variables
%   Outputs
%       J - [2 x 1] Jacobian matrix. First-order partial derivatives
%           of the inequality constraint functions
%       W - [1 x 1] Hessian of the Lagrangian (minus the Hessian of the objective)

validateattributes(theta, {'numeric'}, {'scalar', 'real', '<', pi/2})
validateattributes(Z, {'numeric'}, {'vector', 'real'})

J = [-1; 1];
W = 0;
end
