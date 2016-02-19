function [ alpha ] = maxAlpha( theta, n )
%MAXALPHA Maximize over alpha, player 1's strategy
%   Inputs
%       theta - scalar in [0, pi/2]
%       n - problem size, number of random hyperplane projections
%
%   Outputs
%       alpha - optimal solution

validateattributes(theta, {'numeric'}, {'>=', 0, '<=', pi/2})
validateattributes(n, {'numeric'}, {'integer', 'positive'})
options = optimset('Display','off');

%% Optimizing terms 1 and 2
f = (1-theta/pi).^(1:n);
lb = zeros(n, 1);
ub = ones(n,1);
Aeq = ones(1,n);
beq = 1;
[alpha1, ~, exitflag, ~] = linprog(-f, [], [], Aeq, beq, lb, ub, [], options);
alpha1 = alpha1';
if exitflag ~= 1
    error('Did not converge to alpha solution')
end
[alpha2, ~, exitflag, ~] = linprog(f, [], [], Aeq, beq, lb, ub, [], options);
alpha2 = alpha2';
if exitflag ~= 1
    error('Did not converge to alpha solution')
end
if term1(theta, alpha1) <= term2(theta, alpha1)
    alpha = alpha1;
elseif term2(theta, alpha2) <= term1(theta, alpha2)
    alpha = alpha2;
else
    Aeq = [Aeq;
           (1-theta/pi).^(1:n)];
    beq = [beq;
           cos(theta)];
    [alpha, ~, exitflag, ~] = linprog(-f, [], [], Aeq, beq, lb, ub, [], options);
    alpha = alpha';
    if exitflag ~= 1
        error('Did not converge to alpha solution')
    end
    if abs(term1(theta, alpha) - term2(theta, alpha)) > 1E-6
        error('Terms should be equal')
    end
end

validateattributes(alpha, {'numeric'}, {'vector', 'real', 'nonnegative'})
if abs(sum(alpha)-1) > 1E-12
    error('Alpha must sum to 1')
end
end

