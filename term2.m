function f = term2(theta, alpha)
%TERM2 Value of term 2
slack = 1E-6;
validateattributes(alpha+slack, {'numeric'}, {'nonnegative', 'real', 'vector'})
validateattributes(theta, {'numeric'}, {'real', 'nonnegative', 'vector', '<', pi/2+slack})
if abs(sum(alpha)-1) > slack
    error('Alpha must sum to 1')
end
n = length(alpha)-1;
n0 = 0;
alpha = alpha(:)';
theta = theta(:)';

f = bsxfun(@times, 1-sum(bsxfun(@times, alpha, bsxfun(@power, 1-theta(:)/pi, n0:n)), 2), 1./(1-cos(theta')));

validateattributes(f, {'numeric'}, {'nonnegative', 'real', 'vector'})
if length(f) ~= length(theta)
    error('Input and output dimensions must match')
end
end