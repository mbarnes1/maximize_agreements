function f = term1(theta, alpha)
%TERM1 Value of term 1
slack = 1E-6;
validateattributes(alpha+slack, {'numeric'}, {'nonnegative', 'real', 'vector'})
validateattributes(theta, {'numeric'}, {'real', 'nonnegative', 'vector', '<', pi/2})
if abs(sum(alpha)-1) > slack
    error('Alpha must sum to 1')
end
n = length(alpha)-1;
n0 = 0;
alpha = alpha(:)';
theta = theta(:)';

f = bsxfun(@times, sum(bsxfun(@times, alpha, bsxfun(@power, 1-theta(:)/pi, n0:n)), 2), 1./cos(theta'));

validateattributes(f, {'numeric'}, {'nonnegative', 'real', 'vector'})
if length(f) ~= length(theta)
    error('Input and output dimensions must match')
end
end