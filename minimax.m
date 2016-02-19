%MINIMAX
n = 5;  % number of player 1 moves
ntrials = 50;
options = optimset('TolX', 1E-6, 'TolCon', 1E-6);

bestapprox = 0;
bestalpha = NaN;
besttheta = NaN;

for i = 1:ntrials
    alpha0 = rand(n,1);
    alpha0 = alpha0/sum(alpha0);
    theta = linspace(1E-5, pi/2-1E-5, 100);

    F = @(x) -[term1(theta, x); term2(theta, x)];
    Aeq = ones(1,n);
    beq = 1;
    lb = zeros(n, 1);
    [alpha,fval,maxfval,exitflag, output] = fminimax(F, alpha0, [], [], Aeq, beq, lb, [], [], options);
    disp(['With fminimax achieved approx of ', num2str(min(-fval))])
    if(any(alpha<0))
        alpha = alpha - min(alpha(alpha<0));
    end
    alpha = alpha'/sum(alpha);

    [ theta, approx, duality_gap ] = minTheta(alpha, pi/4);
    disp(['After final theta minimization achieved approximation of ', num2str(approx)])
    if approx > bestapprox
        bestapprox = approx;
        besttheta = theta;
        bestalpha = alpha;
    end
end
disp(['The best approx was ', num2str(bestapprox), ' at theta=', num2str(besttheta), ' and alpha '])
bestalpha



