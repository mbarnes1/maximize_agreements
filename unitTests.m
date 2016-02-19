clear, clc, close all
%% Swamy 2004
alpha = [0 0 1];
[theta, approx, ~] = minTheta(alpha, pi/4);
if abs(approx - 0.75) > 1E-4
    error('Does not achieve same results as Swamy''s analytical solution')
end

%% Charikar et al. 2003
alpha = [0 0 1-.1316 0.1316];
[theta, approx, ~] = minTheta(alpha, pi/4);
if abs(approx - 0.7664) > 1E-4
    error('Does not achieve same results as Charikar et al.''s analytical solution')
end