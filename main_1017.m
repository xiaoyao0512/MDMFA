clear; close all; clc

N = 100;
A = zeros(2, 2, 1);
A(1, 2, 1) = 1;
A(2, 1, 1) = 0;

order = [0,0];

p = 1;

[T, K] = generate_time_series(A, p, N, order);

% scmin=4;
% scmax=512;
% scres=8; % 14 CONSIDER 20 scales
%  
% exponents=linspace(log2(scmin),log2(scmax),scres);
% Scales=round(2.^exponents);
% 
% Q = -5:0.1:5; % 0.5
% 
% [Fxy, Hq, tq, hq, Dq] = MFDFA(T(1,:), Scales, Q); % Hq should be around 0.5;
% 
% plot(Q, Hq);


[causality, pval] = mdmfa(T(1,:), T(2,:), N);