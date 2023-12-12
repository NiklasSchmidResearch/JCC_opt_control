clear all
close all

% This file was used to analyse the fish model
x = zeros(100,1);
L=40;
r=1;
%rng(1234)

sigma = 5;
figure(1)
x_rng = linspace(0,100,100);
K = 1 ./ (1 + exp(-(x_rng-20)/sigma) );
plot(x_rng,K)


M = 10;
d = 0;      % No fishing for now

figure(2)
failure_counter = 0;
for iteration=1:1000
    x(1) = 14;
    for k = 1:100
        v = normrnd(0.2, 0.1^2);
        v = boundValue(v,[0,1]);
        gamma = normrnd(1, 0.6^2);
        delta = normrnd(1.1, 0.2^2);
        delta = boundValue(delta,[0,inf]);
    
        C = delta*d*M*x(k)/L;
        C=0;
    
        K = 1 / (1 + exp(-(x(k)-20)/sigma) );
        R = r*x(k)*(1- x(k)/L);
        R = max(R,0);
        x(k+1) = (1-v)*x(k) + gamma*R*K - C;
        x(k+1) = max(x(k+1),0);
        x(k+1) = round(x(k+1));
        x(k+1);
    end
    plot(x)
    if x(end)<10
        failure_counter = failure_counter +1;
    end
    hold on
end
failure_counter/1000
% Population crashes below 13: We desire to end up with at least 13 in the end of
% the simulation so that population might recover. If we go below, catch
% what is left.

