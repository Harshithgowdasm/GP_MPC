% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda

clear all;
close all;

%% Intializtaion of time step & Prediction horizon
ta = 0; tb = 10; ts = 0.3;
t1 = ta+ts:ts:tb+ts;
t = ta:ts:tb;
N =50; %prediction horizon
t2 = tb+ts : ts : tb+N*ts;

%% Intialization of States and U input as functions for data 
% f = @(t) sin(2*t)+0.00002*randn(size(t));
f = @(t) exp(0.1*t)+0.00002*randn(size(t));
% g= @(t) exp(0.01*t)+0.00002*randn(size(t));
% f = @(t) 10*t+0.00002*randn(size(t));
g = @(t) 10* t+0.00002*randn(size(t));

% utr = 2*ones(51,1);
utr= g(t)'; % input vector uk
xtr = f(t)'; % states vector xk

% ytr = 1*f(t)'+1*utr; 
ytr = f(t1)'; %output discrete xk+1
ntr = size(xtr,1); %number of data points 

xtrc=[xtr,utr]; % aug. xk,uk


%% define the GP hyperparametrs & mean, cov, lik functions
cov = {@covSEard}; sf = 1; ell_1 = 5;
ell_2 = 4;                              % setup the GP
ell = [ell_1;ell_2];
hyp0.cov  = log([ell;sf]);
mean = {@meanConst}; a = 0.3;       % m(x) = a*x+b
hyp0.mean = [a];

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 50;                                   % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});

%% minimization of hyperparametrs 
if Ncg==0
    hyp = hyp0;
else
    hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtrc, ytr); % opt hypers
end

%% prediction of xk+1 for last xk point 
[ymu2, ys2] = gp(hyp, inf, mean, cov, lik, xtrc, ytr,xtrc(end,1:2) );  % predict % Xtr = [X,U)
% [nlZ,dnlZ] = gp(hyp, inf, mean, cov, lik, xtrc, ytr);

%% longterm prediction xk+1 over prediction horizon 
ute = g(t2)';
X_tilda = [xtr,utr];
for i=1:N

    if i==1
        mu_p = [ymu2,ute(i)]; % µ_tilmda = [µt,ut]
        var_p = blkdiag(ys2,0); % Σ_tilda = blkdiag[Σt, 0] 
    else
        mu_p = [mu_up(i-1),ute(i)];
        var_p = blkdiag(var_up(i-1),0);
    end

    Kaa = KernalCov(X_tilda,X_tilda,hyp.cov);

    [mu_up(i),var_up(i)] = Gp_transition_change(mu_p, var_p,Kaa,exp(hyp.lik), X_tilda, ytr(1:end),hyp.cov(1:end-1),hyp.cov(end),utr);


end

%% plots predicted vs original function 
figure, hold on
plot(t,ytr,'Color','r','LineWidth',2)

plot(t2,mu_up,'Color','g','LineWidth',2)
plot(t2,f(t2+ts),'Color','r','LineWidth',2)

hold off



