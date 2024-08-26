function [ymu,ys2,nlZ,dnlZ] = VectorGP(X,Y,Xt,n,N)

%X : vector data with inputs
% Y : vecor data of outputs
% Xt : input vector for prediction 
% n : number of states
% N : no of data sets

Ncg = 50;        % number of conjugate gradient steps


% setup the GP parameters
cov = {@covSEard}; sf = 1; ell_1 = 0.4;   %cov function
hyp0.cov  = log([ell_1;sf]);               % cov hyper parameters
mean = {@meanConst}; a = 1;       % mean function & hyp =a
hyp0.mean = a;

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'};   % Gauss inference
sn = 0.2;        % Gauss dist. stadard deviation
hyp0.lik  = [log([sn])];    % Gauss lik hyper

for i = 0:(n-1)
    xtr = X(i*N+1:i*N+N,1);   % each input(u) data
    ytr = Y(i*N+1:i*N+N,1);    % each ouput(states) data
    xte = Xt(i*N+1:i*N+N,1);    % prediction input

    % gradient decent for optimal hyper parameters
    if Ncg==0
        hyp = hyp0;
    else
        hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtr, ytr); % opt hypers
    end

    % modelling and prediction 
    [ymu{i+1}, ys2{i+1}] = gp(hyp, inf, mean, cov, lik, xtr, ytr, xte);  % predict
    [nlZ{i+1}, dnlZ{i+1}] = gp(hyp, inf, mean, cov, lik, xtr, ytr);
end    