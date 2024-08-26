function [U] = Gp_mpc(xtrc,yk,u0,N,lb,ub)
%% states and input 
xk = xtrc(:,1);
uk = xtrc(:,2);
n = size(xtrc,1);
k = 1:1:n;
%% define the GP hyperparametrs & mean, cov, lik functions
cov = {@covSEard}; sf = 0.9; ell_1 = 1;
ell_2 = 1;                              % setup the GP
ell = [ell_1;ell_2];
hyp0.cov  = log([ell;sf]);
mean = {@meanConst}; a = 0.3;       % m(x) = a*x+b
hyp0.mean = [a];

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 50; % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});

%% minimization of hyperparametrs 
if Ncg==0
    hyp = hyp0;
else
    hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtrc, yk); % opt hypers
end

%% prediction of xk+1 using GP
[ymu2, ys2] = gp(hyp, inf, mean, cov, lik, xtrc, yk,xtrc(end,1:2) );  % predict % Xtr = [X,U)

%% cost function 
u0= u0*ones(1,N); %%intial input u

fun = @(u0) cost(u0,xk,uk,ymu2,ys2,hyp,yk);

%% minimization of cost function 
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    % options.StepTolerance=1.000000e-12;

    A = [];
    b = [];
    Aeq = [];
    beq = [];
    % A = 1*eye(5);
    % b = [3,3,3,4,5];
    [U,fval] = fmincon(fun,u0,A,b,Aeq,beq,lb,ub,[],options)

end