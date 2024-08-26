function [U] = Gp_mpc(xtrc,ytr,u0,N,lb,ub)
% %% Intializtaion of time step & Prediction horizon

xk = xtrc(:,1:end-1);
uk = xtrc(:,end);
yk = ytr;
n = size(xtrc,1);
k = 1:1:n;
%% define the GP hyperparametrs & mean, cov, lik functions
cov = {@covSEard}; 
sf = 0.1; 
for i =1:size(xtrc,2)
    ell(i,:) = 3;                            % setup the GP
end
% ell = [ell_1;ell_2;ell_3];
hyp0.cov  = log([ell;sf]);
mean = {@meanConst}; a = 0.1;       % m(x) = a*x+b
hyp0.mean = [a];

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 100;                                   % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});
hyp_1 = hyp0; n =size(xtrc,2)-1;
%% minimization of hyperparametrs 
for i =1:n
    for j=1:5
        fprintf('Hyperparameters of %d state minimization iteration: %d\n',i,j);
        if Ncg==0
            hyp = hyp0;
        else
            hyp_1 = minimize(hyp_1,'gp', -Ncg, inf, mean, cov, lik, xtrc, ytr(:,i)); % opt hypers
        end 
    end
hyp.state(i) = hyp_1;
[ymu(i), ys(i)] = gp(hyp.state(i), inf, mean, cov, lik, xtrc, ytr(:,i),xtrc(end,:));
end
ys(i+1) = 0;


%% cost function 
u0= u0*ones(1,N); %%intial input u

fun = @(u0) cost(u0,xk,uk,ymu,ys,hyp,yk);

%% minimization of cost function 
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    % options.StepTolerance=1.000000e-12;

    A = [];
    b = [];
    Aeq = [];
    beq = [];
    % A = 1*eye(5);
    % b = [3,3,3,4,5];
    % fn= @(u) fq(yk(end),u,f);
    [U,fval] = fmincon(fun,u0,A,b,Aeq,beq,lb,ub,[],options)

    display('-----------------------------------------')
end