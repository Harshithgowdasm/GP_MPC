function [U] = Gp_mpc(xtrc,ytr,u0,N,lb,ub,Q,R)
% %% Intializtaion of time step & Prediction horizon

xk = xtrc(:,1:end-1);
uk = xtrc(:,end);
yk = ytr;
nx =size(xtrc,2)-1;

%% define the GP hyperparametrs & mean, cov, lik functions
cov = {@covSEard}; 
sf = 0.1; 
for i =1:size(xtrc,2)
    ell(i,:) = 3;                            % setup the GP
end
hyp0.cov  = log([ell;sf]);
mean = {@meanLinear}; a = 0.1;       
hyp0.mean = [a;a;a];  % @meanLinear 
% hyp0.mean = [a];  % @meanConst

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 100;                                   % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});
hyp_1 = hyp0; 

%% optimization of hyperparametrs 
for i =1:nx
    if Ncg==0
        hyp = hyp0;
    else
        hyp_1 = fitGP(hyp_1, inf, mean, cov, lik,xtrc, ytr(:,i),[]);
    end 
    hyp.state(i) = hyp_1;
    [state(i).m, state(i).s2] = gp(hyp.state(i), inf, mean, cov, lik, xtrc, ytr(:,i),xtrc);  % prediction xk+1 for ending data point
    [ymu(i), ys(i)] = gp(hyp.state(i), inf, mean, cov, lik, xtrc, ytr(:,i),xtrc(end,:));
end
ys(i+1) = 0;        % set variance of input u as zero

%% cost function 
u0= u0*ones(1,N); %%intial input u
fun = @(u0) cost(u0,xk,uk,ymu,ys,hyp,yk,Q,R);       % define cost function
%% minimization of cost function using fmincon 
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    [U,fval] = fmincon(fun,u0,A,b,Aeq,beq,lb,ub,[],options)

end