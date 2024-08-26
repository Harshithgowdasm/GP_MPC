clear all;
close all;

% %% Intializtaion of time step & Prediction horizon
% ta = 0; tb = 20; ts = 0.3;
% t1 = ta+ts:ts:tb+ts;
% t = ta:ts:tb;
% N =50; %prediction horizon
% t2 = tb+ts : ts : tb+N*ts;
% 
% %% Intialization of States and U input as functions for data 
% % f = @(t) sin(1*t)+0.00002*randn(size(t));
% % f = @(t) exp(0.08*t)+0.00002*randn(size(t));
% % g= @(t) exp(0.01*t)+0.00002*randn(size(t));
% f1 = @(t) exp(0.1*t)+0.00002*randn(size(t));
% f = @(t) 4*t+0.00002*randn(size(t));
% f2 = @(t) 0.2*t+0.00002*randn(size(t));
% g = @(t) 1* t+0.00002*randn(size(t));
% 
% % utr = 2*ones(51,1);
% utr= g(t)'; % input vector uk
% xtr = [f(t)',f1(t)',f2(t)']; % states vector xk
% 
% % ytr = 1*f(t)'+1*utr; 
% ytr = [f(t1)',f1(t1)',f2(t1)']; %output discrete xk+1
% ntr = size(xtr,1); %number of data points 
% 
% xtrc=[xtr,utr]; % aug. xk,uk
% xte = [f(t2)',f1(t2)',f2(t2)',g(t2)'];

sysc = ss(tf(21,[1.1,1,0]));
T = 0.1;
sysd = c2d(sysc,T,"zoh");
[A,B,C]=ctrbf(sysd.A,sysd.B, sysd.C);
% f = @(x1,x2,u) sysd.A(1,:) *[x1;x2] + sysd.B(1,:)*u + 0.00002*randn(size(u));
% f1 = @(x1,x2,u) sysd.A(2,:) *[x1;x2] + sysd.B(2,:)*u + 0.00002*randn(size(u));

A_d = A+0.1*rand(size(A))*A;
B_d = B+0.1*rand(size(B)).*B;
f = @(x1,x2,u) A_d(1,:) *[x1;x2] + B_d(1,:)*u; %+ 0.02*randn(1);
f1 = @(x1,x2,u) A_d(2,:) *[x1;x2] + B_d(2,:)*u; %+ 0.02*randn(1);

N = 5; %prediction horizon
n = 25; % no of data points
k = 1:1:n+N;
uk_f= randn(50,1);
x1 = 1; x2 = 1;
for i=0:n+N-1
    yk_f(i+1,1) = f(x1(i+1),x2(i+1),uk_f(i+1));
    yk_f(i+1,2) = f1(x1(i+1),x2(i+1),uk_f(i+1));
    x1(i+2) = yk_f(i+1,1);
    x2(i+2) = yk_f(i+1,2);

end
xk_f = [x1(1:n+N)',x2(1:n+N)'];
% save('data_dc_motor_3',"xk","yk","uk")
% load('data_dc_motor.mat');
xk = xk_f(1:n,:);
yk = yk_f(1:n,:);
uk = uk_f(1:n);
xtrc=[xk,uk];
%% define the GP hyperparametrs & mean, cov, lik functions
cov = {@covSEard}; 
sf = 0.1; 
for i =1:size(xtrc,2)
    ell(i,:) = 3;                            % setup the GP
end
hyp0.cov  = log([ell;sf]);
mean = {@meanLinear}; a = 0.1;       % m(x) = a*x+b
hyp0.mean = [a;0.1;a];

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 50;                                   % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});
hyp_1 = hyp0; n =size(xtrc,2)-1;
%% minimization of hyperparametrs 
for i =1:n
    if Ncg==0
        hyp = hyp0;
    else
        hyp_1 = fitGP(hyp_1, inf, mean, cov, lik,xtrc, yk(:,i),[]);
    end 
    hyp.state(i) = hyp_1;
    [ymu(i), ys(i)] = gp(hyp.state(i), inf, mean, cov, lik, xtrc, yk(:,i),xtrc(end,:));
end
ys(i+1) = 0;
%% 

figure(1), hold on
plot(k(1:25),yk(:,1),'Color','r','LineWidth',2)
plot(k(1:25),yk(:,2),'Color','#A2142F','LineWidth',2)
plot(k(1:25),uk,'Color','k','LineWidth',2)

%% longterm prediction xk+1 over prediction horizon 
ute = uk_f(25:end);
for i=1:6

    if i==1
        mu_p = [ymu,ute(i)]; % µ_tilmda = [µt,ut]
        var_p = diag(ys); % Σ_tilda = blkdiag[Σt, 0] 
        mu_pb = mu_p;
        var_pb = var_p;
    else
        mu_p = [mu_a(i-1,:),ute(i)];
        var_p = blkdiag(var_a(:,:,i-1),0);
        % mu_pb = [mu_b(i-1,:),ute(i)];
        % var_pb = blkdiag(var_b(:,:,i-1),0);
    end     

   [mu_a(i,:),var_a(:,:,i)] = Gp_transition_change(mu_p, var_p,hyp,xtrc, yk)
   % [mu_b(i,:),var_b(:,:,i)] = Gp_transition_change_old_old(mu_pb, var_pb,hyp,xtrc, yk);

end

%% plots predicted vs original function 
 plot(k(25:end),mu_a,'Color','g','LineWidth',2);
 % plot(k(25:end),mu_b,'Color','b','LineWidth',2);
 plot(k(25:end),yk_f(25:end,:),'Color','r','LineWidth',2);
 % plot(t2,mu_b,'Color','y','LineStyle','--');

% RMSE = sum(sqrt(sum(([f(t2+ts)',f1(t2+ts)']-mu_a).^2)))

% figure(2), hold on
% plot(t,ytr,'Color','r','LineWidth',2)
% 
% plot(t2,mu_up,'Color','g','LineWidth',2)
% plot(t2,f(t2+ts),'Color','r','LineWidth',2)
% 
% hold off



