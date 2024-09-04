% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda


function [mu_up,var_up] = GP_transition_change(mu_p, var_p,Kaa, sn, xi, Y,cov_hyp,var_x,utr)

% mu_p : input mean of GP dimension is column vector (N)
% var_p : input variance of GP dimension is sqaure matrix (N*N)
% Kaa : Cov. kernal of Training set a, dimension sqaure matrix (N*N)
% Ki : Cov. kernal of training set a with %%%  dimension sqaure matrix (N*N)
% Sn : Noise variance one of hyper parameters,dimension scalar.
% Xi : test input data set, dimension is column vector (N*1)
% Y : Training output data set, dimension is column vector (N*1)
% hyp1 : lengthscale hyper parameter, hyp.cov(1:end-1), dimension scalar, but
% converted into square matrix of exp(hyp1)
% x_var : variance of input data set, exp(hyp.cov(end)) , scalar 

%% mean calculation
% finding lengthscale, sqaure matrix  
    N = size(xi,1); % N=Nd+u(no. training data of augumented matrix of states + actions)
    n =size(var_p,1); % n = no of states + no of actions
    hyp = diag(exp(cov_hyp)); % Lenthscale hyperparamters
    x_var = exp(var_x)^2;   % variance of X data set matrix
    sn2 = max(sn^2,1e-7);   % noise variance

% finding beta = (Ka + σ2)^−1*ya, column vector (N*1) 
    ik = (Kaa+(sn2)*eye(N));
    beta = ik\Y;

iL = inv(hyp);            % iL dimension square matrix n
B = var_p *(iL*iL')+ eye(n);     % B is symmetric square matrix matrix n

for i= 1:N
    vi(i,:) = xi(i,:) - [mu_p]; % dimension scalar
    ka(i,:) = KernalCov(xi(i,:),[mu_p],[cov_hyp;var_x]); 
    iN = vi(i,:)*iL;      % iN scalar
    lb(i,:) = exp(-0.5*sum(iN*pinv(B)*iN',2)); % square matrix converted to row vector N
    c(i,:) = (x_var/sqrt(det(B))) * lb(i,:);  % each row vector appended to sqaure matrix N*N
end
mu_up = beta' * c; % scalar

%% variance calculations 

for i=1:N
    for j=1:N
        R = var_p*(iL*iL'+iL*iL') + eye(n);   % R= Σ(Λ−1a + Λ−1b) + I ,
        Zij =vi(i,:)*iL'*iL+ vi(j,:)*iL'*iL; %zij := Λ−1a*νi + Λ−1b*νj 

        T_ =(pinv(R)*var_p);   %T_ := inverse(Λ−1a + Λ−1b + Σ−1)

        Qij(i,j) = ka(i,:) * ka(j,:) * exp(0.5*Zij * T_ * Zij')/sqrt(det(R)); % Qij=ka(xi,µ(t−1))*kb(xj ,µ(t−1))*exp(-Zij * T_ * Zij) 
    end
end
var_p = beta' * Qij * beta ;

var_f = x_var-sum(diag(ik\Qij))+sn2;
var_up = var_f+var_p-mu_up^2;
end






