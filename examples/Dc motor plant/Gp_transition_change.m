% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda

function [mu_up,var_up] = Gp_transition_change(mu_p, var_p,hyp, xi, Y)

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
    N = size(xi, 1);
    n = size(xi, 2);
    mu_up = zeros(1, n-1);
    var_up = zeros(n-1);

    state = struct('hyp', {}, 'x_var', {}, 'sn', {}, 'K', {}, 'ik', {}, 'beta', {}, ...
                   'iL', {}, 'B', {}, 'k', {}, 'Qij', {});
    vi = xi-mu_p;
    c = zeros(N, 1,size(xi,2)-1);    
    var_p = topdm(var_p);
for k=1:size(xi,2)-1
    for l=1:size(xi,2)-1
        state(l).hyp = diag(exp(hyp.state(l).cov(1:end-1)));  % Lenthscale hyperparamters
        state(l).x_var = exp(hyp.state(l).cov(end))^2;   % variance of X data set matrix
        state(l).sn = max(exp(hyp.state(l).lik)^2,1e-7); % noise variance    
        state(l).K = KernalCov(xi,xi,hyp.state(l).cov);
        state(l).ik = (state(l).K+state(l).sn*eye(N));
        state(l).beta = state(l).ik\Y(:,l);
        state(l).iL = inv(state(l).hyp);            % iL dimension square matrix n
        state(l).B = var_p *(state(l).iL*state(l).iL')+ eye(n);     % B is symmetric square matrix matrix n
        state(l).k = KernalCov(xi,mu_p,hyp.state(l).cov);

        c(:,:,k) = mu(vi,state(k).iL,state(k).B,state(k).x_var);
        state(k).Qij(:,:,l) = var(var_p,state(k).iL,state(l).iL,n,vi,state(k).k,state(l).k);
    end
    mu_up(1,k) = state(k).beta' * c(:,:,k); % scalar
    if ~isreal(mu_up(1,k))
        o9=1;
    end
end
%% variance calculations 
for i = 1:size(xi,2)-1
    for j = 1:size(xi,2)-1
        if i==j
            var_up(i,j) = state(i).x_var-sum(diag(state(i).ik\state(i).Qij(:,:,i)))+state(i).sn+state(i).beta' * state(i).Qij(:,:,i) * state(i).beta-mu_up(1,i)^2;
        elseif j>i
            var_up(i,j) = state(i).beta' * state(i).Qij(:,:,j) * state(j).beta - mu_up(1,i)* mu_up(1,j);
            var_up(j,i) = var_up(i,j);
        end
    end
    % display(var_up);
end
end


function [Qij] = var(var_p,iL1,iL2,n,vi,k1,k2)

    R = var_p*(iL1*iL1'+iL2*iL2') + eye(n);   % R= Σ(Λ−1a + Λ−1b) + I ,
    Zij = vi.*sum(iL1^2,1) + vi.*sum(iL2^2,1); %zij := Λ−1a*νi + Λ−1b*νj 
    T =(R\var_p);   %T_ := inverse(Λ−1a + Λ−1b + Σ−1)
    % % T = inv(iL1*iL1'+ iL2*iL2' + inv(var_p))
    Qij =k1 * k2' .* exp(0.5*Zij * (T * Zij'))/sqrt(det(R)); % Qij=ka(xi,µ(t−1))*kb(xj ,µ(t−1))*exp(-Zij * T_ * Zij)

end

function c = mu(vi, iL1, B1, x_var1)
    % Compute iN1 for all i simultaneously
    iN1 = vi * iL1'; % vi is assumed to be a matrix with rows corresponding to different i

    % Compute lb1 using matrix operations
    lb1 = exp(-0.5*sum((iN1 * pinv(B1)) .* iN1, 2));

    % Compute c using vectorized operations
    c = (x_var1 ./ sqrt((det(B1)))) .* lb1; % Element-wise division and multiplication
end




function [sigma] = topdm(sig) 
 
EPS = 10^-6;        
ZERO = 10^-10;    
sigma = sig; 
[~, err] = cholcov(sigma, 0);           
                               
if (err ~= 0) 
    [v d] = eig(sigma); 
    % display(d);
    d=diag(d);          
    d(d<=ZERO)=EPS; 
     
    d=diag(d);       
    sigma = v*d*v'; 
                    
end 
end