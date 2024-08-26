function [mu_up,var_up] = Gp_transition_change_old(mu_p, var_p,hyp, xi, Y)

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
    vi = zeros(N, size(xi, 2));
    c = zeros(N, 1);    
 
for i=1:size(xi,2)-1
    state(i).hyp = diag(exp(hyp.state(i).cov(1:end-1)));  % Lenthscale hyperparamters
    state(i).x_var = exp(hyp.state(i).cov(end))^2;   % variance of X data set matrix
    state(i).sn = max(exp(hyp.state(i).lik)^2,1e-7); % noise variance    
    state(i).K = KernalCov(xi,xi,hyp.state(i).cov);
    X = [xi(:,i),xi(:,end)];
    cove = [hyp.state(i).cov(i),hyp.state(i).cov(end-1),hyp.state(i).cov(end)];
    state(i).K1 = KernalCov(X,X,cove);
    state(i).ik = (state(i).K+state(i).sn*eye(N));
    state(i).beta = state(i).ik\Y(:,i);
    state(i).iL = inv(state(i).hyp);            % iL dimension square matrix n
    state(i).B = var_p *(state(i).iL*state(i).iL')+ eye(n);     % B is symmetric square matrix matrix n

end
    vi = xi - mu_p;
for k=1:size(xi,2)-1
    for l=1:size(xi,2)-1
        for i= 1:N
            for j=1:N
                if k==1 && i==1
                     % dimension scalar
                    state(l).k(j,:) = KernalCov(xi(j,:),mu_p,hyp.state(l).cov);
                end
                state(k).Qij(i,j,l) = var(var_p,state(k).iL,state(l).iL,n,vi,state(k).k,state(l).k,i,j);
            end
            c(i,:) = mu(vi,state(k).iL,state(k).B,state(k).x_var,i);
        end
    mu_up(1,k) = state(k).beta' * c; % scalar
    end
end
%% variance calculations 
for i = 1:size(xi,2)-1
    for j = 1:size(xi,2)-1
        if i==j
            var_up(i,j) = state(i).x_var-sum(diag(state(i).ik\state(i).Qij(:,:,i)))+state(i).sn+state(i).beta' * state(i).Qij(:,:,i) * state(i).beta-mu_up(1,i)^2;
        else
            var_up(i,j) = state(i).beta' * state(i).Qij(:,:,j) * state(j).beta - mu_up(1,i)* mu_up(1,j);
        end
    end
end
end

function [Qij] = var(var_p,iL1,iL2,n,vi,k1,k2,i,j)

    R = var_p*(iL1*iL1'+iL2*iL2') + eye(n);   % R= Σ(Λ−1a + Λ−1b) + I ,
    Zij =vi(i,:)*iL1'*iL1+ vi(j,:)*iL2'*iL2; %zij := Λ−1a*νi + Λ−1b*νj 
    T =(pinv(R)*var_p);   %T_ := inverse(Λ−1a + Λ−1b + Σ−1)
    Qij = k1(i,:) * k2(j,:) * exp(0.5*Zij * T * Zij')/sqrt(det(R)); % Qij=ka(xi,µ(t−1))*kb(xj ,µ(t−1))*exp(-Zij * T_ * Zij)

end


function [c] = mu(vi,iL1,B1,x_var1,i)
    iN1 = vi(i,:)*iL1;      % iN scalar
    lb1(i,:) = exp(-0.5*sum(iN1*pinv(B1)*iN1',2)); % square matrix converted to row vector N
    c = (x_var1/sqrt(det(B1))) * lb1(i,:);  % each row vector appended to sqaure matrix N*N
    
end 


































































    % vi(i,:) = xi(i,:) - mu_p; % dimension scalar
    % k1(i,:) = KernalCov(xi(i,:),mu_p,hyp_1.cov); 
    % k2(i,:) = KernalCov(xi(i,:),mu_p,hyp_2.cov);
    % iN1 = vi(i,:)*iL1;      % iN scalar
    % iN2 = vi(i,:)*iL2;
    % lb1(i,:) = exp(-0.5*sum(iN1*pinv(B1)*iN1',2)); % square matrix converted to row vector N
    % lb2(i,:) = exp(-0.5*sum(iN2*pinv(B2)*iN2',2));
    % cq1(i,:) = (x_var1/sqrt(det(B1))) * lb1(i,:);  % each row vector appended to sqaure matrix N*N
    % cq2(i,:) = (x_var2/sqrt(det(B2))) * lb2(i,:);



% Qij1 = zeros(N, N);  % m-by-n matrix (assuming m and n are defined)
% Qij12 = zeros(N, N); % m-by-n matrix
% Qij2 = zeros(N, N);  % m-by-n matrix
% for i=1:N
%     for j=1:N
%         Q1(i,j) = var(var_p,iL1,iL1,n,vi,k1,k1,i,j);
%         Q12(i,j) = var(var_p,iL1,iL2,n,vi,k1,k2,i,j);
%         Q2(i,j) = var(var_p,iL2,iL2,n,vi,k2,k2,i,j);
%     end
% end



        % R1 = var_p*(iL1*iL1'+iL1*iL1') + eye(n);   % R= Σ(Λ−1a + Λ−1b) + I ,
        % Zij1 =vi(i,:)*iL1'*iL1+ vi(j,:)*iL1'*iL1; %zij := Λ−1a*νi + Λ−1b*νj 
        % T_1 =(pinv(R1)*var_p);   %T_ := inverse(Λ−1a + Λ−1b + Σ−1)
        % Qij1(i,j) = k1(i,:) * k1(j,:) * exp(0.5*Zij1 * T_1 * Zij1')/sqrt(det(R1)); % Qij=ka(xi,µ(t−1))*kb(xj ,µ(t−1))*exp(-Zij * T_ * Zij) 
        % 
        % R2 = var_p*(iL1*iL1'+iL2*iL2') + eye(n);   % R= Σ(Λ−1a + Λ−1b) + I ,
        % Zij2 =vi(i,:)*iL1'*iL1+ vi(j,:)*iL2'*iL2; %zij := Λ−1a*νi + Λ−1b*νj 
        % T_2 =(pinv(R2)*var_p);   %T_ := inverse(Λ−1a + Λ−1b + Σ−1)
        % Qij12(i,j) = k1(i,:) * k2(j,:) * exp(0.5*Zij2 * T_2 * Zij2')/sqrt(det(R2));
        % 
        % R3 = var_p*(iL2*iL2'+iL2*iL2') + eye(n);   % R= Σ(Λ−1a + Λ−1b) + I ,
        % Zij3 =vi(i,:)*iL2'*iL2+ vi(j,:)*iL2'*iL2; %zij := Λ−1a*νi + Λ−1b*νj 
        % T_3 =(pinv(R3)*var_p);   %T_ := inverse(Λ−1a + Λ−1b + Σ−1)
        % Qij2(i,j) = k2(i,:) * k2(j,:) * exp(0.5*Zij3 * T_3 * Zij3')/sqrt(det(R3));