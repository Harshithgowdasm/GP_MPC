

function [J] = cost(ute,xtr,utr,ymu2,ys2,hyp,ytr)
N = size (ute,2);
X_tilda = [xtr,utr];

%% long term prediction by Moment Matching
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

%% define cost

Q=10*eye(N);
R=0.01*eye(N);    
J = mu_up*Q*mu_up'+ute*R*ute';
end



