

function [J] = cost(ute,xtr,utr,ymu,ys,hyp,ytr)
N = size (ute,2);
xtrc = [xtr,utr];
for i=1:N

   if i==1
        mu_p = [ymu,ute(i)]; % µ_tilmda = [µt,ut]
        var_p = diag(ys); % Σ_tilda = blkdiag[Σt, 0] 
    else
        mu_p = [mu_a(i-1,:),ute(i)];
        var_p = blkdiag(var_a(:,:,i-1),0);
    end     

   [mu_a(i,:),var_a(:,:,i)] = Gp_transition_change(mu_p, var_p,hyp,xtrc, ytr);

end
    % disp("--------")  
Q=1*eye(N);
R=0.1*eye(N);
J = sum(diag(mu_a'*Q*mu_a)) + ute*R*ute';

end



