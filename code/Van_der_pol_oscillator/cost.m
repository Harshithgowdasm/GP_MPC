

function [J] = cost(ute,xtr,utr,ymu,ys,hyp,ytr,Q,R)
N = size (ute,2);
n =size(xtr,2);
xtrc = [xtr,utr];
J =0;

final= [0,0]; % Final set point / target point



for i=1:N
   if i==1
        mu_p = [ymu,ute(i)]; % µ_tilmda = [µt,ut]
        var_p = diag(ys); % Σ_tilda = blkdiag[Σt, 0] 

    else
        mu_p = [mu_a(i-1,:),ute(i)];
        var_p = blkdiag(var_a(:,:,i-1),0);

   end     
   [mu_a(i,:),var_a(:,:,i)] = Gp_transition_change(mu_p, var_p,hyp,xtrc, ytr);
    J = J + (mu_a(i,:)-final)*Q*(mu_a(i,:)-final)'+ ute(i)*R*ute(i);

end 

end



