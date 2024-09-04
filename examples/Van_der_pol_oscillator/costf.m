

function [J] = costf(ute,yk,Q,R,ra)
vdp = VanDerPol_uncertanity();
N = size (ute,2);
x = yk(end,:);
final= [1,0];
J =0;
for i=1:N
    mu_a(i,:) = vdp.f_ud(x(i,:)',ute(i),ra(i)'); %,rand(i,1));
    % mu_a(i,2) = f1(x1(i),x2(i),ute(i)); %,rand(i,1));

    J = J + (mu_a(i,:)-final)*Q*(mu_a(i,:)-final)'+ ute(i)*R*ute(i);
    
    x(i+1,:)=mu_a(i,:);
    % x2(i+1)=mu_a(i,2);

end 
% J = sum(diag(mu_a'*Q*mu_a)) + ute*R*ute';
end



