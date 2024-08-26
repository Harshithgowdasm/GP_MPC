

function [J] = costf(ute,yk,Q,R,f,f1,rand)
N = size (ute,2);
x1(1) = yk(end,1);
x2(1) = yk(end,2);
J =0;
for i=1:N
    mu_a(i,1) = f(x1(i),x2(i),ute(i)); %,rand(i,1));
    mu_a(i,2) = f1(x1(i),x2(i),ute(i)); %,rand(i,1));

    J = J + mu_a(i,:)*Q*mu_a(i,:)'+ ute(i)*R*ute(i);
    
    x1(i+1)=mu_a(i,1);
    x2(i+1)=mu_a(i,2);

end 

end



