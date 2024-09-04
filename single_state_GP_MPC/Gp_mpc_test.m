% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda


clear all;
close all;

%% Intialization of States and U input as functions for data 
N = 10; %prediction horizon
n = 40; % no of data points
k = 1:1:n;
uk = randn(n,1);
% uk = 1*ones(1,n);
% uk = 0.05:0.05:n*0.05;
x = 1;
% f=@(x,u) -0.5*(1 + x)*cos(x)+1.8*cos(x) + 1*u;
f = @(x,u) 0.99*x + 0.1*u;
% f=@(x,u) 1*(-0.5*(1 + x)*x+0.8*x + 2*u);
% f = @(x,u) 1*(x-u)^2+0.1*(x)^2;%+0.0002*randn(1);
for i=0:n-1
    yk(i+1,1) = f(x(i+1),uk(i+1));
    x(i+2) = yk(i+1,1);
end
xk = x(1:n)';
% uk = uk';
xtrc=[xk,uk]; % aug. xk,uk

%% ploting data points

figure(1), hold on
plot(k,yk,'Color','r','LineWidth',2);
plot(k,uk,'Color','k','LineWidth',2);
hold off

%% Model Learning & Generating Controller(MPC)
u0=1*uk(end);
lb = -3*ones(1,N);
ub = 3*ones(1,N);

for i=1:5
    display(i,'iteration');
    [U] = Gp_mpc(xtrc,yk,u0,N,lb,ub);
    u0=U(1);
    xtrc(size(xtrc,1)+1,1:2) = [yk(end),U(1)];
    yk(size(yk,1)+1) = f(yk(end),U(1));
    
    display(u0,'Optimal U');
    display(yk(end),'yk');
    display(xk(end),'xk');
    display('---------------------------');

    k = n:1:size(xtrc,1);
    figure(1), hold on
    plot(k,yk(n:end),'Color','g','LineWidth',1);
    plot(k,xtrc(n:end,2),'Color','b','LineWidth',1);
    legend('yk','uk','opt. Yk','opt. Uk');
    hold off

end
projected_yk = yk(n:end)
optimal_uk = xtrc(n:end,2)