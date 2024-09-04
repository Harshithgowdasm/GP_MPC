% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda

clear all;
close all;

%% Intialization of States and U input as functions for data 
% f = @(t) sin(1*t)+0.00002*randn(size(t));
% f = @(t) exp(0.08*t)+0.00002*randn(size(t));
% g= @(t) exp(0.01*t)+0.00002*randn(size(t));
% f1 = @(t) exp(0.1*t)+0.00002*randn(size(t));
f = @(x,uk) 0.2*x+0.7*uk;
f1 = @(x,uk) 0.7*x+0.2*uk;
% g = @(t) 1* t+0.00002*randn(size(t));
N = 10; %prediction horizon
n = 40; % no of data points
k = 1:1:n;
uk = randn(n,1);
% uk = 1*ones(1,n);
% uk = 0.05:0.05:n*0.05;
x1 = 1; x2 = 1;
for i=0:n-1
    yk(i+1,1) = f(x1(i+1),uk(i+1));
    % yk(i+1,2) = f1(x2(i+1),uk(i+1));
    x1(i+2) = yk(i+1,1);
    % x2(i+2) = yk(i+1,2);

end
xk = [x1(1:n)'];
% xk = [x1(1:n)',x2(1:n)'];
% uk = uk';
xtrc=[xk,uk]; % aug. xk,uk

%% Plotting data points

figure(1), hold on
plot(k,yk,'Color','r','LineWidth',2)
plot(k,uk,'Color','k','LineWidth',2)
hold off
%% input constraints 
u0=1*uk(end);
lb = -10*ones(1,N);
ub = 10*ones(1,N);

%% GP-MPC
for i=1:10
    display(i,'iteration');
    [U] = Gp_mpc(xtrc,yk,u0,N,lb,ub);
    u0=U(1);
    xtrc(size(xtrc,1)+1,:) = [yk(end,:),U(1)];
    % yk(size(yk,1)+1,:) = [f(yk(end,1),U(1)),f1(yk(end,2),U(1))];
    yk(size(yk,1)+1,:) = [f(yk(end,1),U(1))];
    display(u0,'Optimal U');
    display(yk(end,:),'yk');
    display(xk(end,:),'xk')
    display('---------------------------');


    k = n:1:size(xtrc,1);
    figure(1), hold on
    plot(k,yk(n:end,:),'Color','g','LineWidth',1)
    plot(k,xtrc(n:end,end),'Color','b','LineWidth',1)
    % legend('yk1','yk2','uk','opt. Yk1','opt. Yk2','opt. Uk')
    hold off

end
projected_yk = yk(n:end,:)
optimal_uk = xtrc(n:end,end)



% xu=yk(end);
% 
% for j=1:N
%     yuk(j,1) = f(xu(j),U(j));
%     xu(j+1) = yuk(j,1);
% end
% figure(2), hold on
% fprintf('%g ',yuk)
% plot(n+1:1:n+1*N,yuk,'Color','y')
% plot(n+1:1:n+1*N,U,'Color','k')
% 
% % legend('yk','uk','opt. Yk','opt. Uk')
% hold off

