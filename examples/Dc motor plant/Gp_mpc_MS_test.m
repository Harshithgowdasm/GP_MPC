% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda

clear all;
close all;
warning('off', 'all')
%% DC motor plant discretization 

sysc = ss(tf(21,[1.1,1,0]));
T = 0.1;
sysd = c2d(sysc,T,"zoh");
[A,B,C]=ctrbf(sysd.A,sysd.B, sysd.C);
f = @(x1,x2,u) A(1,:) *[x1;x2] + B(1,:)*u; %+ 0.02*randn(1);
f1 = @(x1,x2,u) A(2,:) *[x1;x2] + B(2,:)*u; %+ 0.02*randn(1);

%% Parameter Uncertainty to plant model

% A_d = A-0.1*rand(size(A))*A;
% B_d = B-0.1*rand(size(B)).*B;
% f = @(x1,x2,u) A_d(1,:) *[x1;x2] + B_d(1,:)*u+ 0.005*randn(1);
% f1 = @(x1,x2,u) A_d(2,:) *[x1;x2] + B_d(2,:)*u+ 0.005 *randn(1);
% load('data_dc_Gp_test_final');

%% Intialization of variables 

N = 15; %prediction horizon
n = 20; % no of data points
k = 1:1:n;

%% Generating data-States and U input 
% uk = randn(n,1); 
% x1 = 0; x2 = pi;
% for i=0:n-1
%     yk(i+1,1) = f(x1(i+1),x2(i+1),uk(i+1));
%     yk(i+1,2) = f1(x1(i+1),x2(i+1),uk(i+1));
%     x1(i+2) = yk(i+1,1);
%     x2(i+2) = yk(i+1,2);
% 
% end
% yk2 = yk;
% xk = [x1(1:n)',x2(1:n)'];
% save('data_dc_motor_wo_noise',"xk","yk2","uk")

%% load generated data (if want to load the data)
load("data_dc_motor_wo_noise.mat"); 
yk = yk2;

%% augmentation of xk and uk
xtrc=[xk,uk]; % aug. xk,uk

%% Plotting data points

figure(1),hold on
    title( 'Data points' )
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
set(gca)
    subplot(ax1)
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Angular Velocity','FontSize',14)
        % xlabel('k')
        ylabel('x_{k+1}(1)','FontSize',12)
        grid on
    subplot(ax2)
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','#A2142F')
        grid on
        title( 'angle' )
        % xlabel('k')
        title( 'Angle','FontSize',14)
        ylabel('x_{k+1}(2)','FontSize',12)
        grid on
    subplot(ax3), hold on
        stairs(k(:,1:n)',uk(1:n),'Color','k');
        grid on
        title( 'Input','FontSize',14 );
        xlabel('k','FontSize',12)
        ylabel('U_k','FontSize',12)
        grid on
        hold off
hold off

figure(2)
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
    subplot(ax1)
        plot(k(:,1:n)',yk(1:n,1),'Color','r'); 
        title( 'Angular Velocity','FontSize',12)
        ylabel('x_{k+1}(1)','FontSize',12)
        grid on
    subplot(ax2)
        plot(k(:,1:n)',yk(1:n,2),'Color','#A2142F');
        title( 'Angle','FontSize',14)
        ylabel('x_{k+1}(2)','FontSize',12)
        grid on
    subplot(ax3), hold on
        stairs(k(:,1:n)',uk(1:n),'Color','k');
        title( 'Input','FontSize',12 );
        xlabel('k','FontSize',12)
        ylabel('U_k','FontSize',12)
        grid on
        hold off


%% Input constraints and Weighing matrices Q & R
u0=1*uk(end);
lb = -10*ones(1,N);
ub = 10*ones(1,N);

Q=[805,0;0,10]*eye(2); %data_dc_motor_9 Q=[105,0;0,10]*eye(2)
R=1*eye(1); %data_dc_motor_9 R=1*eye(1)


%% GP-MPC learning 
tic;
for i=1:15
    [U] = Gp_mpc(xtrc,yk,u0,N,lb,ub,Q,R);
    u0=U(1)-eps;                     % optimal control input u first element  
    xtrc(size(xtrc,1)+1,:) = [yk(end,:),U(1)];     % update data set xk & uk
    yk(size(yk,1)+1,:) = [f(yk(end,1),yk(end,2),U(1)),f1(yk(end,1),yk(end,2),U(1))];  % update data set yk

    %% update plot 
    k = 1:1:size(xtrc,1);
    figure(2),hold on
    subplot(ax1), hold on
    plot(k(:,n:end)',yk(n:end,1),'Color','g');
    hold off
    subplot(ax2), hold on
    plot(k(:,n:end)',yk(n:end,2),'Color','g');
    hold off
    subplot(ax3), hold on
    stairs(k(:,n:end)',xtrc(n:end,end),'Color','b')
    hold off
end

%% dispaly the results 

projected_yk = yk(n:end,:)
optimal_uk = xtrc(n:end,end)
elapsed_time = toc


figure(2),hold on
    subplot(ax1), hold on
    legend('Plant Data','GP-MPC')
    hold off
    subplot(ax2), hold on
    legend('Plant Data','GP-MPC')
    hold off
    subplot(ax3), hold on
    legend('Plant Data','GP-MPC')
    hold off
