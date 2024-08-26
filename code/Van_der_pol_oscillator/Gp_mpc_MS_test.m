
clear all;
close all;
warning('off', 'all')
vdp = VanDerPol_uncertanity();
%% Intialization of variables and importing data points  

N = 15; %prediction horizon
n = 30; % no of data points
k = 1:1:n;
load("Van_pol_uncertainty4.mat");
% load('Van_pol_3');

% yk = yk2;
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
        title( 'Displacement','FontSize',12)
        % xlabel('k')
        ylabel('x_{k+1}(1)','FontSize',12)
        grid on
    subplot(ax2)
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','#A2142F')
        grid on
        title( 'Velocity','FontSize',12 )
        
        ylabel('x_{k+1}(2)','FontSize',12)
        grid on
    subplot(ax3), hold on
        stairs(k(:,1:n)',uk(1:n),'Color','k');
        grid on
        title( 'Input','FontSize',12 );
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
        title( 'Displacement','FontSize',12)
        ylabel('x_{k+1}(1)','FontSize',12)
        grid on
    subplot(ax2)
        plot(k(:,1:n)',yk(1:n,2),'Color','#A2142F');
        title( 'Velocity','FontSize',12)
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
Q=[4000,0;0,5500]*eye(2); %data_dc_motor_9 Q=[105,0;0,10]*eye(2)
R=1*eye(1); %data_dc_motor_9 R=1*eye(1)

noise = 0.05*randn(20,1);
%% GP-MPC learning 
tic;
for i=1:20
    display(i,'iteration');
    [U] = Gp_mpc(xtrc,yk,u0,N,lb,ub,Q,R);
    u0=U(1)-eps
    xtrc(size(xtrc,1)+1,:) = [yk(end,:),U(1)];
    yk(size(yk,1)+1,:) = vdp.f_ud(yk(end,:)',U(1),noise(i))';
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
elapsed_time = toc;

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
