
clear all;
close all;
warning('off', 'all')
vdp = VanDerPol_uncertanity();
%% Intialization of variables
N = 5; %prediction horizon
n = 30; % no of data points
k = 1:1:n;
load("Van_pol_uncertainty5.mat");
% load('Van_pol_3');

% yk = yk2;
xtrc=[xk,uk]; % aug. xk,uk
%% Plotting data points

figure(3)
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

Q=[6000,0;0,5500]*eye(2); 
R=1*eye(1); 
%% f-MPC
ra=0*randn(n,1);
tic;
% ra(i:i+5,1)
for i=1:20
    display(i,'iteration');
    u0= u0*ones(1,N); %%intial input u
    fun = @(u0) costf(u0,yk,Q,R,ra(i:i+5,1));
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

    A = [];
    b = []; 
    Aeq = [];
    beq = [];

    [U,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,u0,A,b,Aeq,beq,lb,ub,[],options)
    display('-----------------------------------------')
    u0=U(1)-eps
    xtrc(size(xtrc,1)+1,:) = [yk(end,:),U(1)];
    yk(size(yk,1)+1,:) = vdp.f_ud(yk(end,:)',U(1),ra(i,1))';
    display(u0,'Optimal U');
    display(yk(end,:),'yk');
    display(xk(end,:),'xk');
    display('---------------------------');
    k = 1:1:size(xtrc,1);
    figure(3),hold on
    subplot(ax1), hold on
    plot(k(:,n:end)',yk(n:end,1),'Color','g');
    hold off
    subplot(ax2), hold on
    plot(k(:,n:end)',yk(n:end,2),'Color','g');
    hold off
    subplot(ax3), hold on
    stairs(k(:,n:end)',xtrc(n:end,end),'Color','b');
    hold off
end
figure(3),hold on
    subplot(ax1), hold on
    legend('Plant Data','f-MPC')
    hold off
    subplot(ax2), hold on
    legend('Plant Data','f-MPC')
    hold off
    subplot(ax3), hold on
    legend('Plant Data','f-MPC')
    hold off

projected_yk = yk(n:end,:)
optimal_uk = xtrc(n:end,end)
elapsed_time = toc;

