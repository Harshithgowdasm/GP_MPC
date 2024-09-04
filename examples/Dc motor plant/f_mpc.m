
clear all;
close all;
warning('off', 'all')
%% Intialization of variables
N = 5; %prediction horizon
n = 20; % no of data points
k = 1:1:n;

%% DC motor plant discretization 
sysc = ss(tf(21,[1.1,1,0]));
T = 0.1;
% ra = randn(n,1);
sysd = c2d(sysc,T,"zoh");
[A,B,C]=ctrbf(sysd.A,sysd.B, sysd.C);
f = @(x1,x2,u) A(1,:) *[x1;x2] + B(1,:)*u; %+ 0.02*randn(1);
f1 = @(x1,x2,u) A(2,:) *[x1;x2] + B(2,:)*u;

%% Parameter Uncertainty to plant model
% A_d = A-0.1*rand(size(A))*A;
% B_d = B-0.1*rand(size(B)).*B;
% f = @(x1,x2,u,ra) A_d(1,:) *[x1;x2] + B_d(1,:)*u+ 0.005*ra;
% f1 = @(x1,x2,u,ra) A_d(2,:) *[x1;x2] + B_d(2,:)*u+ 0.005 *ra;

%% Generating data-States and U input 
% uk = randn(n,1); 
% x1 = 1; x2 = pi;
% for i=0:n-1
%     yk(i+1,1) = f(x1(i+1),x2(i+1),uk(i+1),ra(i+1,1));
%     yk(i+1,2) = f1(x1(i+1),x2(i+1),uk(i+1),ra(i+1,1));
%     x1(i+2) = yk(i+1,1);
%     x2(i+2) = yk(i+1,2);
% 
% end
% yk2 = yk;
% xk = [x1(1:n)',x2(1:n)'];
% save('data_dc_motor_6',"xk","yk2","uk")

%% load generated data (if want to load the data)
load("data_dc_motor_noise.mat");
yk = yk2;

%% augmentation of xk and uk
xtrc=[xk,uk]; % aug. xk,uk

%% Plotting data points

figure(3)
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
        title( 'Angle','FontSize',12)
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
Q=[805,0;0,10]*eye(2); %data_dc_motor_2
R=1*eye(1);
%% f-MPC learning
ra=randn(3*N,1);
tic;
for i=1:15
    display(i,'iteration');
    u0= u0*ones(1,N); %%intial input u
    fun = @(u0) costf(u0,yk,Q,R,f,f1,ra(i:i+N-1));    % define cost function
   
    %% minimization of cost function 
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    [U,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,u0,A,b,Aeq,beq,lb,ub,[],options);
    u0=U(1)-eps;            % optimal control input u first element 
    xtrc(size(xtrc,1)+1,:) = [yk(end,:),U(1)];      % update data set xk & uk
    yk(size(yk,1)+1,:) = [f(yk(end,1),yk(end,2),U(1)),f1(yk(end,1),yk(end,2),U(1))];    % update data set yk
    
    %% update plot
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

