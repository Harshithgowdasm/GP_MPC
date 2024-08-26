
clear all;
close all;
warning('off', 'all')
%% Intialization of States and U input as functions for data 

vdp = VanDerPol_uncertanity();
% load('Van_pol_3');
n = 30;
N = 70; %prediction horizon
% Define initial state and control input
x0 = [1; 1];
u = 1*randn(n+1,1);
noise = randn(n+N+1,1);

% Define time span
t_span = 0:vdp.delta_t:n*vdp.delta_t;

% Initialize array to store state trajectory
x_traj = zeros(2, length(n+1));
y_traj = zeros(2, length(n));
x_traj(:, 1) = x0;

% Simulate the system using the Runge-Kutta method
for i = 2:length(t_span)
    x_traj(:, i) = vdp.f_ud(x_traj(:, i-1), u(i),noise(i));
    dx(:,i-1)= x_traj(:, i)-x_traj(:, i-1);
    y_traj(:, i-1) = x_traj(:, i);
end
figure(5), hold on
plot(x_traj(1,:),x_traj(2,:))
scatter(x_traj(1,:),x_traj(2,:), 'k', 'o', 'SizeData', 10);
grid on
xlabel('x_1')
ylabel('x_2')
title('Phase portrait Van der Pol oscillator')
hold off
%% 

xk = x_traj(:,1:end-1)';
yk= y_traj';
uk= u(2:end,:);
% save('Van_pol_uncertainty1',"xk","yk","uk");

k = 1:1:n+N;

xtrc=[xk,uk]; % aug. xk,uk

%% 
x0 = [1; 1];
utk = 1*randn(N+1,1);

% Define time span
t_span = n*vdp.delta_t:vdp.delta_t:(n+N)*vdp.delta_t;

% Initialize array to store state trajectory
xtk = zeros(2, length(N+1));
ytk = zeros(2, length(N));
xtk(:, 1) = xk(end,:);

% Simulate the system using the Runge-Kutta method
for i = 2:length(t_span)
    xtk(:, i) = vdp.f_ud(xtk(:, i-1),utk(i),noise(n+i));
    dx(:,i-1)= xtk(:, i)-xtk(:, i-1);
    ytk(:, i-1) = xtk(:, i);
end
ytk = ytk';
xkrc = [xtk(:,1:end-1)',utk(2:end,:)];
%% Plotting data points

figure(1),hold on
    title( 'Data points' )
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
set(gca);
disp(' '); disp('plot(xk, yk, ''+'')')
    subplot(ax1)
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title('Displacement')
    ylabel('x_{k+1}(1)','FontSize',12)
    subplot(ax2)
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','#A2142F')
        grid on
        title('Velocity')
        ylabel('x_{k+1}(2)','FontSize',12)
    subplot(ax3), hold on
        stairs(k(:,1:n)',uk(1:n),'Color','k');
        grid on
        title( 'Input');
        xlabel('k','FontSize',12)
        ylabel('U_k','FontSize',12)
        hold off
hold off

cov = {@covSEard}; 
sf = 0.1; 
for i =1:size(xtrc,2)
    ell(i,:) = 10*i;                            % setup the GP
end
hyp0.cov  = log([ell;sf]);
mean = {@meanLinear}; a = 0.9;       % m(x) = a*x+b
hyp0.mean = [a;0.1;a];  % @meanLinear 
% hyp0.mean = [a];  % @meanConst

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 100;                                   % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});
hyp_1 = hyp0; nx =size(xtrc,2)-1;
%% minimization of hyperparametrs 
for i =1:nx
    if Ncg==0
        hyp = hyp0;
    else
        hyp_1 = fitGP(hyp_1, inf, mean, cov, lik,xtrc, yk(:,i),[]);
    end 
    hyp.state(i) = hyp_1;
    [ymu(i), ys(i)] = gp(hyp.state(i), inf, mean, cov, lik, xtrc, yk(:,i),xtrc(end,:));
end
ys(i+1) = 0;

for i=1:N

    if i==1
        mu_p = [ymu,utk(i)]; % µ_tilmda = [µt,ut]
        var_p = diag(ys); % Σ_tilda = blkdiag[Σt, 0] 
        mu_pb = mu_p;
        var_pb = var_p;
        xte = [yk(end,:),utk(i)];
    else
        xte = [mu(i-1,:),utk(i)];
        mu_p = [mu_a(i-1,:),utk(i)];
        var_p = blkdiag(var_a(:,:,i-1),0);
    end
    [mu(i,1), vs(i,1)] = gp(hyp.state(1), inf, mean, cov, lik, xtrc, yk(:,1),xte);
    [mu(i,2), vs(i,2)] = gp(hyp.state(2), inf, mean, cov, lik, xtrc, yk(:,2),xte);

   [mu_a(i,:),var_a(:,:,i)] = Gp_transition_change(mu_p, var_p,hyp,xtrc, yk);
   % [mu_b(i,:),var_b(:,:,i)] = Gp_transition_change_old_old(mu_pb, var_pb,hyp,xtrc, yk);

end

pred = mu_a;

figure(2),hold on
    title( 'GP model' )
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
set(gca);
    subplot(ax1)
        f = [mu(:,1)+2*sqrt(vs(:,1)); flipdim(mu(:,1)-2*sqrt(vs(:,1)),1)];
        fill([k(:,n+1:n+N)'; flipdim(k(:,n+1:n+N)',1)], f, [7 7 7]/8);
        hold on; plot(k(:,n+1:n+N)', mu(:,1),'Color','g');
        plot(k', [yk(1:n,1);ytk(:,1)],'Color','r')
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Displacement')
        ylabel('x_{k+1}(1)','FontSize',12)
        legend('Pred.GP-var','Pred.GP mean','Plant','Data Ponits')
        lgd.FontSize = 5;
    subplot(ax2)
        f = [mu(:,2)+2*sqrt(vs(:,2)); flipdim(mu(:,2)-2*sqrt(vs(:,2)),1)];
        fill([k(:,n+1:n+N)'; flipdim(k(:,n+1:n+N)',1)], f, [7 7 7]/8);
        hold on; plot(k(:,n+1:n+N)', mu(:,2),'Color','g');
        plot(k', [yk(1:n,2);ytk(:,2)],'Color','r')
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Velocity')
        legend('Pred.GP-var','Pred.GP mean','Plant','Data Ponits')
        lgd.FontSize = 5;
        ylabel('x_{k+1}(2)','FontSize',12)
    subplot(ax3), hold on
        stairs([k(:,1:n)';n+1],[uk;utk(1,1)],'Color','k');
        stairs(k(:,n+1:n+N)',utk(1:end-1,:),'Color','b');
        grid on
        title( 'Input');
        legend('Data i/p','Test i/p'),
        ylabel('U_k','FontSize',12)
        xlabel('k','FontSize',12)
        hold off
hold off

for i=1:size(ytk,1)
error(i,:)=(ytk(i)-[mu(i,1),mu(i,2)])*100/ytk(i);
end


RMSE1_1 = sqrt(sum((ytk(:,1)-[mu(:,1)]).^2/N))
RMSE21 = sqrt(sum((ytk(:,2)-[mu(:,2)]).^2/N))

avg11= sum(vs(:,2))/N
avg21= sum(vs(:,2))/N





figure(3),hold on
    title( 'GP_transition model' )
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
set(gca);
    subplot(ax1)

        hold on;  
        plot(k(:,n+1:n+N)', pred(:,1),'Color','g');
        plot(k', [yk(1:n,1);ytk(:,1)],'Color','r')
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r');
        grid on;
        title( 'Displacement');
     
        ylabel('x_{k+1}(1)','FontSize',12);
        legend('Predicted','Plant','Data Ponits');
        lgd.FontSize = 10;
    subplot(ax2)

         plot(k(:,n+1:n+N)', pred(:,2),'Color','g');
        hold on; 
        plot(k', [yk(1:n,2);ytk(:,2)],'Color','r');
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','r');
        grid on;
        title( 'Velocity' );
        legend('Predicted','Plant','Data Ponits');
        lgd.FontSize = 10;
        ylabel('x_{k+1}(2)','FontSize',12);
    subplot(ax3), hold on
        stairs([k(:,1:n)';n+1],[uk;utk(1,1)],'Color','k');
        stairs(k(:,n+1:n+N)',utk(1:end-1,:),'Color','b');
        grid on;
        title( 'Input');
        legend('Data i/p','Test i/p');
        ylabel('U_k','FontSize',12);
        xlabel('k','FontSize',12);
        hold off
hold off



RMSE1 = sqrt(sum((ytk(:,1)-[pred(:,1)]).^2/N))
RMSE2 = sqrt(sum((ytk(:,2)-[pred(:,2)]).^2/N))