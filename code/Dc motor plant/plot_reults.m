
clear all;
close all;
warning('off', 'all')
%% DC motor plant discretization and Intialization of variables and Generating data set
sysc = ss(tf(21,[1.1,1,0]));
T = 0.1;
sysd = c2d(sysc,T,"zoh");
[A,B,C]=ctrbf(sysd.A,sysd.B, sysd.C);
A_d = A-0.1*rand(size(A))*A;
B_d = B-0.1*rand(size(B)).*B;
f = @(x1,x2,u) A_d(1,:) *[x1;x2] + B_d(1,:)*u + 0.005*randn(1);
f1 = @(x1,x2,u) A_d(2,:) *[x1;x2] + B_d(2,:)*u + 0.005*randn(1);


N = 60; %prediction horizon
n = 20; % no of data points
k = 1:1:n+N;
% uk = randn(n,1);
% utk = randn(N,1);
load('data_dc_Gp_test_final');
x1 = 0.1; x2 = pi;
for i=0:n-1
    yk(i+1,1) = f(x1(i+1),x2(i+1),uk(i+1));
    yk(i+1,2) = f1(x1(i+1),x2(i+1),uk(i+1));
    x1(i+2) = yk(i+1,1);
    x2(i+2) = yk(i+1,2);

end
yk2 = yk;
xk = [x1(1:n)',x2(1:n)'];
% save('data_dc_Gp_test_final',"uk","utk")
% load("data_dc_motor_9.mat");
% yk = yk2;
xtrc=[xk,uk]; % aug. xk,uk

%% generating test set 
x1 = yk(end,1);
x2 = yk(end,2);
for i=0:N-1
    ytk(i+1,1) = f(x1(i+1),x2(i+1),utk(i+1));
    ytk(i+1,2) = f1(x1(i+1),x2(i+1),utk(i+1));
    x1(i+2) = ytk(i+1,1);
    x2(i+2) = ytk(i+1,2);

end
xtk = [x1(1:N)',x2(1:N)'];
xkrc = [xtk,utk];
%% Plotting data points

figure(1),hold on
    title( 'Data points' )
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
set(gca)
disp(' '); disp('plot(xk, yk, ''+'')')
    subplot(ax1)
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Angular Velocity','FontSize',14 )
    ylabel('x_{k+1}(1)','FontSize',12)
    subplot(ax2)
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','#A2142F')
        grid on
        title( 'Angle','FontSize',14 )
        ylabel('x_{k+1}(2)','FontSize',12)
    subplot(ax3), hold on
        stairs(k(:,1:n)',uk(1:n),'Color','k');
        grid on
        title( 'Input','FontSize',14 );
        xlabel('k','FontSize',12)
        ylabel('U_k','FontSize',12)
        hold off
hold off

%% define the GP hyperparametrs & mean, cov, lik functions

cov = {@covSEard}; 
sf = 0.1; 
for i =1:size(xtrc,2)
    ell(i,:) = 3;                            % setup the GP
end
hyp0.cov  = log([ell;sf]);
mean = {@meanLinear}; a = 0.1;       % m(x) = a*x+b
hyp0.mean = [a;0.1;a];  % @meanLinear 
% hyp0.mean = [a];  % @meanConst

lik = {'likGauss'};    % Gauss likelihoods
inf = {'infGaussLik'}; % inference algs

Ncg = 100;                                   % number of conjugate gradient steps
                  
sn = 0.2;  
hyp0.lik  = [log(sn)];
fprintf('OPT: %s/%s\n',lik{1},inf{1});
hyp_1 = hyp0; 
nx =size(xtrc,2)-1;
%% optimization of hyperparametrs 
for i =1:nx
    if Ncg==0
        hyp = hyp0;
    else
        hyp_1 = fitGP(hyp_1, inf, mean, cov, lik,xtrc, yk(:,i),[]);
    end 
    hyp.state(i) = hyp_1;
    [ymu(i), ys(i)] = gp(hyp.state(i), inf, mean, cov, lik, xtrc, yk(:,i),xtrc(end,:));
end
ys(i+1) = 0;        % set variance of input u as zero

%% Testing GP and Gp transition model(moment matching)

for i=1:N

    if i==1
        mu_p = [ymu,utk(i)]; % µ_tilmda = [µt,ut]
        var_p = diag(ys); % Σ_tilda = blkdiag[Σt, 0] 
        mu_pb = mu_p;
        var_pb = var_p;
        xte = [ymu,utk(i)];
    else
        xte = [mu(i-1,:),utk(i)];
        mu_p = [mu_a(i-1,:),utk(i)];
        var_p = blkdiag(var_a(:,:,i-1),0);
    end     
    [mu(i,1), vs(i,1)] = gp(hyp.state(1), inf, mean, cov, lik, xtrc, yk(:,1),xte);
    [mu(i,2), vs(i,2)] = gp(hyp.state(2), inf, mean, cov, lik, xtrc, yk(:,2),xte);
   % [mu_a(i,:),var_a(:,:,i)] = Gp_transition_change(mu_p, var_p,hyp,xtrc, yk);
   [mu_a(i,:),var_a(:,:,i)] = Gp_transition_change_2_old(mu_p, var_p,hyp,xtrc, yk);

end

pred = mu_a;
%% plot GP testing results

figure(2),hold on
    title( 'GP model' )
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
set(gca)
    subplot(ax1)
        f = [mu(:,1)+2*sqrt(vs(:,1)); flipdim(mu(:,1)-2*sqrt(vs(:,1)),1)];
        fill([k(:,n+1:n+N)'; flipdim(k(:,n+1:n+N)',1)], f, [7 7 7]/8);
        hold on; plot(k(:,n+1:n+N)', mu(:,1),'Color','g');
        plot(k', [yk(1:n,1);ytk(:,1)],'Color','r')
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Angular Velocity')
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
        title( 'Angle')
        legend('Pred.GP-var','Pred.GP mean','Plant','Data Ponits')
        lgd.FontSize = 5;
        ylabel('x_{k+1}(2)','FontSize',12)
    subplot(ax3), hold on
        stairs([k(:,1:n)';n+1],[uk;utk(1,1)],'Color','k');
        stairs(k(:,n+1:n+N)',utk(1:end,:),'Color','b');
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


RMSE11 = sqrt(sum((ytk(:,1)-[mu(:,1)]).^2/N));
RMSE21 = sqrt(sum((ytk(:,2)-[mu(:,2)]).^2/N));

avg11= sum(vs(:,2))/N;
avg21= sum(vs(:,2))/N;

%% plot Moment matching testing results

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
        plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Angular Velocity')
        ylabel('x_{k+1}(1)','FontSize',12)
        legend('Predicted','Plant','Data Ponits')
        lgd.FontSize = 10;
    subplot(ax2)
         plot(k(:,n+1:n+N)', pred(:,2),'Color','g');
        hold on; 
        plot(k', [yk(1:n,2);ytk(:,2)],'Color','r')
        plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12,'Color','r')
        grid on
        title( 'Angle' )
        legend('Predicted','Plant','Data Ponits')
        lgd.FontSize = 10;
        ylabel('x_{k+1}(2)','FontSize',12)
    subplot(ax3), hold on
        stairs([k(:,1:n)';n+1],[uk;utk(1,1)],'Color','k');
        stairs(k(:,n+1:n+N)',utk(1:end,:),'Color','b');
        grid on
        title( 'Input');
        legend('Data i/p','Test i/p')
        ylabel('U_k','FontSize',12)
        xlabel('k','FontSize',12)
        hold off
hold off



RMSE1 = sqrt(sum((ytk(:,1)-[pred(:,1)]).^2)/N);
RMSE2 = sqrt(sum((ytk(:,2)-[pred(:,2)]).^2)/N);