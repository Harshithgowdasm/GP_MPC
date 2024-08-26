% Create an instance of the VanDerPol class
clear all;
% close all;
vdp = VanDerPol_uncertanity();

N = 30;

% Define initial state and control input
x0 = [1; 1];
u = randn(N+1,1);
noise = randn(N+1,1);
% Define time span
t_span = 0:vdp.delta_t:N*vdp.delta_t;

% Initialize array to store state trajectory
x_traj = zeros(2, length(N+1));
y_traj = zeros(2, length(N));
x_traj(:, 1) = x0;

% Simulate the system using the Runge-Kutta method
for i = 2:length(t_span)
    x_traj(:, i) = vdp.f_ud(x_traj(:, i-1), u(i),0);  
    dx(:,i-1)= x_traj(:, i)-x_traj(:, i-1);
    y_traj(:, i-1) = x_traj(:, i);
end
figure(6), hold on
plot(x_traj(1,:),x_traj(2,:))
scatter(x_traj(1,1),x_traj(2,1), 'k', 'o', 'SizeData', 10);
% quiver(x_traj(1,1:end-1),x_traj(2,1:end-1),dx,dx,0)
grid on
xlabel('x_1')
ylabel('x_2')
title('Phase portrait Van der Pol oscillator')
hold off
plot3(x_traj(1,:),x_traj(2,:),u)
% xk = x_traj(:,1:end-1)';
% yk= y_traj';
% uk= u(1:end-1,:);
% save('Van_pol_uncertainty6',"xk","yk","uk");
