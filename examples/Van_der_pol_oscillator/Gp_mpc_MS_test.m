% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Harshith Gowda

clear all;  % Clear all variables from the workspace
close all;  % Close all figure windows

% Add the necessary paths to access the GP function (commented out here)
% addpath('directory\PA_Harshith_Gowda\GP');

warning('off', 'all')  % Disable all warnings
vdp = VanDerPol_uncertanity_dynamics();  % Initialize the Van der Pol dynamics with uncertainty

%% Initialization of variables and importing data points  
N = 5;  % Prediction horizon for MPC
n = 30; % Number of initial data points
k = 1:1:n;  % Time index for plotting

% Load data for Van der Pol oscillator with uncertainty
load("data/Van_pol_uncertainty4.mat");

% Augment the state and control input for training
xtrc = [xk, uk];  % Concatenated state and control input

%% Plotting data points
% Plot the initial data points for displacement, velocity, and input

figure(1), hold on
title('Data points')
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);

% Plot the displacement data points
subplot(ax1)
plot(k(:,1:n)', yk(1:n,1), '+', 'MarkerSize', 12, 'Color', 'r')
grid on
title('Displacement', 'FontSize', 12)
ylabel('x_{k+1}(1)', 'FontSize', 12)

% Plot the velocity data points
subplot(ax2)
plot(k(:,1:n)', yk(1:n,2), '+', 'MarkerSize', 12, 'Color', '#A2142F')
grid on
title('Velocity', 'FontSize', 12)
ylabel('x_{k+1}(2)', 'FontSize', 12)

% Plot the input data points
subplot(ax3), hold on
stairs(k(:,1:n)', uk(1:n), 'Color', 'k');
grid on
title('Input', 'FontSize', 12)
xlabel('k', 'FontSize', 12)
ylabel('U_k', 'FontSize', 12)
hold off
hold off

% Plot the continuous trajectory of the system's displacement, velocity, and input
figure(2)
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);

subplot(ax1)
plot(k(:,1:n)', yk(1:n,1), 'Color', 'r'); 
title('Displacement', 'FontSize', 12)
ylabel('x_{k+1}(1)', 'FontSize', 12)
grid on

subplot(ax2)
plot(k(:,1:n)', yk(1:n,2), 'Color', '#A2142F');
title('Velocity', 'FontSize', 12)
ylabel('x_{k+1}(2)', 'FontSize', 12)
grid on

subplot(ax3), hold on
stairs(k(:,1:n)', uk(1:n), 'Color', 'k');
title('Input', 'FontSize', 12);
xlabel('k', 'FontSize', 12)
ylabel('U_k', 'FontSize', 12)
grid on
hold off

%% Input constraints and Weighing matrices Q & R
% Initialize control input and constraints for MPC
u0 = 1 * uk(end);  % Initial control input
lb = -10 * ones(1, N);  % Lower bound on control input
ub = 10 * ones(1, N);   % Upper bound on control input

% Weighing matrices for the MPC cost function
Q = [4000, 0; 0, 5500] * eye(2);  % State weighting matrix
R = 1 * eye(1);  % Control input weighting matrix

noise = 0.05 * randn(20, 1);  % Additive noise for uncertainty

%% GP-MPC learning 
tic;  % Start timing the execution
for i = 1:20
    display(i, 'iteration');  % Display the current iteration
    [U] = Gp_mpc(xtrc, yk, u0, N, lb, ub, Q, R);  % Solve the GP-MPC problem

    u0 = U(1) - eps;  % Update the control input for the next iteration
    xtrc(size(xtrc, 1) + 1, :) = [yk(end, :), U(1)];  % Update the state-control history
    yk(size(yk, 1) + 1, :) = vdp.f_ud(yk(end, :)', U(1), noise(i))';  % Simulate the system with the new control input

    %% Update plot 
    k = 1:1:size(xtrc, 1);  % Update the time index
    figure(2), hold on
    
    % Update the displacement plot with GP-MPC results
    subplot(ax1), hold on
    plot(k(:,n:end)', yk(n:end, 1), 'Color', 'g');
    hold off
    
    % Update the velocity plot with GP-MPC results
    subplot(ax2), hold on
    plot(k(:,n:end)', yk(n:end, 2), 'Color', 'g');
    hold off
    
    % Update the input plot with GP-MPC results
    subplot(ax3), hold on
    stairs(k(:,n:end)', xtrc(n:end, end), 'Color', 'b')
    hold off
end

%% Display the results 
projected_yk = yk(n:end, :);  % Extract the predicted states
optimal_uk = xtrc(n:end, end);  % Extract the optimal control inputs
elapsed_time = toc;  % Measure the elapsed time for the entire process

% Finalize the plot by adding legends
figure(2), hold on
subplot(ax1), hold on
legend('Plant Data', 'GP-MPC')
hold off

subplot(ax2), hold on
legend('Plant Data', 'GP-MPC')
hold off

subplot(ax3), hold on
legend('Plant Data', 'GP-MPC')
hold off
