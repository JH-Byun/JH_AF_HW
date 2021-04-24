%% IMM_DEMO  Tracking a Target with Simple Manouvers demonstration
clear all
close all
clc

% 
% Simple demonstration for linear IMM using the following models:
%  1. Standard Wiener process velocity model
%  2. Standard Wiener process acceleration model
%
% The measurement model is linear, which gives noisy measurements 
% of target's position.
% 
% Copyright (C) 2007-2008 Jouni Hartikainen
% This software is distributed under the GNU General Public 
%

% (ADDED) add source 
addpath('source');
%% Simple demonstration for IMM using velocity and acceleration models

% Dimensionality of the state space
dims = 4;

% Space for models
ind = cell(1,3);
A = cell(1,3);
Q = cell(1,3);
H = cell(1,3);
R = cell(1,3);

% Stepsize
dt = 0.1;

% Transition matrix for the continous-time velocity model.
F = [0 0 1 0;
      0 0 0 1;
      0 0 0 0;
      0 0 0 0];

% Noise effect matrix for the continous-time system.
L = [0 0;
      0 0;
      1 0;
      0 1];

% model 1
ind{1} = [1 2 3 4]';

% Process noise variance 1
q1 = 0.2^2;
Qc1 = diag([q1 q1]);

% Discretization of the continous-time system. 1
[A{1},Q{1}] = lti_disc(F,L,Qc1,dt);

% Measurement model 1
H{1} = [1 0 0 0;
        0 1 0 0];

% model 2
ind{2} = [1 2 3 4]';

% Process noise variance 2
q2 = 3^2;
Qc2 = diag([q2 q2]);

% Discretization of the continous-time system 2
[A{2},Q{2}] = lti_disc(F,L,Qc2,dt);

% Measurement model 2
H{2} = [1 0 0 0;
        0 1 0 0];
    
% model 3
ind{3} = [1 2 3 4]';

% Process noise variance 3
q3 = 5^2;
Qc3 = diag([q3 q3]);

% Discretization of the continous-time system 2
[A{3},Q{3}] = lti_disc(F,L,Qc3,dt);

% Measurement model 2
H{3} = [1 0 0 0;
        0 1 0 0];
    
hdims = 2;

% Variance in the measurements.
r1 = 1;
r2 = 1;
R{1} = diag([r1 r2]);

r1 = 1;
r2 = 1;
R{2} = diag([r1 r2]);

r1 = 1;
r2 = 1;
R{3} = diag([r1 r2]);

% Generate the data.
n = 200;
Y = zeros(hdims,n);
X_r = zeros(dims,n);
X_r(:,1) = [250 100 10 0]';
mstate = zeros(1,n);

mu_ip{1} = [0.3 0.3 0.4];
% (ADDED) initial distribution
mu_ip{2} = [0.1 0.4 0.5];
mu_ip{3} = [0.1 0.8 0.1];
mu_ip{4} = [0.2 0.1 0.7];
mu_ip{5} = [0.9 0.05 0.05];
mu_0j = mu_ip;

% case 1 : ideal
p_ij = [0.98 0.01 0.01;...
    0.01 0.98 0.01;...
    0.01 0.01 0.98];

% Forced mode transitions 
qstate(1:50) = 0.2;
qstate(51:70) = 1.5;
qstate(71:120) = 3;
qstate(121:150) = 4;
qstate(151:200) = 5;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %
% random seed - do not change      %
rng(15);
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %

% Acceleration model
for i = 2:n
    [A_temp,Q_temp] = lti_disc(F,L,diag([qstate(i)^2 qstate(i)^2]),dt);
    X_r(ind{1},i) = A{1}*X_r(ind{1},i-1) + gauss_rnd(zeros(size(A{1},1),1), Q_temp);
end

% Generate the measurements.
for i = 1:n
    Y(:,i) = H{1}*X_r(ind{1},i) + gauss_rnd(zeros(size(Y,1),1), R{1});
end

%% main algorithm

m = [250 100 10 0]';
P = diag([5 5 1 1]);

%%%% Space for the estimates %%%% 

% Overall estimates of IMM filter
MM = cell(5,n);
PP = cell(5,n);
% Model-conditioned estimates of IMM
MM_i = cell(3,n);
PP_i = cell(3,n);

% Model probabilities 
MU = cell(5,n);

%%%% Prior estimates %%%%

% IMM
x_ip = [250 100 10 0]';
P_ip = diag([5 5 1 1]);

% (ADDED) IMM initialization
for k = 1:5
    MU{k,1} = mu_ip.';
end

% Filtering steps.
for k = 1:5
    MM_i{1,1} = x_ip;
    MM_i{2,1} = x_ip;
    MM_i{3,1} = x_ip;
    PP_i{1,1} = P_ip;
    PP_i{2,1} = P_ip;
    PP_i{3,1} = P_ip;
    c_j = mu_0j{k}.';
    for i = 1:size(Y,2)
        % IMM - 
        if i ~= 1
            MM_temp{1} = MM_i{1,i-1}; MM_temp{2} = MM_i{2,i-1}; MM_temp{3} = MM_i{3,i-1};
            PP_temp{1} = PP_i{1,i-1}; PP_temp{2} = PP_i{2,i-1}; PP_temp{3} = PP_i{3,i-1};
            [MM_temp, PP_temp, c_j, MM{k,i}, PP{k,i}] = imm_predict(MM_temp, PP_temp, MU{k,i-1}, p_ij, ind, dims, A, Q); % (ADDED) IMM prediction
            MM_i{1,i} = MM_temp{1}; MM_i{2,i} = MM_temp{2}; MM_i{3,i} = MM_temp{3};
            PP_i{1,i} = PP_temp{1}; PP_i{2,i} = PP_temp{2}; PP_i{3,i} = PP_temp{3};
        end
        MM_temp{1} = MM_i{1,i};  MM_temp{2} = MM_i{2,i}; MM_temp{3} = MM_i{3,i};
        PP_temp{1} = PP_i{1,i}; PP_temp{2} = PP_i{2,i}; PP_temp{3} = PP_i{3,i};
        [MM_temp, PP_temp, MU{k,i}, MM{k,i}, PP{k,i}] = imm_update(MM_temp, PP_temp, c_j, ind, dims, Y(:,i), H, R); % (ADDED) IMM update
        MM_i{1,i} = MM_temp{1}; MM_i{2,i} = MM_temp{2}; MM_i{3,i} = MM_temp{3};
        PP_i{1,i} = PP_temp{1}; PP_i{2,i} = PP_temp{2}; PP_i{3,i} = PP_temp{3};
    end
end

%% Calculate the MSEs
% (ADDED)
sum{k} = cell(5,1);
for k = 1:5
    sum{k} = zeros(4,1);
end
% (ADDED) get summation of the squares of the estimation errors from each
% time step
for i = 1:size(Y,2)
    for k = 1:5
        sum{k} = sum{k} + (X_r(:,i) - MM{k,i}).^2;
    end
end
% (ADDED) calculate and print MSE (Mean-squared error)
MSE = cell(5,1);
for k = 1:5
    MSE{k} = sum{k}/size(Y,2);
end
fprintf('MSE(IMM_case 1) - x: %f y: %f v_x: %f v_y: %f \n', MSE{1}(1), MSE{1}(2), MSE{1}(3), MSE{1}(4));
fprintf('MSE(IMM_case 2) - x: %f y: %f v_x: %f v_y: %f \n', MSE{2}(1), MSE{2}(2), MSE{2}(3), MSE{2}(4));
fprintf('MSE(IMM_case 3) - x: %f y: %f v_x: %f v_y: %f \n', MSE{3}(1), MSE{3}(2), MSE{3}(3), MSE{3}(4));
fprintf('MSE(IMM_case 4) - x: %f y: %f v_x: %f v_y: %f \n', MSE{4}(1), MSE{4}(2), MSE{4}(3), MSE{4}(4));
fprintf('MSE(IMM_case 5) - x: %f y: %f v_x: %f v_y: %f \n', MSE{5}(1), MSE{5}(2), MSE{5}(3), MSE{5}(4));

%% (ADDED) Plot estimation result
tspan = dt*(1:n); % (ADDED) time span
err_arr = cell(5,1);
MU_arr = cell(3,1);
for k = 1:5
    for i = 1:size(Y,2)
        err_arr{k}(:,i) = X_r(:,i) - MM{k,i};
        MU_arr{k}(:,i) = MU{k,i};
    end
end

figure;
subplot(4,1,1)
plot(tspan, err_arr{1}(1,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{2}(1,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{3}(1,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{4}(1,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{5}(1,:), 'LineWidth', 2)
ylabel('error [m]');
grid on
legend('case 1','case 2','case 3','case 4','case 5');
title('Position x error');
subplot(4,1,2)
plot(tspan, err_arr{1}(2,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{2}(2,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{3}(2,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{4}(2,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{5}(2,:), 'LineWidth', 2)
grid on
ylabel('error [m]');
legend('case 1','case 2','case 3','case 4','case 5');
title('Position y error');
subplot(4,1,3)
plot(tspan, err_arr{1}(3,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{2}(3,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{3}(3,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{4}(3,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{5}(3,:), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
legend('case 1','case 2','case 3','case 4','case 5');
title('Velocity x error');
subplot(4,1,4)
plot(tspan, err_arr{1}(4,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{2}(4,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{3}(4,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{4}(4,:), 'LineWidth', 2)
hold on
plot(tspan, err_arr{5}(4,:), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
xlabel('time [sec]');
legend('case 1','case 2','case 3','case 4','case 5');
title('Velocity y error');

%% (ADDED) Model Probability
figure;
subplot(3,1,1)
plot(tspan, MU_arr{1}(1,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{2}(1,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{3}(1,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{4}(1,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{5}(1,:), 'LineWidth', 2)
grid on
ylabel('\mu_1 [-]');
ylim([0 1]);
legend('case 1','case 2', 'case 3', 'case 4', 'case 5');
title('model 1 probability plot');
subplot(3,1,2)
plot(tspan, MU_arr{1}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{2}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{3}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{4}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{5}(2,:), 'LineWidth', 2)
grid on
ylabel('\mu_2 [-]');
ylim([0 1]);
legend('case 1','case 2', 'case 3', 'case 4', 'case 5');
title('model 2 probability plot');
subplot(3,1,3)
plot(tspan, MU_arr{1}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{2}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{3}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{4}(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU_arr{5}(2,:), 'LineWidth', 2)
grid on
xlabel('time [sec]');
ylabel('\mu_2 [-]');
ylim([0 1]);
legend('case 1','case 2', 'case 3', 'case 4', 'case 5');
title('model 2 probability plot');
