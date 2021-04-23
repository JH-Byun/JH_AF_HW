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
dims = 6;

% Space for models
ind = cell(1,2);
A = cell(1,2);
Q = cell(1,2);
H = cell(1,2);
R = cell(1,2);

% Stepsize
dt = 0.1;

ind{2} = [1 2 3 4 5 6]';

% Transition matrix for the continous-time acceleration model.
F2 = [0 0 1 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 0 0 0;
      0 0 0 0 0 0];

% Noise effect matrix for the continous-time system.
L2 = [0 0;
      0 0;
      0 0;
      0 0;
      1 0;
      0 1];

% Process noise variance
q2 = 1.00;
Qc2 = diag([q2 q2]);

% Discretization of the continous-time system.
[A{2},Q{2}] = lti_disc(F2,L2,Qc2,dt);

ind{1} = [1 2 3 4]';

% Transition matrix for the continous-time velocity model.
F1 = [0 0 1 0;
      0 0 0 1;
      0 0 0 0;
      0 0 0 0];

% Noise effect matrix for the continous-time system.
L1 = [0 0;
      0 0;
      1 0;
      0 1];

  % Process noise variance
q1 = .01;
Qc1 = diag([q1 q1]);

% Discretization of the continous-time system.
[A{1},Q{1}] = lti_disc(F1,L1,Qc1,dt);

H{1} = [1 0 0 0;
        0 1 0 0];

% Measurement models.
H{2} = [1 0 0 0 0 0;
        0 1 0 0 0 0];
hdims = 2;

% Variance in the measurements.
r1 = .1;
r2 = .1;
R{1} = diag([r1 r2]);

r1 = .1;
r2 = .1;
R{2} = diag([r1 r2]);

% Generate the data.
n = 200;
Y = zeros(hdims,n);
X_r = zeros(dims,n);
X_r(:,1) = [0 0 0 -1 0 0]';
mstate = zeros(1,n);

mu_ip = [0.9 0.1];
mu_0j = mu_ip;

% case 1 : ideal
p_ij{1} = [0.98 0.02;
        0.02 0.98];

% (ADDED) case 2 ~ 5    
p_ij{2} = [0.9 0.1;...
    0.1 0.9];
p_ij{3} = [0.8 0.2;...
    0.2 0.8];
p_ij{4} = [0.75 0.25;...
    0.25 0.75];
p_ij{3} = [0.75 0.25;...
    0.25 0.75];

% Forced mode transitions 
mstate(1:50) = 2;
mstate(51:70) = 2;
mstate(71:120) = 2;
mstate(121:150) = 1;
mstate(151:200) = 1;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %
% random seed - do not change      %
rng(15);
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %

% Acceleration model
for i = 2:n
    st = mstate(i);
    X_r(ind{st},i) = A{st}*X_r(ind{st},i-1) + gauss_rnd(zeros(size(A{st},1),1), Q{st});
end

% Generate the measurements.
for i = 1:n
    Y(:,i) = H{mstate(i)}*X_r(ind{mstate(i)},i) + gauss_rnd(zeros(size(Y,1),1), R{mstate(i)});
end


%% main algorithm

m = [0 0 0 -1 0 0]';
P = diag([0.1^2 0.1^2 0.1^2 0.1^2 0.5^2 0.5^2]);

%%%% Space for the estimates %%%% 

% Overall estimates of IMM filter
MM = cell(5,1);
PP = cell(5,1);
for k = 1:5
    MM{k} = zeros(size(m,1), size(Y,2));
    PP{k} = zeros(size(m,1), size(m,1), size(Y,2));
end

% Model-conditioned estimates of IMM
MM_i = cell(2,n);
PP_i = cell(2,n);

% Model probabilities 
MU = cell(5,1);
for k = 1:5
    MU{k} = zeros(2,size(Y,2));
end


%%%% Prior estimates %%%%

% IMM
x_ip{1} = [0 0 0 -1]';
P_ip{1} = diag([0.1^2 0.1^2 0.1^2 0.1^2]);
x_ip{2} = [0 0 0 -1 0 0]';
P_ip{2} = diag([0.1^2 0.1^2 0.1^2 0.1^2 0.5^2 0.5^2]);

% (ADDED) IMM initialization
MM_i{1,1} = x_ip{1};
MM_i{2,1} = x_ip{2};
PP_i{1,1} = P_ip{1};
PP_i{2,1} = P_ip{2};
for k = 1:5
    MU{k}(:,1) = mu_ip.';
end
c_j = mu_0j.';

% Filtering steps.
for i = 1:size(Y,2)
    for k = 1:5
        % case k
        if i ~= 1
            MM_temp{1} = MM_i{1,i-1};  MM_temp{2} = MM_i{2,i-1};
            PP_temp{1} = PP_i{1,i-1}; PP_temp{2} = PP_i{2,i-1};
            [MM_temp, PP_temp, c_j, MM{k}(:,i), PP{k}(:,:,i)] = imm_predict(MM_temp, PP_temp, MU(:,i-1), p_ij{k}, ind, dims, A, Q); % (ADDED) IMM prediction
            MM_i{1,i} = MM_temp{1}; MM_i{2,i} = MM_temp{2};
            PP_i{1,i} = PP_temp{1}; PP_i{2,i} = PP_temp{2};
        end
        MM_temp{1} = MM_i{1,i};  MM_temp{2} = MM_i{2,i};
        PP_temp{1} = PP_i{1,i}; PP_temp{2} = PP_i{2,i};
        [MM_temp, PP_temp, MU(:,i), MM{k}(:,i), PP{k}(:,:,i)] = imm_update(MM_temp, PP_temp, c_j, ind, dims, Y(:,i), H, R); % (ADDED) IMM update
        MM_i{1,i} = MM_temp{1}; MM_i{2,i} = MM_temp{2};
        PP_i{1,i} = PP_temp{1}; PP_i{2,i} = PP_temp{2};
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
    for j = 1:4
        for k = 1:5
            sum{k}(j,1) = sum{k}(j,1) + (X_r(j,i) - MM{k}(j,i))^2;
        end
    end
end
% (ADDED) calculate and print MSE (Mean-squared error)
MSE = cell(5,1);
for k = 1:5
    MSE{k} = sum{k}(j,1)/size(Y,2);
end
fprintf('MSE(IMM_case 1) - x: %f y: %f v_x: %f v_y: %f \n', MSE{1}(1), MSE{1}(2), MSE{1}(3), MSE{1}(4));
fprintf('MSE(IMM_case 2) - x: %f y: %f v_x: %f v_y: %f \n', MSE{2}(1), MSE{2}(2), MSE{2}(3), MSE{2}(4));
fprintf('MSE(IMM_case 3) - x: %f y: %f v_x: %f v_y: %f \n', MSE{3}(1), MSE{3}(2), MSE{3}(3), MSE{3}(4));
fprintf('MSE(IMM_case 4) - x: %f y: %f v_x: %f v_y: %f \n', MSE{4}(1), MSE{4}(2), MSE{4}(3), MSE{4}(4));
fprintf('MSE(IMM_case 5) - x: %f y: %f v_x: %f v_y: %f \n', MSE{5}(1), MSE{5}(2), MSE{5}(3), MSE{5}(4));

%% (ADDED) Plot estimation result
tspan = dt*(1:n); % (ADDED) time span
figure;
subplot(4,1,1)
plot(tspan, X_r(1,:) - MM1(1,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(1,:) - MM2(1,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(1,:) - MMMMAE(1,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(1,:) - MM(1,1:n), 'LineWidth', 2)
ylabel('error [m]');
grid on
legend('KF1','KF2','MMAE','IMM');
title('Position x error');
subplot(4,1,2)
plot(tspan, X_r(2,:) - MM1(2,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(2,:) - MM2(2,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(2,:) - MMMMAE(2,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(2,:) - MM(2,1:n), 'LineWidth', 2)
grid on
ylabel('error [m]');
legend('KF1','KF2','MMAE','IMM');
title('Position y error');
subplot(4,1,3)
plot(tspan, X_r(3,:) - MM1(3,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(3,:) - MM2(3,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(3,:) - MMMMAE(3,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(3,:) - MM(3,1:n), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
legend('KF1','KF2','MMAE','IMM');
title('Velocity x error');
subplot(4,1,4)
plot(tspan, X_r(4,:) - MM1(4,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(4,:) - MM2(4,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(4,:) - MMMMAE(4,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(4,:) - MM(4,1:n), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
xlabel('time [sec]');
legend('KF1','KF2','MMAE','IMM');
title('Velocity y error');

%% (ADDED) Model Probability
figure;
subplot(2,1,1)
plot(tspan, MUMMAE(1,:), 'LineWidth', 2)
hold on
plot(tspan, MUMMAE(2,:), 'LineWidth', 2)
grid on
ylabel('\mu [-]');
legend('model 1','model 2');
title('MMAE model probability plot');
subplot(2,1,2)
plot(tspan, MU(1,:), 'LineWidth', 2)
hold on
plot(tspan, MU(2,:), 'LineWidth', 2)
grid on
ylabel('\mu [-]');
xlabel('time [sec]');
legend('model 1','model 2');
title('IMM model probability plot');