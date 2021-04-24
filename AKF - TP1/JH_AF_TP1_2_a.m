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

mu_ip = [0.3 0.3 0.4];
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

% (ADDED) Overall estimates of MMAE filter
MMMMAE = zeros(size(m,1), size(Y,2));
PPMMAE = zeros(size(m,1), size(m,1), size(Y,2));
% (ADDED) Model-conditioned predictions of MMAE
MMMMAE_i = cell(3,n);
PPMMAE_i = cell(3,n);

% (ADDED) Model probabilities
pdf = zeros(3,size(Y,2));
MUMMAE = zeros(3,size(Y,2));

% Overall estimates of IMM filter
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(m,1), size(Y,2));
% Model-conditioned estimates of IMM
MM_i = cell(3,n);
PP_i = cell(3,n);

% Model probabilities 
MU = zeros(3,size(Y,2));

%%%% Prior estimates %%%%

% (MODIFIED) IMM & MMAE
x_ip = [250 100 10 0]';
P_ip = diag([5 5 1 1]);

% (ADDED) MMAE initialization
MMMMAE_i{1,1} = x_ip;
MMMMAE_i{2,1} = x_ip;
MMMMAE_i{3,1} = x_ip;
PPMMAE_i{1,1} = P_ip;
PPMMAE_i{2,1} = P_ip;
PPMMAE_i{3,1} = P_ip;
pdf(1,1) = gauss_pdf(Y(:,1), H{1}*MMMMAE_i{1,1}, H{1}*PPMMAE_i{1,1}*H{1}.'+R{1});
pdf(2,1) = gauss_pdf(Y(:,1), H{2}*MMMMAE_i{2,1}, H{2}*PPMMAE_i{2,1}*H{2}.'+R{2});
pdf(3,1) = gauss_pdf(Y(:,1), H{3}*MMMMAE_i{3,1}, H{3}*PPMMAE_i{3,1}*H{3}.'+R{3});
MUMMAE(:,1) = mu_ip.';

% (ADDED) IMM initialization
MM_i{1,1} = x_ip;
MM_i{2,1} = x_ip;
MM_i{3,1} = x_ip;
PP_i{1,1} = P_ip;
PP_i{2,1} = P_ip;
PP_i{3,1} = P_ip;
MU(:,1) = mu_ip.';
c_j = mu_0j.';

% Filtering steps.
for i = 1:size(Y,2)
    % (ADDED) MMAE
    if i ~= 1
        [MMMMAE_i{1,i}, PPMMAE_i{1,i}] = kf_predict(MMMMAE_i{1,i-1}, PPMMAE_i{1,i-1}, A{1}, Q{1}); % KF_1 prediction
        [MMMMAE_i{2,i}, PPMMAE_i{2,i}] = kf_predict(MMMMAE_i{2,i-1}, PPMMAE_i{2,i-1}, A{2}, Q{2}); % KF_2 prediction
        [MMMMAE_i{3,i}, PPMMAE_i{3,i}] = kf_predict(MMMMAE_i{3,i-1}, PPMMAE_i{3,i-1}, A{3}, Q{3}); % KF_3 prediction
        
        pdf(1,i) = gauss_pdf(Y(:,i),H{1}*MMMMAE_i{1,i},H{1}*PPMMAE_i{1,i}*H{1}.'+R{1})*pdf(1,i-1); % p(Z_i*|alpha_1)
        pdf(2,i) = gauss_pdf(Y(:,i),H{2}*MMMMAE_i{2,i},H{2}*PPMMAE_i{2,i}*H{2}.'+R{2})*pdf(2,i-1); % p(Z_i*|alpha_2)
        pdf(3,i) = gauss_pdf(Y(:,i),H{3}*MMMMAE_i{3,i},H{3}*PPMMAE_i{3,i}*H{3}.'+R{3})*pdf(3,i-1); % p(Z_i*|alpha_3)
        
        % probabilities
        MUMMAE(1,i) = pdf(1,i)*0.3/(pdf(1,i)*0.3+pdf(2,i)*0.3+pdf(3,i)*0.4);
        MUMMAE(2,i) = pdf(2,i)*0.3/(pdf(1,i)*0.3+pdf(2,i)*0.3+pdf(3,i)*0.4);
        MUMMAE(3,i) = pdf(3,i)*0.4/(pdf(1,i)*0.3+pdf(2,i)*0.3+pdf(3,i)*0.4);
    end

    % (ADDED) KF update
    [MMMMAE_i{1,i}, PPMMAE_i{1,i}] = kf_update(MMMMAE_i{1,i}, PPMMAE_i{1,i}, Y(:,i), H{1}, R{1});
    [MMMMAE_i{2,i}, PPMMAE_i{2,i}] = kf_update(MMMMAE_i{2,i}, PPMMAE_i{2,i}, Y(:,i), H{2}, R{2});
    [MMMMAE_i{3,i}, PPMMAE_i{3,i}] = kf_update(MMMMAE_i{3,i}, PPMMAE_i{3,i}, Y(:,i), H{3}, R{3});
    
    % (ADDED) MMAE update
    MMMMAE(:,i) = MUMMAE(1,i)*MMMMAE_i{1,i}(:,1) + MUMMAE(2,i)*MMMMAE_i{2,i}(:,1) + MUMMAE(3,i)*MMMMAE_i{3,i}(:,1);
    
    % IMM - 
    if i ~= 1
        MM_temp{1} = MM_i{1,i-1}; MM_temp{2} = MM_i{2,i-1}; MM_temp{3} = MM_i{3,i-1};
        PP_temp{1} = PP_i{1,i-1}; PP_temp{2} = PP_i{2,i-1}; PP_temp{3} = PP_i{3,i-1};
        [MM_temp, PP_temp, c_j, MM(:,i), PP(:,:,i)] = imm_predict(MM_temp, PP_temp, MU(:,i-1), p_ij, ind, dims, A, Q); % (ADDED) IMM prediction
        MM_i{1,i} = MM_temp{1}; MM_i{2,i} = MM_temp{2}; MM_i{3,i} = MM_temp{3};
        PP_i{1,i} = PP_temp{1}; PP_i{2,i} = PP_temp{2}; PP_i{3,i} = PP_temp{3};
    end
    MM_temp{1} = MM_i{1,i};  MM_temp{2} = MM_i{2,i}; MM_temp{3} = MM_i{3,i};
    PP_temp{1} = PP_i{1,i}; PP_temp{2} = PP_i{2,i}; PP_temp{3} = PP_i{3,i};
    [MM_temp, PP_temp, MU(:,i), MM(:,i), PP(:,:,i)] = imm_update(MM_temp, PP_temp, c_j, ind, dims, Y(:,i), H, R); % (ADDED) IMM update
    MM_i{1,i} = MM_temp{1}; MM_i{2,i} = MM_temp{2}; MM_i{3,i} = MM_temp{3};
    PP_i{1,i} = PP_temp{1}; PP_i{2,i} = PP_temp{2}; PP_i{3,i} = PP_temp{3};
end

%% Calculate the MSEs
% (ADDED)
sumMMAE = zeros(4,1);
sum = zeros(4,1);
% (ADDED) get summation of the squares of the estimation errors from each
% time step
for i = 1:size(Y,2)
    for j = 1:4
        sumMMAE(j,1) = sumMMAE(j,1) + (X_r(j,i) - MMMMAE(j,i))^2;
        sum(j,1) = sum(j,1) + (X_r(j,i) - MM(j,i))^2; 
    end
end
% (ADDED) calculate and print MSE (Mean-squared error)
MSEMMAE = sumMMAE/size(Y,2);
fprintf('MSE(MMAE) - x: %f y: %f v_x: %f v_y: %f \n', MSEMMAE(1), MSEMMAE(2), MSEMMAE(3), MSEMMAE(4));
MSE = sum/size(Y,2);
fprintf('MSE(IMM) - x: %f y: %f v_x: %f v_y: %f \n', MSE(1), MSE(2), MSE(3), MSE(4));

%% (ADDED) Plot estimation result
tspan = dt*(1:n); % (ADDED) time span
figure;
subplot(4,1,1)
plot(tspan, X_r(1,:) - MMMMAE(1,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(1,:) - MM(1,1:n), 'LineWidth', 2)
ylabel('error [m]');
grid on
legend('MMAE','IMM');
title('Position x error');
subplot(4,1,2)
plot(tspan, X_r(2,:) - MMMMAE(2,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(2,:) - MM(2,1:n), 'LineWidth', 2)
grid on
ylabel('error [m]');
title('Position y error');
subplot(4,1,3)
plot(tspan, X_r(3,:) - MMMMAE(3,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(3,:) - MM(3,1:n), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
title('Velocity x error');
subplot(4,1,4)
plot(tspan, X_r(4,:) - MMMMAE(4,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(4,:) - MM(4,1:n), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
xlabel('time [sec]');
title('Velocity y error');

%% (ADDED) Model Probability
figure;
subplot(2,1,1)
plot(tspan, MUMMAE(1,:), 'LineWidth', 2)
hold on
plot(tspan, MUMMAE(2,:), 'LineWidth', 2)
hold on
plot(tspan, MUMMAE(3,:), 'LineWidth', 2)
grid on
ylabel('\mu [-]');
legend('model 1','model 2','model 3');
title('MMAE model probability plot');
subplot(2,1,2)
plot(tspan, MU(1,:), 'LineWidth', 2)
hold on
plot(tspan, MU(2,:), 'LineWidth', 2)
hold on
plot(tspan, MU(3,:), 'LineWidth', 2)
grid on
ylabel('\mu [-]');
xlabel('time [sec]');
title('IMM model probability plot');