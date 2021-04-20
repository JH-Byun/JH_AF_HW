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
p_ij = [0.98 0.02;
        0.02 0.98];

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
P = diag([0.3 0.3 0.3 0.3 0.7 0.7]);

%%%% Space for the estimates %%%% 

% KF using model 1 
MM1 = zeros(size(A{1},1), size(Y,2));
PP1 = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% KF using model 2
MM2 = zeros(size(A{2},1), size(Y,2));
PP2 = zeros(size(A{2},1), size(A{2},1), size(Y,2));

% Overall estimates of IMM filter
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(m,1), size(Y,2));
% Model-conditioned estimates of IMM
MM_i = cell(2,n);
PP_i = cell(2,n);

% Model probabilities 
MU = zeros(2,size(Y,2));

%%%% Prior estimates %%%%

% KF1
M1 = [0 0 0 -1]';
P1 = diag([0.1 0.1 0.1 0.1]);

% KF2
M2 = [0 0 0 -1 0 0]';
P2 = diag([0.1 0.1 0.1 0.1 0.5 0.5]);

% IMM
x_ip{1} = [0 0 0 -1]';
P_ip{1} = diag([0.1 0.1 0.1 0.1]);
x_ip{2} = [0 0 0 -1 0 0]';
P_ip{2} = diag([0.1 0.1 0.1 0.1 0.5 0.5]);

% (ADDED) estimated value storage
X_KF_1 = zeros(4,200);
X_KF_2 = zeros(6,200);
X_MMAE = zeros(6,200);
X_IMM = zeros(4,200);

% (ADDED) state covariance storage
P_KF_1 = zeros(4,4,200);
P_KF_2 = zeros(6,6,200);
P_MMAE = zeros(6,6,200);
P_IMM = zeros(6,6,200);

% (ADDED) initial estimation
X_KF_1(:,1) = M1;
P_KF_1(:,:,1) = P1;
X_KF_2(:,1) = M2;
P_KF_2(:,:,1) = P2;

% Filtering steps.
for i = 1:size(Y,2)
    % KF with model 1
    [X_KF_1(:,i), P_KF_1(:,:,i)] = kf_predict(X_KF_1(:,i), P_KF_1(:,:,i), A{1}, Q{1}); % (ADDED) KF_1 prediction
    [X_KF_1(:,i+1), P_KF_1(:,:,i+1)] = kf_update(X_KF_1(:,i), P_KF_1(:,:,i), Y(:,i), H{1}, R{1}); % (ADDED) KF_1 update
    
    % KF with model 2
    [X_KF_2(:,i), P_KF_2(:,:,i)] = kf_predict(X_KF_2(:,i), P_KF_2(:,:,i), A{2}, Q{2}); % (ADDED) KF_2 prediction
    [X_KF_2(:,i+1), P_KF_2(:,:,i+1)] = kf_update(X_KF_2(:,i), P_KF_2(:,:,i), Y(:,i), H{2}, R{2}); % (ADDED) KF_2 update

    % MMAE


    % IMM
    
end


% Calculate the MSEs

% (ADDED) Plot results
tspan = dt*(1:200); % (ADDED) time span
figure;
subplot(4,1,1)
plot(tspan, X_r(1,:) - X_KF_1(1,1:200), 'LineWidth', 2)
hold on
plot(tspan, X_r(1,:) - X_KF_2(1,1:200), 'LineWidth', 2)
ylabel('error [m]');
grid on
legend('KF1','KF2','MMAE','IMM');
title('Position x error');
subplot(4,1,2)
plot(tspan, X_r(2,:) - X_KF_1(2,1:200), 'LineWidth', 2)
hold on
plot(tspan, X_r(2,:) - X_KF_2(2,1:200), 'LineWidth', 2)
grid on
ylabel('error [m]');
legend('KF1','KF2','MMAE','IMM');
title('Position y error');
subplot(4,1,3)
plot(tspan, X_r(3,:) - X_KF_1(3,1:200), 'LineWidth', 2)
hold on
plot(tspan, X_r(3,:) - X_KF_2(3,1:200), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
legend('KF1','KF2','MMAE','IMM');
title('Velocity x error');
subplot(4,1,4)
plot(tspan, X_r(4,:) - X_KF_1(4,1:200), 'LineWidth', 2)
hold on
plot(tspan, X_r(4,:) - X_KF_2(4,1:200), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
xlabel('time [sec]');
legend('KF1','KF2','MMAE','IMM');
title('Velocity y error');