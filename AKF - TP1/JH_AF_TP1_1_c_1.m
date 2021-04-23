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

mu_ip = [0.1 0.9];
mu_0j = mu_ip;

% case 1 : ideal
p_ij = [0.98 0.02;
        0.02 0.98];

% (MODIFIED) Forced mode transitions 
mstate(1:50) = 2;
mstate(51:70) = 1;
mstate(71:120) = 2;
mstate(121:150) = 1;
mstate(151:200) = 2;

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

% KF using model 1 
MM1 = zeros(size(A{1},1), size(Y,2));
PP1 = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% KF using model 2
MM2 = zeros(size(A{2},1), size(Y,2));
PP2 = zeros(size(A{2},1), size(A{2},1), size(Y,2));

% (ADDED) Overall estimates of MMAE filter
MMMMAE = zeros(size(m,1), size(Y,2));
PPMMAE = zeros(size(m,1), size(m,1), size(Y,2));
% (ADDED) Model-conditioned predictions of MMAE
MMMMAE_i = cell(2,n);
PPMMAE_i = cell(2,n);

% (ADDED) Model probabilities
pdf = zeros(2,size(Y,2));
MUMMAE = zeros(2,size(Y,2));

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
% (ADDED) initialization
MM1(:,1) = M1;
PP1(:,:,1) = P1;

% KF2
M2 = [0 0 0 -1 0 0]';
P2 = diag([0.1 0.1 0.1 0.1 0.5 0.5]);
% (ADDED) initialization
MM2(:,1) = M2;
PP2(:,:,1) = P2;

% (MODIFIED) IMM & MMAE
x_ip{1} = [0 0 0 -1]';
P_ip{1} = diag([0.1^2 0.1^2 0.1^2 0.1^2]);
x_ip{2} = [0 0 0 -1 0 0]';
P_ip{2} = diag([0.1^2 0.1^2 0.1^2 0.1^2 0.5^2 0.5^2]);

% (ADDED) MMAE initialization
MMMMAE_i{1,1} = x_ip{1};
MMMMAE_i{2,1} = x_ip{2};
PPMMAE_i{1,1} = P_ip{1};
PPMMAE_i{2,1} = P_ip{2};
pdf(1,1) = gauss_pdf(Y(:,1), H{1}*x_ip{1}, H{1}*P_ip{1}*H{1}.'+R{1});
pdf(2,1) = gauss_pdf(Y(:,1), H{2}*x_ip{2}, H{2}*P_ip{2}*H{2}.'+R{2});
MUMMAE = mu_ip.';

% (ADDED) IMM initialization
MM_i{1,1} = x_ip{1};
MM_i{2,1} = x_ip{2};
PP_i{1,1} = P_ip{1};
PP_i{2,1} = P_ip{2};
MU(:,1) = mu_ip.';
c_j = mu_0j.';

% Filtering steps.
for i = 1:size(Y,2)
    % KF with model 1
    if i ~= 1
        [MM1(:,i), PP1(:,:,i)] = kf_predict(MM1(:,i-1), PP1(:,:,i-1), A{1}, Q{1}); % (ADDED) KF_1 prediction
    end
    [MM1(:,i), PP1(:,:,i)] = kf_update(MM1(:,i), PP1(:,:,i), Y(:,i), H{1}, R{1}); % (ADDED) KF_1 update
    
    % KF with model 2
    if i ~= 1
        [MM2(:,i), PP2(:,:,i)] = kf_predict(MM2(:,i-1), PP2(:,:,i-1), A{2}, Q{2}); % (ADDED) KF_2 prediction
    end
    [MM2(:,i), PP2(:,:,i)] = kf_update(MM2(:,i), PP2(:,:,i), Y(:,i), H{2}, R{2}); % (ADDED) KF_2 update

    % (ADDED) MMAE
    if i ~= 1
        [MMMMAE_i{1,i}, PPMMAE_i{1,i}] = kf_predict(MMMMAE_i{1,i-1}, PPMMAE_i{1,i-1}, A{1}, Q{1}); % KF_1 prediction
        [MMMMAE_i{2,i}, PPMMAE_i{2,i}] = kf_predict(MMMMAE_i{2,i-1}, PPMMAE_i{2,i-1}, A{2}, Q{2}); % KF_2 prediction
        
        pdf(1,i) = gauss_pdf(Y(:,i),H{1}*MMMMAE_i{1,i},H{1}*PPMMAE_i{1,i}*H{1}.'+R{1})*pdf(1,i-1); % p(Z_i*|alpha_1)
        pdf(2,i) = gauss_pdf(Y(:,i),H{2}*MMMMAE_i{2,i},H{2}*PPMMAE_i{2,i}*H{2}.'+R{2})*pdf(2,i-1); % p(Z_i*|alpha_2)

        % probabilities
        MUMMAE(1,i) = pdf(1,i)*0.5/(pdf(1,i)*0.5+pdf(2,i)*0.5);
        MUMMAE(2,i) = pdf(2,i)*0.5/(pdf(1,i)*0.5+pdf(2,i)*0.5); 
    end

    % (ADDED) KF update
    [MMMMAE_i{1,i}, PPMMAE_i{1,i}] = kf_update(MMMMAE_i{1,i}, PPMMAE_i{1,i}, Y(:,i), H{1}, R{1});
    [MMMMAE_i{2,i}, PPMMAE_i{2,i}] = kf_update(MMMMAE_i{2,i}, PPMMAE_i{2,i}, Y(:,i), H{2}, R{2});
    
    % (ADDED) MMAE update
    MMMMAE(1:4,i) = MUMMAE(1,i)*MMMMAE_i{1,i}(1:4,1) + MUMMAE(2,i)*MMMMAE_i{2,i}(1:4,1);
    MMMMAE(5:6,i) = MMMMAE_i{2,i}(5:6,1);
    
    % IMM
    if i ~= 1
        MM_temp2{1} = MM_i{1,i-1};  MM_temp2{2} = MM_i{2,i-1};
        PP_temp2{1} = PP_i{1,i-1}; PP_temp2{2} = PP_i{2,i-1};
        [MM_temp2, PP_temp2, c_j, MM(:,i), PP(:,:,i)] = imm_predict(MM_temp2, PP_temp2, MU(:,i-1), p_ij, ind, dims, A, Q); % (ADDED) IMM prediction
        MM_i{1,i} = MM_temp2{1}; MM_i{2,i} = MM_temp2{2};
        PP_i{1,i} = PP_temp2{1}; PP_i{2,i} = PP_temp2{2};
    end
    MM_temp2{1} = MM_i{1,i};  MM_temp2{2} = MM_i{2,i};
    PP_temp2{1} = PP_i{1,i}; PP_temp2{2} = PP_i{2,i};
    [MM_temp2, PP_temp2, MU(:,i), MM(:,i), PP(:,:,i)] = imm_update(MM_temp2, PP_temp2, c_j, ind, dims, Y(:,i), H, R); % (ADDED) IMM update
    MM_i{1,i} = MM_temp2{1}; MM_i{2,i} = MM_temp2{2};
    PP_i{1,i} = PP_temp2{1}; PP_i{2,i} = PP_temp2{2};
end

%% Calculate the MSEs
% (ADDED)
sum1 = zeros(4,1);
sum2 = zeros(4,1);
sumMMAE = zeros(4,1);
sum = zeros(4,1);
% (ADDED) get summation of the squares of the estimation errors from each
% time step
for i = 1:size(Y,2)
    for j = 1:4
        sum1(j,1) = sum1(j,1) + (X_r(j,i) - MM1(j,i))^2;
        sum2(j,1) = sum2(j,1) + (X_r(j,i) - MM2(j,i))^2;
        sumMMAE(j,1) = sumMMAE(j,1) + (X_r(j,i) - MMMMAE(j,i))^2;
        sum(j,1) = sum(j,1) + (X_r(j,i) - MM(j,i))^2; 
    end
end
% (ADDED) calculate and print MSE (Mean-squared error)
MSE1 = sum1/size(Y,2);
fprintf('MSE(KF_1) - x: %f y: %f v_x: %f v_y: %f \n', MSE1(1), MSE1(2), MSE1(3), MSE1(4));
MSE2 = sum2/size(Y,2);
fprintf('MSE(KF_2) - x: %f y: %f v_x: %f v_y: %f \n', MSE2(1), MSE2(2), MSE2(3), MSE2(4));
MSEMMAE = sumMMAE/size(Y,2);
fprintf('MSE(MMAE) - x: %f y: %f v_x: %f v_y: %f \n', MSEMMAE(1), MSEMMAE(2), MSEMMAE(3), MSEMMAE(4));
MSE = sum/size(Y,2);
fprintf('MSE(IMM) - x: %f y: %f v_x: %f v_y: %f \n', MSE(1), MSE(2), MSE(3), MSE(4));

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