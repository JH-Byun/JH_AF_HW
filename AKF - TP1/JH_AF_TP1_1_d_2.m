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
mstate(1:50) = 1;
mstate(51:70) = 1;
mstate(71:120) = 2;
mstate(121:150) = 2;
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

% (ADDED) Overall estimates of MMAE filter
MMMMAE = zeros(size(m,1), size(Y,2));
PPMMAE = zeros(size(m,1), size(m,1), size(Y,2));
% (ADDED) Model-conditioned predictions of MMAE
MMMMAE_i = cell(2,n);
PPMMAE_i = cell(2,n);

% (ADDED) Model probabilities
pdf = zeros(2,size(Y,2));
MUMMAE = zeros(2,size(Y,2));

% (ADDED) Overall estimates of MMAE (modified) filter
MMMMAE2 = zeros(size(m,1), size(Y,2));
PPMMAE2 = zeros(size(m,1), size(m,1), size(Y,2));
% (ADDED) Model-conditioned predictions of MMAE (modified)
MMMMAE2_i = cell(2,n);
PPMMAE2_i = cell(2,n);

% (ADDED) Model probabilities
pdf2 = zeros(2,size(Y,2));
MUMMAE2 = zeros(2,size(Y,2));

%%%% Prior estimates %%%%

% (MODIFIED) MMAE
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
MUMMAE(:,1) = mu_ip.';

% (ADDED) MMAE (modified) initialization
MMMMAE2_i{1,1} = x_ip{1};
MMMMAE2_i{2,1} = x_ip{2};
PPMMAE2_i{1,1} = P_ip{1};
PPMMAE2_i{2,1} = P_ip{2};
pdf2(1,1) = gauss_pdf(Y(:,1), H{1}*x_ip{1}, H{1}*P_ip{1}*H{1}.'+R{1});
pdf2(2,1) = gauss_pdf(Y(:,1), H{2}*x_ip{2}, H{2}*P_ip{2}*H{2}.'+R{2});
MUMMAE2(:,1) = mu_ip.';

% Filtering steps.
for i = 1:size(Y,2)
    % (ADDED) MMAE
    if i ~= 1
        [MMMMAE_i{1,i}, PPMMAE_i{1,i}] = kf_predict(MMMMAE_i{1,i-1}, PPMMAE_i{1,i-1}, A{1}, Q{1}); % KF_1 prediction
        [MMMMAE_i{2,i}, PPMMAE_i{2,i}] = kf_predict(MMMMAE_i{2,i-1}, PPMMAE_i{2,i-1}, A{2}, Q{2}); % KF_2 prediction
          
        pdf(1,i) = gauss_pdf(Y(:,i),H{1}*MMMMAE_i{1,i},H{1}*PPMMAE_i{1,i}*H{1}.'+R{1})*pdf(1,i-1); % p(Z_i*|alpha_1)
        pdf(2,i) = gauss_pdf(Y(:,i),H{2}*MMMMAE_i{2,i},H{2}*PPMMAE_i{2,i}*H{2}.'+R{2})*pdf(2,i-1); % p(Z_i*|alpha_2)
        
        % probabilities
        p_1 = 0.5;
        p_2 = 0.5;
        MUMMAE(1,i) = pdf(1,i)*p_1/(pdf(1,i)*p_1+pdf(2,i)*p_2);
        MUMMAE(2,i) = pdf(2,i)*p_2/(pdf(1,i)*p_1+pdf(2,i)*p_2);
    end

    % (ADDED) KF update
    [MMMMAE_i{1,i}, PPMMAE_i{1,i}] = kf_update(MMMMAE_i{1,i}, PPMMAE_i{1,i}, Y(:,i), H{1}, R{1});
    [MMMMAE_i{2,i}, PPMMAE_i{2,i}] = kf_update(MMMMAE_i{2,i}, PPMMAE_i{2,i}, Y(:,i), H{2}, R{2});
    
    % (ADDED) MMAE update
    MMMMAE(1:4,i) = MUMMAE(1,i)*MMMMAE_i{1,i}(1:4,1) + MUMMAE(2,i)*MMMMAE_i{2,i}(1:4,1);
    MMMMAE(5:6,i) = MMMMAE_i{2,i}(5:6,1);
    
   % (ADDED) MMAE (modified)
    if i ~= 1
        % (ADDED) periodically reinitialization
        if rem(i, 30) == 1
            [MMMMAE2_i{1,i}, PPMMAE2_i{1,i}] = kf_predict(MMMMAE2(1:4,i-1), PPMMAE2_i{1,i-1}, A{1}, Q{1}); % KF_1 prediction
            [MMMMAE2_i{2,i}, PPMMAE2_i{2,i}] = kf_predict(MMMMAE2(:,i-1), PPMMAE2_i{2,i-1}, A{2}, Q{2}); % KF_2 prediction
            pdf2(1,i) = gauss_pdf(Y(:,i),H{1}*MMMMAE2_i{1,i},H{1}*PPMMAE2_i{1,i}*H{1}.'+R{1});
            pdf2(2,i) = gauss_pdf(Y(:,i),H{2}*MMMMAE2_i{2,i},H{2}*PPMMAE2_i{2,i}*H{2}.'+R{2});
        else
            [MMMMAE2_i{1,i}, PPMMAE2_i{1,i}] = kf_predict(MMMMAE2_i{1,i-1}, PPMMAE2_i{1,i-1}, A{1}, Q{1}); % KF_1 prediction
            [MMMMAE2_i{2,i}, PPMMAE2_i{2,i}] = kf_predict(MMMMAE2_i{2,i-1}, PPMMAE2_i{2,i-1}, A{2}, Q{2}); % KF_2 prediction
            pdf2(1,i) = gauss_pdf(Y(:,i),H{1}*MMMMAE2_i{1,i},H{1}*PPMMAE2_i{1,i}*H{1}.'+R{1})*pdf2(1,i-1); % p(Z_i*|alpha_1)
            pdf2(2,i) = gauss_pdf(Y(:,i),H{2}*MMMMAE2_i{2,i},H{2}*PPMMAE2_i{2,i}*H{2}.'+R{2})*pdf2(2,i-1); % p(Z_i*|alpha_2)
        end

%         pdf2(1,i) = gauss_pdf(Y(:,i),H{1}*MMMMAE2_i{1,i},H{1}*PPMMAE2_i{1,i}*H{1}.'+R{1})*pdf2(1,i-1); % p(Z_i*|alpha_1)
%         pdf2(2,i) = gauss_pdf(Y(:,i),H{2}*MMMMAE2_i{2,i},H{2}*PPMMAE2_i{2,i}*H{2}.'+R{2})*pdf2(2,i-1); % p(Z_i*|alpha_2)
        MUMMAE2(1,i) = pdf2(1,i)*p_1/(pdf2(1,i)*p_1+pdf2(2,i)*p_2);
        MUMMAE2(2,i) = 1 - MUMMAE2(1,i);
    end

    % (ADDED) KF update
    [MMMMAE2_i{1,i}, PPMMAE2_i{1,i}] = kf_update(MMMMAE2_i{1,i}, PPMMAE2_i{1,i}, Y(:,i), H{1}, R{1});
    [MMMMAE2_i{2,i}, PPMMAE2_i{2,i}] = kf_update(MMMMAE2_i{2,i}, PPMMAE2_i{2,i}, Y(:,i), H{2}, R{2});
    
    % (ADDED) MMAE update
    MMMMAE2(1:4,i) = MUMMAE2(1,i)*MMMMAE2_i{1,i}(1:4,1) + MUMMAE2(2,i)*MMMMAE2_i{2,i}(1:4,1);
    MMMMAE2(5:6,i) = MMMMAE2_i{2,i}(5:6,1);
end

%% Calculate the MSEs
% (ADDED)
sumMMAE = zeros(4,1);
sumMMAE2 = zeros(4,1);
% (ADDED) get summation of the squares of the estimation errors from each
% time step
for i = 1:size(Y,2)
    for j = 1:4
        sumMMAE(j,1) = sumMMAE(j,1) + (X_r(j,i) - MMMMAE(j,i))^2;
        sumMMAE2(j,1) = sumMMAE2(j,1) + (X_r(j,i) - MMMMAE2(j,i))^2;
    end
end
% (ADDED) calculate and print MSE (Mean-squared error)
MSEMMAE = sumMMAE/size(Y,2);
fprintf('MSE(MMAE) - x: %f y: %f v_x: %f v_y: %f \n', MSEMMAE(1), MSEMMAE(2), MSEMMAE(3), MSEMMAE(4));
MSEMMAE2 = sumMMAE2/size(Y,2);
fprintf('MSE(MMAE (modified)) - x: %f y: %f v_x: %f v_y: %f \n', MSEMMAE2(1), MSEMMAE2(2), MSEMMAE2(3), MSEMMAE2(4));

%% (ADDED) Plot estimation result
tspan = dt*(1:n); % (ADDED) time span
figure;
subplot(4,1,1)
plot(tspan, X_r(1,:) - MMMMAE(1,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(1,:) - MMMMAE2(1,1:n), 'LineWidth', 2)
ylabel('error [m]');
grid on
legend('MMAE','MMAE (modified)');
title('Position x error');
subplot(4,1,2)
plot(tspan, X_r(2,:) - MMMMAE(2,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(2,:) - MMMMAE2(2,1:n), 'LineWidth', 2)
grid on
ylabel('error [m]');
title('Position y error');
subplot(4,1,3)
plot(tspan, X_r(3,:) - MMMMAE(3,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(3,:) - MMMMAE2(3,1:n), 'LineWidth', 2)
grid on
ylabel('error [m/2]');
title('Velocity x error');
subplot(4,1,4)
plot(tspan, X_r(4,:) - MMMMAE(4,1:n), 'LineWidth', 2)
hold on
plot(tspan, X_r(4,:) - MMMMAE2(4,1:n), 'LineWidth', 2)
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
grid on
ylabel('\mu [-]');
legend('model 1','model 2');
title('MMAE model probability plot');
subplot(2,1,2)
plot(tspan, MUMMAE2(1,:), 'LineWidth', 2)
hold on
plot(tspan, MUMMAE2(2,:), 'LineWidth', 2)
grid on
ylabel('\mu [-]');
xlabel('time [sec]');
title('MMAE (modified) model probability plot');