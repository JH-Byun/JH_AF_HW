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


% Filtering steps.
for i = 1:size(Y,2)
    % KF with model 1


    % KF with model 2

    
    % IMM

end


% Calculate the MSEs

