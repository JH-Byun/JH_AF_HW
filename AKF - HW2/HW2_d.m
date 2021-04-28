clear all; close all; clc
%% HW - 2_d (JH version)

%% parameters
T = 3;
u = 1;
Q = diag([4 4 1 1]);
R = diag([900 900]);
h_a = 0.9*pi; % heading angle
x_0 = zeros(4,1);
xhat_0 = x_0;
P_0 = eye(4);

% model
F = [1 0 T 0;...
    0 1 0 T;...
    0 0 1 0;...
    0 0 0 1];
G = [0;0;T*sin(h_a);T*cos(h_a)];
H = [eye(2) zeros(2,2)];

%% simulation
S = eye(4);
L = eye(4);
Sbar = L.'*S*L;
theta_m = 1;
thetaspan = 0:0.0001:theta_m;

for i = 1:length(thetaspan)
    sol = idare(F.',eye(4),Q,H.'*inv(R)*H-thetaspan(i)*Sbar,zeros(4,4),eye(4));
end