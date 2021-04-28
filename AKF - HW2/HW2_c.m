clear all; close all; clc
%% HW - 2_c (JH version)

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
t_f = 300;
tspan = 0:T:t_f;
theta = 0.0005;
x = zeros(4, length(tspan));
x(:,1) = x_0;
y = zeros(2, length(tspan));
v_0 = mvnrnd(zeros(2,1), R).';
y(:,1) = H*x_0 + v_0;

% Kalman filter (KF)
MM_KF = zeros(4, length(tspan));
PP_KF = zeros(4,4,length(tspan));
MM_KF(:,1) = xhat_0;
PP_KF(:,:,1) = P_0;

% H_inf filter
MM_H_inf = zeros(4, length(tspan));
PP_H_inf = zeros(4,4,length(tspan));
MM_H_inf(:,1) = xhat_0;
PP_H_inf(:,:,1) = P_0;

for i = 1:length(tspan)-1
    % simulation
    w_i = mvnrnd(zeros(4,1), Q).';
    x(:,i+1) = F*x(:,i) + G*u + w_i;
    v_i = mvnrnd(zeros(2,1), R).';
    y(:,i+1) = H*x(:,i+1) + v_i;
    
    % kalman filter
        % update
        K_KF = PP_KF(:,:,i)*H.'*inv(H*PP_KF(:,:,i)*H.'+R);
        MM_KF(:,i) = MM_KF(:,i) + K_KF*(y(:,i)-H*MM_KF(:,i));
        PP_KF(:,:,i) = (eye(4) - K_KF*H)*PP_KF(:,:,i);
    
        % prediction
        PP_KF(:,:,i+1) = F*PP_KF(:,:,i)*F.' + Q;
        MM_KF(:,i+1) = F*MM_KF(:,i) + G*u;
        
        % loop eigenvalues
        eig_KF(:,i) = eig((eye(4)-K_KF*H)*F);
        
    % H_inf filter
        % filter
        K_H_inf = PP_H_inf(:,:,i)*inv(eye(4)-theta*Sbar*PP_H_inf(:,:,i)+H.'*inv(R)*H*PP_H_inf(:,:,i))*H.'*inv(R);
        MM_H_inf(:,i+1) = F*MM_H_inf(:,i)+ G*u + F*K_H_inf*(y(:,i)-H*MM_H_inf(:,i));
        PP_H_inf(:,:,i+1) = F*PP_H_inf(:,:,i)*inv(eye(4)-theta*Sbar*PP_H_inf(:,:,i)+H.'*inv(R)*H*PP_H_inf(:,:,i))*F.' + Q;
        
        % loop eigenvalues
        eig_H_inf(:,i) = eig((eye(4)-K_H_inf*H)*F);
end

%% plot eigenvalues
r=1;
x0=0;
y0=0;
theta = linspace(0,2*pi,100);

eig_KF_real = real(eig_KF);
eig_KF_imag = imag(eig_KF);
figure;
plot(x0 + r*cos(theta),y0 + r*sin(theta),'--')
hold on
plot(eig_KF_real(1,:), eig_KF_imag(1,:))
hold on
plot(eig_KF_real(2,:), eig_KF_imag(2,:))
hold on
plot(eig_KF_real(3,:), eig_KF_imag(3,:))
hold on
plot(eig_KF_real(4,:), eig_KF_imag(4,:))
grid on
xlim([-2 2]);
ylim([-2 2]);
xlabel('Real axis');
ylabel('Imag. axis');
title('KF closed-loop poles');

eig_H_inf_real = real(eig_H_inf);
eig_H_inf_imag = imag(eig_H_inf);
figure;
plot(x0 + r*cos(theta),y0 + r*sin(theta),'--')
hold on
plot(eig_H_inf_real(1,:), eig_H_inf_imag(1,:))
hold on
plot(eig_H_inf_real(2,:), eig_H_inf_imag(2,:))
hold on
plot(eig_H_inf_real(3,:), eig_H_inf_imag(3,:))
hold on
plot(eig_H_inf_real(4,:), eig_H_inf_imag(4,:))
grid on
xlim([-2 2]);
ylim([-2 2]);
xlabel('Real axis');
ylabel('Imag. axis');
title('H_{inf} closed-loop poles');