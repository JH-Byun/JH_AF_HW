clear all; close all; clc
%% HW - 2_b (JH version)

%% parameters
T = 3;
u = 2;
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
        MM_KF(:,i+1) = F*MM_KF(:,i) + G*1;
        
    % H_inf filter
    K_H_inf = PP_H_inf(:,:,i)*inv(eye(4)-theta*Sbar*PP_H_inf(:,:,i)+H.'*inv(R)*H*PP_H_inf(:,:,i))*H.'*inv(R);
    MM_H_inf(:,i+1) = F*MM_H_inf(:,i)+ G*1 + F*K_H_inf*(y(:,i)-H*MM_H_inf(:,i));
    PP_H_inf(:,:,i+1) = F*PP_H_inf(:,:,i)*inv(eye(4)-theta*Sbar*PP_H_inf(:,:,i)+H.'*inv(R)*H*PP_H_inf(:,:,i))*F.' + Q;
end

figure;
subplot(2,2,1)
plot(tspan, x(1,:) - MM_KF(1,:), 'LineWidth', 2);
hold on
plot(tspan, x(1,:) - MM_H_inf(1,:), 'LineWidth', 2);
grid on
title('north position');
ylabel('x [m]');
legend('KF','H_{inf}');
subplot(2,2,2)
plot(tspan, x(2,:) - MM_KF(2,:), 'LineWidth', 2);
hold on
plot(tspan, x(2,:) - MM_H_inf(2,:), 'LineWidth', 2);
grid on
title('east position');
ylabel('y [m]');
legend('KF','H_{inf}');
subplot(2,2,3)
plot(tspan, x(3,:) - MM_KF(3,:), 'LineWidth', 2);
hold on
plot(tspan, x(3,:) - MM_H_inf(3,:), 'LineWidth', 2);
grid on
title('north velocity');
ylabel('v_x [m/s]');
legend('KF','H_{inf}');
subplot(2,2,4)
plot(tspan, x(4,:) - MM_KF(4,:), 'LineWidth', 2);
hold on
plot(tspan, x(4,:) - MM_H_inf(4,:), 'LineWidth', 2);
grid on
title('east velocity');
ylabel('v_y [m/s]');
legend('KF','H_{inf}');

%% calculate RMS
sum_KF = zeros(4,1);
sum_H_inf = zeros(4,1);

for i = 1:length(tspan)
    sum_KF = sum_KF + (x(:,i) - MM_KF(:,i)).^2;
    sum_H_inf = sum_H_inf + (x(:,i) - MM_H_inf(:,i)).^2;
end

for j = 1:4
    RMS_KF(j,1) = sqrt(sum_KF(j,1)/size(x,2));
    RMS_H_inf(j,1) = sqrt(sum_H_inf(j,1)/size(x,2));
end

fprintf('x(RMS): %f (KF), %f (H_{inf}) \n', RMS_KF(1,1), RMS_H_inf(1,1));
fprintf('y(RMS): %f (KF), %f (H_{inf}) \n', RMS_KF(2,1), RMS_H_inf(2,1));
fprintf('v_x(RMS): %f (KF), %f (H_{inf}) \n', RMS_KF(3,1), RMS_H_inf(3,1));
fprintf('v_y(RMS): %f (KF), %f (H_{inf}) \n', RMS_KF(4,1), RMS_H_inf(4,1));