clc; clear all; close all;

load 2_d_UKFPF.mat

e_x_d = X - xhat;
% ESS_d = ESS;
time_spent_d = time_spent;



MC = size(e_x_d,1);
simLen = size(e_x_d,2);

for i = 1:simLen
    RMS_d(:,i) = sqrt(sum(e_x_d(:,i).^2)/MC);
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% UKFPF
figure('DefaultAxesFontSize',default_font_size);
plot(1:simLen,RMS_d,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF + PF: RMSE (N = 500)')
%%
subplot(2,1,2)
plot(1:simLen,ESS_100,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: ESS (N = 500)')

% UKF
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_500,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: RMSE (N = 500)')
subplot(2,1,2)
plot(1:simLen,ESS_500,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: ESS (N = 500)')

% SIS
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_1000,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: RMSE (N = 1000)')
subplot(2,1,2)
plot(1:simLen,ESS_1000,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: ESS (N = 1000)')

% computation time comparison
figure('DefaultAxesFontSize',default_font_size);
plot(1:MC, time_spent_100,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_500,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_1000,'LineWidth',default_line_width)
grid on
xlabel('MC iteration step [-]')
title('Computation time comparison: particle numbers N')
legend({'N = 100','N = 500','N = 1000'},'Location', 'best')